import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pathlib import Path
from typing import Optional
import h5py
from ios import (
    lsm_metadata,
    read_image,
    read_mask,
    find_file_paths,
    find_sample_dirs,
    write_dict_to_hdf_group,
    write_dataframe_to_hdf_group,
)
from fitting import recovery_fit
from calculator import phair_double_normalization, subsample
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
import time


def log_message(message: str, log_dict: dict, sample_name: str):
    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S.%f", time.localtime())
    log_dict[timestamp] = f"{sample_name} - {message}"


class FrapSample:
    def __init__(self, root_dir: Path):
        self.logs = {}
        self.root = root_dir

        self.image_paths: dict = find_file_paths(self.root)

        self.metadata = lsm_metadata(self.image_paths["lsm_path"])
        self.metadata.update({"sample_name": self.root.name, "sample_path": self.root})

        # Prep vars for images and masks for later population
        self.image_sequence: NDArray = np.array([])
        self.masks: dict[str, NDArray] = {
            "irradiated": np.array([]),
            "reference": np.array([]),
            "background": np.array([]),
        }

        self.data = {}

        self.populate_images_and_masks()

    def populate_images_and_masks(self):

        self.image_sequence = read_image(self.image_paths["image_path"])

        if not self.image_sequence.ndim == 3:
            raise ValueError(
                f"Expected image sequence to have 3 dimensions (time, x, y), but got {self.image_sequence.shape}"
            )
        if not self.image_sequence.shape[0] == self.metadata["time_size"]:
            raise ValueError(
                f"Images does not have the same time shape.{self.image_sequence.shape[0]}"
            )

        self.masks["irradiated"] = read_mask(self.image_paths["irradiated_mask_path"])
        log_message("Irradiated mask loaded successfully.", self.logs, self.root.name)
        self.masks["reference"] = read_mask(self.image_paths["correction_mask_path"])
        log_message("Reference mask loaded successfully.", self.logs, self.root.name)
        self.masks["background"] = read_mask(self.image_paths["background_mask_path"])
        log_message("Background mask loaded successfully.", self.logs, self.root.name)

        self.process_data()

    def process_data(self):

        # Initial check for corect image population
        if self.image_sequence is None:
            raise ValueError("Image sequence has not been populated.")
        if any(mask.size == 0 for mask in self.masks.values()):
            raise ValueError("One or more masks have not been populated.")

        normalized_intensity_dict = phair_double_normalization(
            image_sequence=self.image_sequence,
            irradiated_mask=self.masks["irradiated"],
            reference_mask=self.masks["reference"],
            background_mask=self.masks["background"],
            pre_bleach_frames=self.metadata["pre_bleach_frame_count"],
        )

        log_message("Data normalized", self.logs, self.root.name)

        if len(self.metadata["time"]) != len(
            normalized_intensity_dict["intensity_normalized"]
        ):
            raise ValueError(
                f"Time points and normalized intensity data do not have the same length: {len(self.metadata['time'])} vs {len(self.data['normalized_intensity'])}"
            )

        self.data.update(
            {
                "time": self.metadata["time"],
                "index": np.arange(len(self.metadata["time"]))
                - self.metadata["pre_bleach_frame_count"],
                "normalized_intensity": normalized_intensity_dict[
                    "intensity_normalized"
                ],
                "all_intensities": normalized_intensity_dict,
            }
        )

        self.data["fit_params"], self.data["fit_values"] = recovery_fit(
            time=self.data["time"],
            intensity=self.data["normalized_intensity"],
            model="one_phase",
            low_limit_index=self.metadata["pre_bleach_frame_count"],
        )

        self.fit_type = "one_phase"
        log_message(
            f"{self.fit_type} recovery fit completed", self.logs, self.root.name
        )

    def generate_dict_data(self):
        dict_data = []

        # Iterate though time points and gather all relevant data into a list of dictionaries
        for idx, time_point in enumerate(self.data["time"]):
            dict_data.append(
                {
                    "sample_name": self.metadata["sample_name"],
                    "time": time_point,
                    "index": self.data["index"][idx],
                    "normalized_intensity": self.data["normalized_intensity"][idx],
                    "intensity_irradiated_raw": self.data["all_intensities"][
                        "intensity_irradiated_raw"
                    ][idx],
                    "intensity_reference_raw": self.data["all_intensities"][
                        "intensity_reference_raw"
                    ][idx],
                    "intensity_background_raw": self.data["all_intensities"][
                        "intensity_background_raw"
                    ][idx],
                    "time_delta": self.metadata["mean_timestamp_difference"],
                }
            )
        return dict_data

    def generate_report(self, axes=None):

        required_subplots = 5
        if axes is None:
            fig, axes = plt.subplots(1, required_subplots, figsize=(16, 6))
            show_plot = True
        else:
            show_plot = False

        if len(axes) != required_subplots:
            raise ValueError(
                f"Expected axes to have shape {required_subplots}, but got {axes.shape}"
            )

        sns.lineplot(
            x=self.data["time"],
            y=self.data["normalized_intensity"],
            ax=axes[0],
            markers=True,
            dashes=True,
            # palette="gray",
        )
        sns.lineplot(
            x=self.data["fit_values"]["time"],
            y=self.data["fit_values"]["fitted_intensity"],
            ax=axes[0],
            label=f"Fit ({self.fit_type})",
            # palette="red",
        )
        axes[0].set_title(
            f"{self.root.name} Recovery Curve (tau = {self.data['fit_params']['tau']:.2f} s)"
        )
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Normalized Intensity")
        axes[0].set_ylim(0, 1.2)

        sns.lineplot(
            x=self.data["time"],
            y=self.data["all_intensities"]["intensity_irradiated_raw"],
            ax=axes[1],
            label="Irradiated",
        )
        sns.lineplot(
            x=self.data["time"],
            y=self.data["all_intensities"]["intensity_reference_raw"],
            ax=axes[1],
            label="Reference",
        )
        sns.lineplot(
            x=self.data["time"],
            y=self.data["all_intensities"]["intensity_background_raw"],
            ax=axes[1],
            label="Background",
        )
        axes[1].set_title("Raw Intensity Curves")

        pre_bleach_frames = self.metadata["pre_bleach_frame_count"]
        pre_bleach_vmax = np.percentile(
            self.image_sequence[: pre_bleach_frames - 2], 99.5
        )
        axes[2].imshow(
            self.image_sequence[pre_bleach_frames - 1],
            cmap="gray",
            vmax=pre_bleach_vmax,
        )
        axes[2].contour(self.masks["irradiated"], colors="royalblue", linewidths=0.5)
        axes[2].contour(self.masks["reference"], colors="orange", linewidths=0.5)
        axes[2].contour(self.masks["background"], colors="green", linewidths=0.5)
        axes[2].set_title("Pre-bleach Image")

        axes[3].imshow(
            self.image_sequence[pre_bleach_frames],
            cmap="gray",
            vmax=pre_bleach_vmax,
        )
        axes[3].contour(self.masks["irradiated"], colors="r", linewidths=0.5)
        axes[3].set_title("Post-bleach Image")

        axes[4].imshow(self.image_sequence[-1], cmap="gray", vmax=pre_bleach_vmax)
        axes[4].contour(self.masks["irradiated"], colors="r", linewidths=0.5)
        axes[4].set_title("Final Image")

        for i in [2, 3, 4]:
            axes[i].axis("off")

        if show_plot:
            sns.despine()
            plt.tight_layout()
            plt.show()

    def export_to_hdf(self, hdf_path: Path, experiment_name: str):
        hdf_path = Path(hdf_path)

        with h5py.File(hdf_path, "a") as hdf_file:
            samples_group = hdf_file.require_group(experiment_name).require_group(
                "FrapSamples"
            )

            if self.metadata["sample_name"] in samples_group:
                del samples_group[self.metadata["sample_name"]]
            sample_group = samples_group.create_group(self.metadata["sample_name"])

            images_group = sample_group.create_group("images")
            images_group.create_dataset("image_sequence", data=self.image_sequence)
            images_group.create_dataset(
                "mask_irradiated", data=self.masks["irradiated"].astype(np.uint8)
            )
            images_group.create_dataset(
                "mask_background", data=self.masks["background"].astype(np.uint8)
            )
            images_group.create_dataset(
                "mask_reference", data=self.masks["reference"].astype(np.uint8)
            )

            metadata_group = sample_group.create_group("meta_data")
            write_dict_to_hdf_group(metadata_group, self.metadata)

            intensity_data_group = sample_group.create_group("intensity_data")
            intensity_data = {
                "time": self.data["time"],
                "index": self.data["index"],
                "normalized_intensity": self.data["normalized_intensity"],
                "intensity_irradiated_raw": self.data["all_intensities"][
                    "intensity_irradiated_raw"
                ],
                "intensity_reference_raw": self.data["all_intensities"][
                    "intensity_reference_raw"
                ],
                "intensity_background_raw": self.data["all_intensities"][
                    "intensity_background_raw"
                ],
            }
            write_dict_to_hdf_group(intensity_data_group, intensity_data)

            fit_data_group = sample_group.create_group("fit_data")
            fit_data = {
                "fit_type": self.fit_type,
                "fit_params": self.data["fit_params"],
                "fit_values": self.data["fit_values"],
            }
            write_dict_to_hdf_group(fit_data_group, fit_data)


class FrapExperiment:
    def __init__(self, root_dir: Path, name=None):
        self.root = root_dir
        self.name = root_dir.name if name is None else name
        self.logs = {}
        self.sample_paths = find_sample_dirs(root_dir)
        self.samples = []
        self.metadata = {
            "experiment_name": self.root.name,
            "experiment_path": self.root,
        }
        self.fit_params = {}
        self.fit_values = {}

        self.df_averaged = pd.DataFrame()
        self.timepoint_dicts = []

        for sample_path in self.sample_paths:
            log_message("Loading sample...", self.logs, sample_path.name)
            # try:
            self.samples.append(FrapSample(sample_path))
            self.logs.update({k: v for k, v in self.samples[-1].logs.items()})
            log_message("Sample loaded successfully.", self.logs, sample_path.name)
            # # except Exception as error:
            #     log_message(
            #         f"Failed to load sample: {error}", self.logs, sample_path.name
            #     )
            #
        for ts, msg in self.logs.items():
            print(f"{ts} - {msg}")

        self.process_sample_data()

        self.overview_figure_samples()

    def get_dataframes(self):
        df_meta = pd.DataFrame([fs.metadata for fs in self.samples])

        df_fit = pd.DataFrame([fs.data["fit_params"] for fs in self.samples])

        df_intensity = pd.DataFrame([fs.data["all_intensities"] for fs in self.samples])

        return df_meta, df_fit, df_intensity

    def generate_experiment_df(self):
        """ """
        all_data = []
        for sample in self.samples:
            sample_dict_data = sample.generate_dict_data()
            all_data.extend(sample_dict_data)

        self.timepoint_dicts = all_data
        df = pd.DataFrame(all_data)
        if df.empty:
            log_message(
                "No data available to create DataFrame.", self.logs, self.root.name
            )
            return pd.DataFrame(), pd.DataFrame()

        # Mean of all non-string columns grouped by index (time point) to get the average recovery curve for the experiment
        columns_to_average = df.select_dtypes(include=np.number).columns.tolist()
        df_averaged = df.groupby("index")[columns_to_average].mean()

        return df, df_averaged

    def process_sample_data(self):
        df, df_averaged = self.generate_experiment_df()

        if df.empty or df_averaged.empty:
            log_message(
                "No data available to process sample data.", self.logs, self.root.name
            )
            return None

        number_of_pre_bleach_frames = df_averaged[df_averaged.index < 0][
            "normalized_intensity"
        ].count()  # shape[0])

        # Get the count of samples over time and get the index (real index) where there are less than 70% of the maximum samples left
        quantity_threshold = 0.7
        max_samples = df.groupby("index").count()["normalized_intensity"].max()
        count_over_time = df.groupby("index").count()["normalized_intensity"]
        if count_over_time[count_over_time < quantity_threshold * max_samples].empty:
            index_threshold = None
        else:
            index_threshold = count_over_time[
                count_over_time < quantity_threshold * max_samples
            ].index[0]
        self.df_averaged = df_averaged
        self.fit_params, self.fit_values = recovery_fit(
            df_averaged["time"],
            df_averaged["normalized_intensity"],
            low_limit_index=number_of_pre_bleach_frames,
            high_limit_index=index_threshold,
            model="one_phase_fixed_zero",
        )
        if len(df_averaged["time"]) == len(self.fit_values["fitted_intensity"]):
            self.df_averaged["fitted_intensity"] = self.fit_values["fitted_intensity"]
        else:
            print("They dont match")
            # Where time of df_averaged < 0, set fitted intensity to nan then shift the fitted intensity so they start at the first positive time point of the averaged dataframe, then assign the shifted fitted intensity to the averaged dataframe, expect the df_averaged to be larger that the fitted intensity values since the fitted values only start at the first positive time point of the averaged dataframe, so the first part of the fitted intensity values will be nan and then the rest will be the shifted fitted intensity values
            fit_time = self.fit_values["time"]
            fit_intensity = self.fit_values["fitted_intensity"]
            fit_intensity_shifted = np.full_like(df_averaged["time"], np.nan)
            # fit_intensity_shifted[df_averaged["time"] >= fit_time[0]] =

            print("Data frame / series lengths and heads:")
            print(f"len() df_averaged time: {len(df_averaged['time'])}")
            print(f"len() fit_time: {len(fit_time)}")
            print("df_averaged head:")
            print(df_averaged.head(7))
            print(fit_intensity_shifted[:10])

            if len(df_averaged["time"]) == len(fit_intensity_shifted):
                # self.df_averaged["fitted_intensity"] = fit_intensity_shifted
                pass
            else:
                print("They still dont match")
            # TODO: Find a way to shift the fitted intensity values to match the time points of the averaged dataframe without using interpolation, maybe by finding the index of the time point in the averaged dataframe that is closest to each time point in the fitted values and assigning the fitted intensity value to that index in the averaged dataframe

    def generate_report(self, ax=None):

        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(12, 6))
            save_plot = True
        else:
            save_plot = False

        df, df_averaged = self.generate_experiment_df()
        max_time = df_averaged["time"].max()
        next_lower_five = max_time - (max_time % 20)
        time_delta_mean = df_averaged["time_delta"].mean()
        time_points = np.arange(0, next_lower_five, 20)
        time_point_labels = np.round(time_points * time_delta_mean, 2)

        # Make lineplot with index and rename axis with averaged time values
        sns.lineplot(
            df,
            x="index",
            y="normalized_intensity",
            ax=ax,
            label="Individual Samples",
            errorbar="sd",
        )
        sns.lineplot(
            x=self.fit_values["time"].index,
            y=self.fit_values["fitted_intensity"],
            ax=ax,
        )
        ax.axvline(0, color="red", linestyle="--", label="Bleach Time")
        ax.set_xlim(None, max(self.fit_values["time"].index))
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Normalized Intensity")
        ax.set_title("FRAP Recovery Curves")
        ax.set_xticks(time_points)
        ax.set_xticklabels(time_point_labels)
        # ax.set_xticks(np.arange(0, next_lower_five + 5, 5))
        # ax.set_xticklabels(
        #     np.round(np.arange(0, next_lower_five + 5, 5) * time_delta_mean, 2)
        # )
        # # ax.set_xticks(subsample(np.asarray(df_averaged.index), 10))
        # ax.set_xticklabels(subsample(np.asarray(df_averaged["time"].round(2)), 10))

        if save_plot:
            plt.tight_layout()
            fig.savefig(self.root / f"{self.root.name}_experiment_report.pdf")
            plt.show()

    def overview_figure_samples(self):

        number_of_samples = len(self.samples)

        fig = plt.figure(figsize=(21, 3 * number_of_samples))
        gs = gridspec.GridSpec(
            number_of_samples, 7, figure=fig, height_ratios=[1] * number_of_samples
        )
        for i, sample in enumerate(self.samples):
            # try:
            ax0 = fig.add_subplot(gs[i, :2])
            ax1 = fig.add_subplot(gs[i, 2:4])
            ax2 = fig.add_subplot(gs[i, 4])
            ax3 = fig.add_subplot(gs[i, 5])
            ax4 = fig.add_subplot(gs[i, 6])

            sample.generate_report(axes=[ax0, ax1, ax2, ax3, ax4])
            # except Exception as e:
            #     print(f"Error processing sample at {sample.root}: {e}")

        sns.despine()
        plt.tight_layout()
        plt.savefig(self.root / f"{self.name}_overview.pdf")
        plt.close(fig)

    def export_to_hdf(self, hdf_path: Path, experiment_name: Optional[str] = None):

        hdf_path = Path(hdf_path)
        experiment_name = experiment_name or self.metadata["experiment_name"]

        with h5py.File(hdf_path, "a") as hdf_file:
            if experiment_name in hdf_file:
                del hdf_file[experiment_name]

            experiment_group = hdf_file.create_group(experiment_name)

            metadata_group = experiment_group.create_group("meta_data")
            write_dict_to_hdf_group(metadata_group, self.metadata)

            experiment_group.create_group("FrapSamples")

        for sample in self.samples:
            sample.export_to_hdf(hdf_path=hdf_path, experiment_name=experiment_name)

        df, df_averaged = self.generate_experiment_df()

        with h5py.File(hdf_path, "a") as hdf_file:
            experiment_group = hdf_file[experiment_name]
            composed_data_group = experiment_group.create_group("ComposedData")
            data_group = composed_data_group.create_group("Data")

            write_dataframe_to_hdf_group(data_group.create_group("df"), df)
            write_dataframe_to_hdf_group(
                data_group.create_group("df_averaged"), df_averaged
            )
            write_dict_to_hdf_group(
                data_group.create_group("fit_params"), self.fit_params
            )
            write_dict_to_hdf_group(
                data_group.create_group("fit_values"), self.fit_values
            )
