import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pathlib import Path
from ios import lsm_metadata, read_image, read_mask, find_file_paths, find_sample_dirs
from fitting import recovery_fit
from calculator import phair_double_normalization, subsample
import seaborn as sns
import matplotlib.pyplot as plt
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
            palette="gray",
        )
        sns.lineplot(
            x=self.data["time"],
            y=self.data["fit_values"]["fitted_intensity"],
            ax=axes[0],
            label=f"Fit ({self.fit_type})",
            palette="red",
        )
        axes[0].set_title(
            f"FRAP Recovery Curve (tau = {self.data['fit_params']['tau']:.2f} s)"
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
        axes[2].contour(self.masks["irradiated"], colors="r", linewidths=0.5)
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
        pass


class FrapExperiment:
    def __init__(self, root_dir: Path):
        self.root = root_dir
        self.logs = {}
        self.sample_paths = find_sample_dirs(root_dir)
        self.samples = []
        self.metadata = {
            "experiment_name": self.root.name,
            "experiment_path": self.root,
        }
        self.fit_params = {}
        self.fit_values = {}

        for sample_path in self.sample_paths:
            log_message("Loading sample...", self.logs, sample_path.name)
            try:
                self.samples.append(FrapSample(sample_path))
                self.logs.update({k: v for k, v in self.samples[-1].logs.items()})
                log_message("Sample loaded successfully.", self.logs, sample_path.name)
            except Exception as error:
                log_message(
                    f"Failed to load sample: {error}", self.logs, sample_path.name
                )

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

        df = pd.DataFrame(all_data)

        df_averaged = df.groupby("index").mean()

        return df, df_averaged

    def process_sample_data(self):
        df, df_averaged = self.generate_experiment_df()

        number_of_pre_bleach_frames = df_averaged[df_averaged.index < 0][
            "normalized_intensity"
        ].count()  # shape[0])

        # Get the count of samples over time and get the index (real index) where there are less than 70% of the maximum samples left
        quantity_threshold = 0.7
        max_samples = df.groupby("index").count()["normalized_intensity"].max()
        count_over_time = df.groupby("index").count()["normalized_intensity"]
        index_threshold = count_over_time[
            count_over_time < quantity_threshold * max_samples
        ].index[0]

        print(
            df_averaged.loc[index_threshold],
            " is the time point where there are less than 70% of the maximum samples left",
        )
        self.fit_params, self.fit_values = recovery_fit(
            df_averaged["time"],
            df_averaged["normalized_intensity"],
            low_limit_index=number_of_pre_bleach_frames,
            high_limit_index=index_threshold,
        )
        print(self.fit_params)

    def generate_report(self):
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        df, df_averaged = self.generate_experiment_df()

        # Make lineplot with index and rename axis with averaged time values
        sns.lineplot(
            df,
            x="index",
            y="normalized_intensity",
            ax=axes[0],
            label="Individual Samples",
        )
        sns.lineplot(
            x=self.fit_values["time"].index,
            y=self.fit_values["fitted_intensity"],
            ax=axes[0],
        )
        axes[0].axvline(0, color="red", linestyle="--", label="Bleach Time")
        axes[0].set_xlim(None, max(self.fit_values["time"].index))
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Normalized Intensity")
        axes[0].set_title("FRAP Recovery Curves")
        axes[0].set_xticks(subsample(np.asarray(df_averaged.index), 10))
        axes[0].set_xticklabels(subsample(np.asarray(df_averaged["time"].round(2)), 10))

        # On the second subplot, plot the number of samoles contributing to each time point as a bar plot
        sns.barplot(
            data=df.groupby("index").count().reset_index(),
            x="index",
            y="normalized_intensity",
            ax=axes[1],
        )
        axes[1].set_xlabel("Time (s)")
        axes[1].set_ylabel("Number of Samples")
        axes[1].set_title("Sample Contribution Over Time")

        plt.show()
