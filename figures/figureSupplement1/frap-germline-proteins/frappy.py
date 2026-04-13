import numpy as np
from numpy.typing import NDArray
from pathlib import Path
from ios import lsm_metadata, read_image, read_mask, find_file_paths
from fitting import recovery_fit
from calculator import phair_double_normalization
import seaborn as sns
import matplotlib.pyplot as plt


class FrapSample:
    def __init__(self, root_dir: Path):
        self.root = root_dir

        self.image_paths: dict = find_file_paths(self.root)

        self.metadata = lsm_metadata(self.image_paths["lsm_path"])

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
        self.masks["reference"] = read_mask(self.image_paths["correction_mask_path"])
        self.masks["background"] = read_mask(self.image_paths["background_mask_path"])

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

        if len(self.metadata["time"]) != len(
            normalized_intensity_dict["intensity_normalized"]
        ):
            raise ValueError(
                f"Time points and normalized intensity data do not have the same length: {len(self.metadata['time'])} vs {len(self.data['normalized_intensity'])}"
            )

        self.data.update(
            {
                "time": self.metadata["time"],
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
