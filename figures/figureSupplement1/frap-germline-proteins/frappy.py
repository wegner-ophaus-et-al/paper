import numpy as np
from numpy.typing import NDArray
from pathlib import Path
from io import lsm_metadata, read_image, read_mask, find_file_paths


class FrapSample:
    def __init__(self, root_dir: Path):
        self.root = root_dir

        self.image_paths = find_file_paths(self.root)

        self.metadata = lsm_metadata(self.image_paths["lsm_path"])

        # Prep vars for images and masks for later population
        self.image_sequence = None
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
