import sys

sys.path.append("/Users/icb_remote/Documents/JW/py/packages/")

from gamgee.instance import Marker
import tifffile as tiff
import numpy as np
import re
from datetime import datetime
from pathlib import Path


class TheCell:
    def __init__(
        self, path: str | Path, full_auto=False, conditions=[], model_handler=None
    ):
        self.logs = {}
        path = Path(path) if isinstance(path, str) else path
        self.path = path
        full_name = path.name
        uid, file_name = full_name.split("__", 1)
        file_name_splits = file_name.split("_")
        # Find date structure in the splits (YYYYMMDD or YYYY-MM-DD)

        date_str: str = ""
        cell_condition: str = ""
        for split in file_name_splits:
            if re.match(r"\d{8}", split) or re.match(r"\d{4}-\d{2}-\d{2}", split):
                date = (
                    datetime.strptime(split, "%Y%m%d")
                    if len(split) == 8
                    else datetime.strptime(split, "%Y-%m-%d")
                )
                # Date to iso
                date_str = date.isoformat()
                continue
            for condition in conditions:
                if condition.lower() in split.lower():
                    cell_condition = condition

        self.uid = uid
        self.acquisition_date = date_str
        self.condition = cell_condition
        self.name = file_name

        self.log(
            f"Initialized cell with UID: {self.uid}, acquisition date: {self.acquisition_date}, condition: {self.condition}"
        )

        self.markers = self._populate_markers(model_handler)

    def log(self, message: str):
        self.logs[datetime.now().isoformat()] = f"THECELL: {message}"

    def _populate_markers(self, model_handler):
        images_dir = self.path / "images"
        return {
            "dnd1": Marker(
                image_path=images_dir / "dnd1.tif",
                parent_name=self.name,
                parent_id=self.uid,
                model_handler=model_handler,
                compartment="granules",
            ),
            "gra": Marker(
                image_path=images_dir / "granulito.tif",
                parent_name=self.name,
                parent_id=self.uid,
                model_handler=model_handler,
                compartment="granules",
            ),
            "nucleus": Marker(
                image_path=images_dir / "nls.tif",
                parent_name=self.name,
                parent_id=self.uid,
                model_handler=model_handler,
                compartment="nucleus",
            ),
            "cell": Marker(
                image_path=images_dir / "nls.tif",
                parent_name=self.name,
                parent_id=self.uid,
                model_handler=model_handler,
                compartment="cell",
            ),
        }

    def write_segmentations(self):
        """
        Write segmenatations from TheCell to disk. Funtion mean to be used befor screening and refining segmentations
        """
        mask_dir = self.path / "masks"

        for marker_name, marker in self.markers.items():
            output_path = mask_dir / f"{marker_name}.tif"
            tiff.imwrite(output_path, marker.segmentation.astype(np.uint16))

    def read_segmentations(self):
        """
        Read segmentations from disk to TheCell. Function meant to be used after screening and refining segmentations
        """
        mask_dir = self.path / "masks"

        for marker_name, marker in self.markers.items():
            segmentation_path = mask_dir / f"{marker_name}.tif"
            if segmentation_path.exists():
                marker.segmentation = tiff.imread(segmentation_path)

    def plot_markers_on_axis(self, ax: np.ndarray, blind=True):

        if ax.shape != (3,):
            return

        ax[0].imshow(self.markers["nucleus"].raw_image, cmap="gray")
        ax[1].imshow(self.markers["dnd1"].raw_image, cmap="gray")
        ax[2].imshow(self.markers["gra"].raw_image, cmap="gray")

        ax[0].contour(
            self.markers["nucleus"].segmentation.astype(bool),
            colors="r",
            linewidths=0.5,
        )
        ax[0].contour(
            self.markers["cell"].segmentation.astype(bool), colors="b", linewidths=0.5
        )

        ax[1].imshow(
            self.markers["dnd1"].segmentation,
            cmap="nipy_spectral",
            alpha=self.markers["dnd1"].segmentation > 0,
        )

        ax[2].imshow(
            self.markers["gra"].segmentation,
            cmap="nipy_spectral",
            alpha=self.markers["gra"].segmentation > 0,
        )

        ax[0].text(
            -0.1,
            0.5,
            f"{self.uid}-{self.condition if not blind else ''}",
            ha="center",
            va="center",
            rotation="vertical",
            transform=ax[0].transAxes,
        )
