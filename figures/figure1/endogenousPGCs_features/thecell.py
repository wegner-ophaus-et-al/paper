import sys
import copy

sys.path.append("/Users/icb_remote/Documents/JW/py/packages/")

from gamgee.instance import Marker
import gamgee.features as features
import tifffile as tiff
import numpy as np
import re
from datetime import datetime
from pathlib import Path


class TheCell:
    def __init__(
        self,
        path: str | Path,
        full_auto=False,
        conditions=[],
        model_handler=None,
        model_handler_id: str | None = None,
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
        dev_stage: str = ""
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

            if "hpf" in split.lower():
                dev_stage = split

        self.uid = uid
        self.acquisition_date = date_str
        self.condition = cell_condition
        self.stage = dev_stage
        self.name = file_name

        self.log(
            f"Initialized cell with UID: {self.uid}, acquisition date: {self.acquisition_date}, condition: {self.condition}"
        )

        self.model_handler_id = model_handler_id
        self.markers = self._populate_markers(model_handler)

    def attach_model_handler(self, model_handler, model_handler_id: str | None = None):
        if model_handler_id is not None:
            self.model_handler_id = model_handler_id
        for marker in self.markers.values():
            if hasattr(marker, "sam_model"):
                marker.set_sam_model(model_handler)

    def __getstate__(self):
        state = self.__dict__.copy()
        markers = state.get("markers")
        if markers:
            stripped_markers = {}
            for name, marker in markers.items():
                marker_copy = copy.copy(marker)
                # for k, v in marker_copy.__dict__.items():
                #     print(f"{k}\t {v}")
                #     print("\n")
                if hasattr(marker_copy, "sam_model"):
                    marker_copy.sam_model = None
                stripped_markers[name] = marker_copy
            state["markers"] = stripped_markers
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if not hasattr(self, "model_handler_id"):
            self.model_handler_id = None

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

    def plot_markers_on_axis(
        self,
        ax: np.ndarray,
        blind=True,
        granule_percentile=99,
        segmentation_cmap="nipy_spectral",
        granule_alpha=0.8,
    ):

        if ax.shape != (3,):
            return

        ax[0].imshow(self.markers["nucleus"].raw_image, cmap="gray")
        ax[1].imshow(
            self.markers["dnd1"].raw_image,
            cmap="gray",
            vmax=np.percentile(self.markers["dnd1"].raw_image, granule_percentile),
        )
        ax[2].imshow(
            self.markers["gra"].raw_image,
            cmap="gray",
            vmax=np.percentile(self.markers["gra"].raw_image, granule_percentile),
        )

        ax[0].contour(
            self.markers["nucleus"].segmentation.astype(bool),
            colors="#f0be4a",
            linewidths=0.5,
        )
        ax[0].contour(
            self.markers["cell"].segmentation.astype(bool),
            colors="#3f83bf",
            linewidths=0.5,
        )

        ax[1].imshow(
            self.markers["dnd1"].segmentation,
            cmap=segmentation_cmap,
            alpha=(self.markers["dnd1"].segmentation > 0).astype(float) * granule_alpha,
        )

        ax[2].imshow(
            self.markers["gra"].segmentation,
            cmap=segmentation_cmap,
            alpha=(self.markers["gra"].segmentation > 0).astype(float) * granule_alpha,
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

    def clean_segmentations(self):

        # Delete segmentations outside of the cell segmentation
        for marker_name in ["dnd1", "gra", "nucleus"]:
            self.markers[marker_name].segmentation[
                self.markers["cell"].segmentation == 0
            ] = 0

    def compute_features(self):
        data_collector = {
            "uid": self.uid,
            "condition": self.condition,
            "stage": self.stage,
        }

        # Pass marker features to data collector
        data_collector.update(
            {
                "granule_features": {
                    "dnd1": self.markers["dnd1"].get_features(),
                    "gra": self.markers["gra"].get_features(),
                }
            }
        )

        # Get basic cell and nucleus features
        data_collector.update(
            {
                "cell": {
                    "Area": np.sum(self.markers["cell"].segmentation > 0),
                    "SphericalVolume": features.spherical_volume(
                        self.markers["cell"].segmentation > 0
                    ),
                },
                "nucleus": {
                    "Area": np.sum(self.markers["nucleus"].segmentation > 0),
                    "SphericalVolume": features.spherical_volume(
                        self.markers["nucleus"].segmentation > 0
                    ),
                },
            }
        )

        # Get nuclear distance features
        for marker_name in ["dnd1", "gra"]:
            data_collector["granule_features"][marker_name][
                "NuclearDistanceFeatures"
            ] = features.nuclear_distance_features(
                self.markers[marker_name].segmentation,
                self.markers["nucleus"].segmentation,
            )

        # Get dnd1 and gra co-localization features
        coloc_features = {
            "MandersPercentile": features.manders_across_percentiles(
                self.markers["dnd1"].raw_image,
                self.markers["gra"].raw_image,
                mask=self.markers["cell"].segmentation > 0,
            ),
            "MandersSegmentationsM1": features.manders(
                self.markers["dnd1"].raw_image,
                self.markers["gra"].segmentation > 0,
                mask=self.markers["cell"].segmentation > 0,
            ),
            "MandersSegmentationsM2": features.manders(
                self.markers["gra"].raw_image,
                self.markers["dnd1"].segmentation > 0,
                mask=self.markers["cell"].segmentation > 0,
            ),
            "Pearson": features.pearson(
                self.markers["dnd1"].raw_image,
                self.markers["gra"].raw_image,
                mask=self.markers["cell"].segmentation > 0,
            ),
            "IoU": features.iou(
                self.markers["dnd1"].segmentation > 0,
                self.markers["gra"].segmentation > 0,
            ),
        }
        data_collector["granule_features"]["colocalization"] = coloc_features

        self.features = data_collector
