from pathlib import Path
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.morphology import disk
from skimage.filters import rank, gaussian
from skimage import measure
import numpy as np
import matplotlib.pyplot as plt
from utils import get_pseudo_nucleus_mask, get_measurements


class Cell:
    def __init__(self, path: Path):
        self.path = path
        self.uid, self.name = self.path.name.split("__", 1)
        self.rna_probe, self.condition, self.sample_id = self.name.split("-", 2)
        self.acquisition_date = None
        self.images = {}
        self.segmentations = {}

        self.measurements = []
        self.load_images()

    def load_images(self):
        img = tiff.imread(self.path / "original_raw" / f"{self.name}.lsm")

        self.images = {"granule": img[0], "pla": img[1], self.rna_probe: img[2]}

    def clean_segmentations(self):
        # Delete segmentations outside the cell mask
        self.segmentations["granule"][self.segmentations["cell"] == 0] = 0

    def generate_various_masks(self, granule_periphery_distance=10):
        self.segmentations["granule"] = measure.label(self.segmentations["granule"])
        granule_dilated = (
            ndi.distance_transform_edt(self.segmentations["granule"] == 0)
            < granule_periphery_distance
        )
        self.segmentations["granule_periphery"] = np.bitwise_xor(
            granule_dilated, self.segmentations["granule"].astype(bool)
        )
        self.segmentations["pseudo_nucleus"] = get_pseudo_nucleus_mask(
            self.segmentations["granule"],
            self.images["granule"],
            self.segmentations["cell"],
        )

        self.segmentations["cytoplasm"] = self.segmentations["cell"].astype(bool) & ~(
            granule_dilated | self.segmentations["pseudo_nucleus"].astype(bool)
        )

        self.segmentations["nucleo_cytoplasm"] = (
            self.segmentations["cell"].astype(bool) & ~granule_dilated
        )

    def eval(self):
        measurements_pre_dict = {
            "uid": self.uid,
            "rna_type": self.rna_probe,
            "condition": self.condition,
            "sample_id": self.sample_id,
        }

        intensity_parameters = {
            "sum": np.sum,
            "mean": np.mean,
            "median": np.median,
        }

        collected_measurements = []

        for segmentation_name, segmentation in self.segmentations.items():
            for image_name, image in self.images.items():
                for measurement_name, measurement_func in intensity_parameters.items():
                    measurement_value = get_measurements(
                        measurement_func, image, segmentation
                    )

                    data_point_dict = measurements_pre_dict.copy()
                    data_point_dict.update(
                        {
                            "segmentation_name": segmentation_name,
                            "segmentation_area": np.sum(segmentation > 0),
                            "image_name": image_name
                            if image_name != self.rna_probe
                            else "rna_type",
                            "measurement_name": measurement_name,
                            "value": measurement_value,
                        }
                    )
                    collected_measurements.append(data_point_dict)
        return collected_measurements

    def write_images(self):
        image_path = self.path / "images"
        image_path.mkdir(exist_ok=True)
        for name, image in self.images.items():
            tiff.imwrite(image_path / f"{name}.tif", image)

        segmentation_path = self.path / "segmentations"
        segmentation_path.mkdir(exist_ok=True)
        for name, segmentation in self.segmentations.items():
            tiff.imwrite(segmentation_path / f"{name}.tif", segmentation)

    def read_images(self):
        for name in ["granule", "pla", self.rna_probe]:
            self.images[name] = tiff.imread(self.path / "images" / f"{name}.tif")

        for name in ["granule", "cell"]:
            self.segmentations[name] = measure.label(
                tiff.imread(self.path / "segmentations" / f"{name}.tif")
            )

    def segment(self, segmentation_model):
        if self.images.get("granule", None) is not None:
            self.segmentations["granule"] = segmentation_model.segment(
                self.images["granule"]
            )
        else:
            self.load_images()
            self.segmentations["granule"] = segmentation_model.segment(
                self.images["granule"]
            )

        cell_img = self.images.get("granule", None)

        # Subtract the background using a rank mean filter
        # footprint = disk(20)
        # normal_result = rank.mean(cell_img, footprint=footprint)
        # cell_img = cell_img - normal_result
        cell_img = gaussian(cell_img, sigma=3)
        cell_mask = cell_img > np.percentile(cell_img, 45)
        cell_mask = ndi.binary_opening(cell_mask)
        cell_mask = ndi.binary_fill_holes(cell_mask)
        self.segmentations["cell"] = cell_mask.astype(np.uint8)

    def plot_images(self, ax):

        if ax is None:
            fig, ax = plt.subplots(1, 3, figsize=(15, 5))
            ax_provided = False
        else:
            ax_provided = True

        ax[0].imshow(
            self.images["granule"],
            cmap="gray",
            vmin=0,
            vmax=np.percentile(self.images["granule"], 95),
        )
        ax[0].set_title("Granule")
        ax[0].axis("off")
        ax[0].contour(
            self.segmentations["granule"].astype(bool), colors="r", linewidths=0.5
        )
        ax[0].contour(
            self.segmentations["cell"].astype(bool), colors="b", linewidths=0.5
        )
        ax[0].text(
            -0.05,
            0.5,  # x slightly outside the left edge, y centered
            f"{self.uid}-{self.condition}",
            transform=ax[0].transAxes,
            rotation=90,
            va="center",
            ha="right",
        )
        ax[1].imshow(
            self.images["pla"],
            cmap="inferno",
            vmin=0,
            vmax=np.percentile(self.images["pla"], 85),
        )
        ax[1].set_title("PLA")
        ax[1].axis("off")
        ax[1].contour(
            self.segmentations["granule"].astype(bool), colors="r", linewidths=0.5
        )
        ax[2].imshow(
            self.images[self.rna_probe],
            cmap="viridis",
            vmin=0,
            vmax=np.percentile(self.images[self.rna_probe], 95),
        )
        ax[2].set_title(self.rna_probe, fontstyle="italic")
        ax[2].axis("off")

        if not ax_provided:
            plt.show()

    def plot_segmentations(self, ax):
        if ax is None:
            fig, ax = plt.subplots(1, len(self.segmentations), figsize=(15, 5))
            ax_provided = False
        else:
            ax_provided = True

        for i, (name, segmentation) in enumerate(self.segmentations.items()):
            if i == 0:
                ax[i].set_ylabel(self.uid)
            ax[i].imshow(self.images["granule"], cmap="gray")
            ax[i].contour(segmentation.astype(bool), colors="r", linewidths=0.5)
            ax[i].imshow(segmentation, cmap="tab20", alpha=0.5)
            ax[i].set_title(name)
            ax[i].axis("off")

        if not ax_provided:
            plt.show()
