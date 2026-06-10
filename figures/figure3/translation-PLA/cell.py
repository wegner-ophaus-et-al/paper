from pathlib import Path
import tifffile as tiff
from scipy import ndimage as ndi
from skimage.morphology import disk
from skimage.filters import rank, gaussian
import numpy as np
import matplotlib.pyplot as plt


class Cell:
    def __init__(self, path: Path):
        self.path = path
        self.uid, self.name = self.path.name.split("__", 1)
        self.rna_probe, self.condition, self.sample_id = self.name.split("-", 2)
        self.acquisition_date = None
        self.images = {}
        self.segmentations = {}

        self.load_images()

    def load_images(self):
        img = tiff.imread(self.path / "original_raw" / f"{self.name}.lsm")

        self.images = {"granule": img[0], "pla": img[1], self.rna_probe: img[2]}

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
            self.segmentations[name] = tiff.imread(
                self.path / "segmentations" / f"{name}.tif"
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
        cell_img = gaussian(cell_img, sigma=1.5)
        cell_mask = cell_img > np.percentile(cell_img, 20)
        cell_mask = ndi.binary_opening(cell_mask)
        self.segmentations["cell"] = cell_mask.astype(np.uint8)

    def plot_images(self, ax):

        if ax is None:
            fig, ax = plt.subplots(1, 3, figsize=(15, 5))
            ax_provided = False
        else:
            ax_provided = True

        ax[0].imshow(self.images["granule"], cmap="gray")
        ax[0].set_title("Granule")
        ax[0].axis("off")
        ax[0].contour(
            self.segmentations["granule"].astype(bool), colors="r", linewidths=0.5
        )
        ax[0].contour(
            self.segmentations["cell"].astype(bool), colors="b", linewidths=0.5
        )
        ax[1].imshow(self.images["pla"], cmap="magma")
        ax[1].set_title("PLA")
        ax[1].axis("off")
        ax[1].contour(
            self.segmentations["granule"].astype(bool), colors="r", linewidths=0.5
        )
        ax[2].imshow(self.images[self.rna_probe], cmap="viridis")
        ax[2].set_title(self.rna_probe, fontstyle="italic")
        ax[2].axis("off")

        if not ax_provided:
            plt.show()
