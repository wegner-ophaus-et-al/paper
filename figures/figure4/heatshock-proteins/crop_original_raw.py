import numpy as np
from pathlib import Path
import tifffile as tiff
import hashlib
import shutil


def get_uid_from_image(img):
    hash_sha = hashlib.sha256()
    hash_sha.update(img.tobytes())
    return hash_sha.hexdigest()[:8]


def normalize_image(img, bit_depth=16):
    img_norm = img / img.max()
    return (img_norm * (2**bit_depth - 1)).astype(f"uint{bit_depth}")


def crop_image_with_mask(image, mask):
    mask = mask > 0

    # Find bouding box of non zero pixels in the mask
    y_indices, x_indices = np.where(mask)
    y_min, y_max = y_indices.min(), y_indices.max()
    x_min, x_max = x_indices.min(), x_indices.max()

    return image[:, :, y_min : y_max + 1, x_min : x_max + 1]


projection = "max"  # or "avg", "min", "sum"

root = Path("/Users/julian/local_files/20260505_heat-shock/")

constucts = ["F617", "F618"]

output_folder_name = "cropped_images"

channel_resolve_map = {0: "nls", 1: "stress", 2: "gra"}


for construct_name in constucts:
    construct_folder = root / construct_name

    sample_folders = [
        p
        for p in construct_folder.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    ]
    for sample_folder in sample_folders:
        sample_id, sampel_name = sample_folder.name.split("__", 1)
        masks_dir = sample_folder / "masks"
        masks_path_list = [
            mf for mf in masks_dir.glob("*.tif") if not mf.name.startswith(".")
        ]

        original_raw_path = (sample_folder / "original_raw").glob("*.tif").__next__()
        img = tiff.imread(original_raw_path)

        for rolling_idx, mask_path in enumerate(masks_path_list):
            mask_img = tiff.imread(mask_path)
            crop_img = crop_image_with_mask(img, mask_img)

            uid = get_uid_from_image(crop_img)
            out_folder = (
                root
                / output_folder_name
                / construct_name
                / f"{uid}__{sampel_name}-{rolling_idx}"
            )
            img_out_dir = out_folder / "images"
            masks_out_dir = out_folder / "masks"

            # Make a few dirs
            out_folder.mkdir(parents=True, exist_ok=True)
            (out_folder / "original_raw").mkdir(exist_ok=True)
            img_out_dir.mkdir(exist_ok=True)
            masks_out_dir.mkdir(exist_ok=True)

            # copy original raw image
            shutil.copy2(
                original_raw_path, out_folder / "original_raw" / original_raw_path.name
            )
            shutil.copy2(mask_path, masks_out_dir / mask_path.name)

            tiff.imwrite(img_out_dir / "raw_crop.tif", crop_img)

            if projection == "max":
                img_proj = np.max(crop_img, axis=0)
            elif projection == "avg":
                img_proj = np.mean(crop_img, axis=0)
            elif projection == "min":
                img_proj = np.min(crop_img, axis=0)
            elif projection == "sum":
                img_proj = np.sum(crop_img, axis=0)
            else:
                raise ValueError(f"Unknown projection type: {projection}")

            for ch_idx, name in channel_resolve_map.items():
                tiff.imwrite(
                    img_out_dir / f"{name}.tif", normalize_image(img_proj[ch_idx])
                )
