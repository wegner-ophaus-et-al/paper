import numpy as np
import tifffile as tiff
import pandas as pd
from pathlib import Path
from typing import Optional
import h5py


def lsm_metadata(lsm_file_path: Path, raw_time=False):
    """
    Read meta_data from a Zeiss LSM file, extracting bleach event timestamps and calculating relative time, bleach duration, mean timestamp difference, and pre-bleach frame count.
    Parameters:
    lsm_file_path (Path): Path to the LSM file.
    raw_time (bool): If True, include raw timestamps and bleach event times in the output

    Returns:
    meta_data (dict): A dictionary containing relative time, bleach duration, mean timestamp difference, pre-bleach frame count, and optionally raw timestamps and bleach event times.
    """
    meta_data = {}

    with tiff.TiffFile(lsm_file_path) as lsm_file:
        lsm_meta = lsm_file.lsm_metadata if lsm_file.lsm_metadata else {}

        # Take base measurements
        meta_data.update(
            {
                "x_size": lsm_meta.get("DimensionX", None),
                "y_size": lsm_meta.get("DimensionY", None),
                "time_size": lsm_meta.get("DimensionTime", None),
                "dtype_raw": lsm_meta.get("DataType", None),
                "voxel_size_x": lsm_meta.get("VoxelSizeX", None),
                "voxel_size_y": lsm_meta.get("VoxelSizeY", None),
            }
        )

        bleach_start_time = 0
        bleach_end_time = 0
        for stamp, onoff, __ in lsm_meta["EventList"]:
            if onoff == 2:
                bleach_start_time = stamp
            elif onoff == 3:
                bleach_end_time = stamp

        relative_timestamps = lsm_meta["TimeStamps"] - bleach_end_time
        bleach_duration = bleach_end_time - bleach_start_time

        positive_timestamps = relative_timestamps[relative_timestamps > 0]
        mean_timestamp_difference = np.diff(positive_timestamps).mean()

        count_pre_bleach_frames = np.sum(relative_timestamps < 0)

    if raw_time:
        meta_data["raw_time"] = lsm_meta["TimeStamps"]
        meta_data["bleach_start_time"] = bleach_start_time
        meta_data["bleach_end_time"] = bleach_end_time

    meta_data.update(
        {
            "time": relative_timestamps,
            "bleach_duration": bleach_duration,
            "mean_timestamp_difference": mean_timestamp_difference,
            "pre_bleach_frame_count": count_pre_bleach_frames,
        }
    )

    return meta_data


def normalize_image(image, bit_depth=16, saturation_percentile=100):
    """
    Normalize image data to a specified bit depth, with optional saturation based on a percentile.
    Parameters:
    image (ndarray): Input image data as a NumPy array.
    bit_depth (int): Desired bit depth for normalization (default is 16).
    saturation_percentile (float): Percentile for saturation (default is 100 which is the max value).

    Returns:
    normalized_image (ndarray): Normalized image data as a NumPy array.
    """
    max_value = (
        np.percentile(image, saturation_percentile)
        if saturation_percentile < 100
        else image.max()
    )
    normalized_image = np.clip(image, 0, max_value) / max_value * (2**bit_depth - 1)
    if bit_depth == 16:
        return normalized_image.astype(np.uint16)
    elif bit_depth == 8:
        return normalized_image.astype(np.uint8)
    elif bit_depth == 12:
        print(
            "Warning: 12-bit images are often stored in 16-bit containers. Normalizing to 12 bits will still return a 16-bit image with values scaled to the 12-bit range."
        )
        return normalized_image.astype(
            np.uint16
        )  # 12-bit images are often stored in 16-bit containers
    elif bit_depth == 32:
        return normalized_image.astype(np.float32)
    elif bit_depth == 64:
        return normalized_image.astype(np.float64)
    else:
        raise ValueError(
            "Unsupported bit depth. Supported values are 8, 16, 32, and 64."
        )


def read_image(image_file_path: Path, saturation_percentile=100):
    """
    Read image data from drift corrected image and normalize in 16 bit
    """
    return normalize_image(
        tiff.imread(image_file_path),
        bit_depth=16,
        saturation_percentile=saturation_percentile,
    )


def read_mask(mask_file_path: Path):
    """
    Read mask data from a TIFF file, returning it as a binary NumPy array.
    Parameters:
    mask_file_path (Path): Path to the mask TIFF file.

    Returns:
    mask (ndarray): Binary mask as a NumPy array where non-zero values are set to 1.
    """
    return (tiff.imread(mask_file_path) > 0).astype(np.uint8)


def find_file_paths(dir: Path):
    """
    Search for all files needed for a FRAP analysis with Phair double normalization
    """

    file_paths: dict[str, Optional[Path]] = {
        "lsm_path": None,
        "image_path": None,
        "irradiated_mask_path": None,
        "correction_mask_path": None,
        "background_mask_path": None,
    }

    file_paths["lsm_path"] = next(dir.glob("*.lsm"))
    tif_path = next(dir.glob("*.tif"))
    if tif_path.stem != file_paths["lsm_path"].stem:
        print(
            f"WARNING: LSM and TIFF file names do not match: {file_paths['lsm_path'].name} vs {tif_path.name}"
        )

    file_paths["image_path"] = tif_path

    masks_subfolder = dir / "masks"

    if not masks_subfolder.exists() or not masks_subfolder.is_dir():
        # Search for subdolder that contains "mask" in its name
        masks_subfolder = next(
            (
                sub
                for sub in dir.iterdir()
                if sub.is_dir() and "mask" in sub.name.lower()
            ),
            None,
        )

    if masks_subfolder is None:
        raise FileNotFoundError(f"No 'masks' subfolder found in {dir}")

    irradiated_mask_path = next(masks_subfolder.glob("*irrad*.tif"), None)
    correction_mask_path = next(masks_subfolder.glob("*correction*.tif"), None)
    background_mask_path = next(masks_subfolder.glob("*background*.tif"), None)

    for mask_type, mask_path in zip(
        ["irradiated_mask_path", "correction_mask_path", "background_mask_path"],
        [irradiated_mask_path, correction_mask_path, background_mask_path],
    ):
        if mask_path is None or not mask_path.exists():
            print(
                f"WARNING: No {mask_type} found in {masks_subfolder} or does not exist."
            )
        else:
            file_paths[mask_type] = mask_path

    for key, path in file_paths.items():
        if path is None:
            raise FileNotFoundError(f"{key} is required but was not found in {dir}")

    return file_paths
