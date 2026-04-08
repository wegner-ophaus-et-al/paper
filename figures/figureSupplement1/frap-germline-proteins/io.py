import numpy as np
from fitting import recovery_fit
import tifffile as tiff
import pandas as pd
from pathlib import Path


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
        lsm_meta = lsm_file.lsm_metadata

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
