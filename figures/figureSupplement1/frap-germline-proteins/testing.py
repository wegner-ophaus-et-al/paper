import tifffile as tiff
import numpy as np
from pathlib import Path

root = Path(
    "/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_full_vasa/10hpf/20220725_p3"
)

lsm_file_path = Path()
raw_time = False


def lsm_metadata(root: Path, raw_time=False):
    for lsm_path in root.glob("*.lsm"):
        if not lsm_path.name.startswith("."):
            lsm_file_path = lsm_path

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
