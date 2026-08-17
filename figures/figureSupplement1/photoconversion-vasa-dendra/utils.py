import numpy as np
from matplotlib.patches import Rectangle
from tifffile import TiffFile


def scale_bar(ax, length, thickness, pixel_size, pad=0.03, color="white"):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    sx, sy = np.sign(x1 - x0), np.sign(y1 - y0)
    px = length / pixel_size
    x = x1 - sx * (pad * abs(x1 - x0) + px)
    y = y0 + sy * pad * abs(y1 - y0)
    ax.add_patch(Rectangle((x, y), sx * px, sy * thickness, color=color, zorder=10))


def lsm_pixel_size(path):
    """Return pixel size in microns, asserting square pixels."""
    with TiffFile(path) as tif:
        m = tif.lsm_metadata
        x, y = m["VoxelSizeX"] * 1e6, m["VoxelSizeY"] * 1e6
    assert abs(x - y) < 1e-9, f"non-square pixels: {x} != {y}"
    return x


def bleach_times(lsm_meta):
    start_time = end_time = None
    for stamp, event_type, _ in lsm_meta["EventList"]:
        if event_type == 2:
            start_time = stamp
        elif event_type == 3:
            end_time = stamp
    return start_time, end_time


def get_lsm_meta(path):
    with TiffFile(path) as tif:
        return tif.lsm_metadata
