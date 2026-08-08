import struct
from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import tifffile as tiff
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import median_abs_deviation
from skimage.measure import label
from skimage.registration import phase_cross_correlation
from skimage.feature import blob_dog
from skimage.draw import disk, rectangle, polygon

LSM_PATH = (
    Path(__file__).parent / "data" / "p11_directHit_bigGranule_same_cell_as_p10.lsm"
)
RESULTS_PATH = Path(__file__).parent / "results.csv"

MAX_DISTANCE_PX = 15  # max px a granule may move between adjacent frames to keep its ID

RECTANGLE, ELLIPSE, CLOSED_POLYLINE, OPEN_POLYLINE, CIRCLE = 18, 19, 20, 21, 24


def read_lsm(lsm_path):
    with tiff.TiffFile(lsm_path) as tif:
        stack = tif.asarray()  # (T, C, Y, X)
        lsm_meta = tif.lsm_metadata
    return stack, lsm_meta


def drift_correct(stack):
    """Register every frame's channel 1 to frame 0's channel 1; apply the same shift to both channels."""
    reference = stack[0, 0]
    corrected = np.empty_like(stack)
    corrected[0] = stack[0]
    shifts = [(0.0, 0.0)]
    for t in range(1, stack.shape[0]):
        shift, _, _ = phase_cross_correlation(reference, stack[t, 0])
        shifts.append(tuple(shift))
        for c in range(stack.shape[1]):
            corrected[t, c] = ndimage.shift(stack[t, c], shift)
    return corrected, shifts


def bleach_times(lsm_meta):
    start_time = end_time = None
    for stamp, event_type, _ in lsm_meta["EventList"]:
        if event_type == 2:
            start_time = stamp
        elif event_type == 3:
            end_time = stamp
    return start_time, end_time


def _to_px(value, voxel_size):
    """User-drawn overlay shapes (VectorOverlay/BleachRoi) are stored directly in pixel
    coordinates; the scan-restriction ROI (Roi) is stored in physical units (meters) and
    needs dividing by the pixel size instead.
    """
    return value if voxel_size is None else value / voxel_size


def _read_shape_header(fh):
    """Read the common per-shape header fields, up to (not including) the shape-specific geometry.

    Field layout per Bio-Formats' ZeissLSMReader.parseOverlays.
    """
    shape_type = struct.unpack("<i", fh.read(4))[0]
    fh.read(4)  # block length
    fh.read(4)  # line width (stored as int32)
    fh.read(4)  # measurements
    fh.read(16)  # text offset x, y (2 x float64)
    fh.read(4)  # color
    fh.read(4)  # valid flag
    fh.read(4)  # knot width
    fh.read(4)  # catch area
    fh.read(4 * 13)  # font height/width/escapement/orientation/weight/italic/
    #                   underlined/strikeout/charset/output precision/clip precision/
    #                   quality/pitch-and-family
    fh.read(64)  # font name
    fh.read(2)  # enabled (int16)
    fh.read(4)  # moveable
    fh.read(34)  # reserved
    return shape_type


def read_bleach_roi_mask(
    lsm_path, offset, width, height, voxel_size_x=None, voxel_size_y=None
):
    """Parse the LSM bleach ROI overlay at `offset` into a boolean (height, width) mask.

    Format reverse-engineered from Bio-Formats' ZeissLSMReader.parseOverlays. Only
    rectangle/ellipse/circle/polyline shapes are supported; anything else raises.
    """
    if not offset:
        return None

    mask = np.zeros((height, width), dtype=bool)
    with open(lsm_path, "rb") as fh:
        fh.seek(offset)
        number_of_shapes = struct.unpack("<i", fh.read(4))[0]
        size = struct.unpack("<i", fh.read(4))[0]
        if size <= 194:
            return mask
        fh.read(20)
        fh.read(4)  # valid flag
        fh.read(164)

        for _ in range(number_of_shapes):
            shape_type = _read_shape_header(fh)

            if shape_type == RECTANGLE:
                fh.read(4)
                x0, y0, x1, y1 = struct.unpack("<4d", fh.read(32))
                r0, r1 = sorted((_to_px(y0, voxel_size_y), _to_px(y1, voxel_size_y)))
                c0, c1 = sorted((_to_px(x0, voxel_size_x), _to_px(x1, voxel_size_x)))
                rr, cc = rectangle((r0, c0), end=(r1, c1), shape=mask.shape)
                mask[rr, cc] = True

            elif shape_type == CIRCLE:
                fh.read(4)
                cx, cy, ex, ey = struct.unpack("<4d", fh.read(32))
                cx_px, cy_px = _to_px(cx, voxel_size_x), _to_px(cy, voxel_size_y)
                ex_px, ey_px = _to_px(ex, voxel_size_x), _to_px(ey, voxel_size_y)
                radius = np.hypot(ex_px - cx_px, ey_px - cy_px)
                rr, cc = disk((cy_px, cx_px), radius, shape=mask.shape)
                mask[rr, cc] = True

            elif shape_type in (ELLIPSE, CLOSED_POLYLINE, OPEN_POLYLINE):
                n_knots = struct.unpack("<i", fh.read(4))[0]
                coords = struct.unpack(f"<{2 * n_knots}d", fh.read(16 * n_knots))
                xs = [_to_px(coords[2 * i], voxel_size_x) for i in range(n_knots)]
                ys = [_to_px(coords[2 * i + 1], voxel_size_y) for i in range(n_knots)]
                rr, cc = polygon(ys, xs, shape=mask.shape)
                mask[rr, cc] = True

            else:
                raise NotImplementedError(
                    f"Unsupported bleach ROI shape type: {shape_type}"
                )

    return mask


def detect_spots(frame, k=5, min_sigma=1, max_sigma=5, threshold=0.1):
    """DoG candidate spots (skimage.feature.blob_dog), kept if intensity > median + k*MAD."""
    frame = frame.astype(float)
    normalized = frame / frame.max()
    blobs = blob_dog(
        normalized, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold
    )
    if len(blobs) == 0:
        return blobs

    median = np.median(frame)
    mad = median_abs_deviation(frame, axis=None)
    rows = blobs[:, 0].astype(int).clip(0, frame.shape[0] - 1)
    cols = blobs[:, 1].astype(int).clip(0, frame.shape[1] - 1)
    intensities = frame[rows, cols]
    return blobs[intensities > median + k * mad]


def measure_spot_intensities(frame, blobs):
    means = []
    for y, x, sigma in blobs:
        radius = sigma * np.sqrt(2)
        rr, cc = disk((y, x), radius, shape=frame.shape)
        means.append(frame[rr, cc].mean())
    return means


def link_granules(
    prev_blobs, prev_ids, curr_blobs, next_id, max_distance_px=MAX_DISTANCE_PX
):
    """Assign granule IDs to curr_blobs via optimal one-to-one matching (Hungarian
    algorithm) against prev_blobs, gated by max_distance_px. Returns (curr_ids, next_id).
    """
    n_curr = len(curr_blobs)
    curr_ids = np.full(n_curr, -1, dtype=int)

    if len(prev_blobs) > 0 and n_curr > 0:
        cost = cdist(prev_blobs[:, :2], curr_blobs[:, :2])
        prev_idx, curr_idx = linear_sum_assignment(cost)
        for p, c in zip(prev_idx, curr_idx):
            if cost[p, c] <= max_distance_px:
                curr_ids[c] = prev_ids[p]

    unmatched = curr_ids == -1
    n_new = unmatched.sum()
    curr_ids[unmatched] = np.arange(next_id, next_id + n_new)
    next_id += n_new

    return curr_ids, next_id


def draw_blob_contours(ax, blobs, color="lime"):
    """Overlay each blob's detection circle (y, x, sigma) on ax."""
    for y, x, sigma in blobs:
        ax.add_patch(
            Circle((x, y), sigma * np.sqrt(2), color=color, fill=False, linewidth=0.5)
        )


def main():
    stack, lsm_meta = read_lsm(LSM_PATH)
    corrected, shifts = drift_correct(stack)
    tiff.imwrite(
        LSM_PATH.with_suffix(".corrected.tif"),
        corrected.astype(np.float32),
        imagej=True,
    )
    bleach_start, bleach_end = bleach_times(lsm_meta)
    frame_before_bleach = np.searchsorted(lsm_meta["TimeStamps"], bleach_start) - 1
    frame_after_bleach = np.searchsorted(lsm_meta["TimeStamps"], bleach_end) + 1
    height, width = stack.shape[-2:]
    voxel_size_x, voxel_size_y = lsm_meta["VoxelSizeX"], lsm_meta["VoxelSizeY"]

    # try, in order: the dedicated bleach ROI, the user-drawn overlay (both in pixel
    # units), and finally the hardware scan-restriction ROI (physical units, meters) --
    # some acquisitions only populate one of these with the actual bleached region.
    bleach_mask = None
    for offset_key, unit in (
        ("OffsetBleachRoi", "px"),
        ("OffsetVectorOverlay", "px"),
        ("OffsetRoi", "physical"),
    ):
        vx, vy = (None, None) if unit == "px" else (voxel_size_x, voxel_size_y)
        bleach_mask = read_bleach_roi_mask(
            LSM_PATH, lsm_meta.get(offset_key, 0), width, height, vx, vy
        )
        if bleach_mask is not None and bleach_mask.any():
            break

    timestamps = np.asarray(lsm_meta["TimeStamps"])

    rows = []
    prev_blobs = np.empty((0, 3))
    prev_ids = np.empty((0,), dtype=int)
    next_id = 0

    blob_segmentations = {}
    blobs_by_frame = {}
    for t in range(corrected.shape[0]):
        ch1_frame = corrected[t, 0]
        ch2_frame = corrected[t, 1]
        blobs = detect_spots(ch1_frame, k=10, min_sigma=5, max_sigma=10, threshold=0.1)
        intensities = measure_spot_intensities(ch2_frame, blobs)
        curr_ids, next_id = link_granules(prev_blobs, prev_ids, blobs, next_id)

        for (y, x, sigma), mean_intensity, granule_id in zip(
            blobs, intensities, curr_ids
        ):
            rows.append(
                {
                    "frame": t,
                    "time": timestamps[t],
                    "granule_id": granule_id,
                    "y": y,
                    "x": x,
                    "radius": sigma * np.sqrt(2),
                    "mean_intensity_ch2": mean_intensity,
                }
            )

        blob_segmentations[t] = np.zeros(ch1_frame.shape, dtype=int)
        blob_segmentations[t][blobs[:, 0].astype(int), blobs[:, 1].astype(int)] = (
            curr_ids
        )
        blobs_by_frame[t] = blobs
        prev_blobs, prev_ids = blobs, curr_ids

    results = pd.DataFrame(rows)
    results.to_csv(RESULTS_PATH, index=False)

    print(f"bleach start: {bleach_start}, end: {bleach_end}")
    print(
        f"bleach mask: {'none' if bleach_mask is None else f'{bleach_mask.sum()} px'}"
    )
    print(f"drift shifts (first 5): {shifts[:5]}")
    print(
        f"detected {len(results)} spot-frame measurements across {corrected.shape[0]} frames"
    )
    print(f"wrote {RESULTS_PATH}")

    sns.set_style("ticks")

    plt.rcParams.update(
        {
            "font.size": 6,  # default text
            "font.family": "Arial",
            "axes.titlesize": 8,  # subplot titles
            "axes.labelsize": 6,  # x/y axis labels
            "axes.labelweight": "bold",  # x/y axis labels
            "xtick.labelsize": 6,  # x tick labels
            "ytick.labelsize": 6,  # y tick labels
            "legend.fontsize": 6,  # legend
            "figure.titlesize": 8,  # suptitle
        }
    )
    fig = plt.figure(figsize=(4.5, 2))

    ncols = 5
    nrows = 2

    gs = GridSpec(nrows, ncols, figure=fig, hspace=0.3, wspace=0.3)

    ax_unconverted0 = fig.add_subplot(gs[0, 0])
    ax_unconverted1 = fig.add_subplot(gs[0, 1])
    ax_unconverted2 = fig.add_subplot(gs[0, 2])
    ax_plot_activated = fig.add_subplot(gs[0, 3:])

    ax_converted0 = fig.add_subplot(gs[1, 0])
    ax_converted1 = fig.add_subplot(gs[1, 1])
    ax_converted2 = fig.add_subplot(gs[1, 2])
    ax_plot_bleached = fig.add_subplot(gs[1, 3:])

    last_frame = corrected.shape[0] - 1

    # Plot frame before bleach
    ax_unconverted0.imshow(corrected[frame_before_bleach, 0], cmap="gray")
    ax_unconverted0.set_title(f"Frame {frame_before_bleach} (ch1)")
    draw_blob_contours(ax_unconverted0, blobs_by_frame[frame_before_bleach])

    ax_converted0.imshow(corrected[frame_before_bleach, 1], cmap="magma")
    ax_converted0.set_title(f"Frame {frame_before_bleach} (ch2)")
    draw_blob_contours(ax_converted0, blobs_by_frame[frame_before_bleach])

    # Plot frame after bleach
    ax_unconverted1.imshow(corrected[frame_after_bleach, 0], cmap="gray")
    ax_unconverted1.set_title(f"Frame {frame_after_bleach} (ch1)")
    draw_blob_contours(ax_unconverted1, blobs_by_frame[frame_after_bleach])

    ax_converted1.imshow(corrected[frame_after_bleach, 1], cmap="magma")
    ax_converted1.set_title(f"Frame {frame_after_bleach} (ch2)")
    draw_blob_contours(ax_converted1, blobs_by_frame[frame_after_bleach])

    # Plot the last frame
    ax_unconverted2.imshow(corrected[-1, 0], cmap="gray")
    ax_unconverted2.set_title(f"Frame {last_frame} (ch1)")
    draw_blob_contours(ax_unconverted2, blobs_by_frame[last_frame])

    ax_converted2.imshow(corrected[-1, 1], cmap="magma")
    ax_converted2.set_title(f"Frame {last_frame} (ch2)")
    draw_blob_contours(ax_converted2, blobs_by_frame[last_frame])

    # Plot the mean intensity of detected spots over time and the bleach time as blue fill
    sns.lineplot(
        data=results,
        x="time",
        y="mean_intensity_ch2",
        hue="granule_id",
        ax=ax_plot_activated,
        ci=None,
    )
    plt.show()


if __name__ == "__main__":
    main()
