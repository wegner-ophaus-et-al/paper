import tifffile as tiff
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from skimage.measure import label
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from utils import scale_bar, lsm_pixel_size, bleach_times, get_lsm_meta


COLORS = {
    "cyan": (0, 1, 1),
    "magenta": (1, 0, 1),
    "yellow": (1, 1, 0),
    "grey": (1, 1, 1),
}


def _cmap(name):
    return LinearSegmentedColormap.from_list(name, [(0, 0, 0), COLORS[name]])


img = tiff.imread("data/p11_directHit_bigGranule_same_cell_as_p10_dc_imagej.tif")
lsm_path = "data/p11_directHit_bigGranule_same_cell_as_p10.lsm"
img_uc = img[:, 0]
img_conv = img[:, 1]

mask_path = "data/masks/granule_drift_corrected.tif"
mask = label(tiff.imread(mask_path))


lsm_meta = get_lsm_meta(lsm_path)
timestamps = lsm_meta["TimeStamps"]
bleach_start, bleach_end = bleach_times(lsm_meta)
frame_before_bleach = np.searchsorted(lsm_meta["TimeStamps"], bleach_start) - 1
frame_after_bleach = np.searchsorted(lsm_meta["TimeStamps"], bleach_end)

timepoints = np.array(timestamps) - bleach_end
bleach_start = bleach_start - bleach_end
bleach_end = bleach_end - bleach_end  # This will be 0, but keeping for clarity


print(f"Bleach start: {bleach_start:.3f}s, bleach end: {bleach_end:.3f}s")
print(
    f"Frame before bleach: {frame_before_bleach}, frame after bleach: {frame_after_bleach}"
)

results = []
for tp, (uc_frame, conv_frame) in enumerate(zip(img_uc, img_conv)):
    for label_index in range(mask.max()):
        if label_index == 0:
            continue
        results.append(
            {
                "frame": tp,
                "time": timepoints[tp],
                "label_index": label_index,
                "sum_intensity_unconverted": np.sum(uc_frame[mask == label_index]),
                "sum_intensity_converted": np.sum(conv_frame[mask == label_index]),
                "mean_intensity_unconverted": np.mean(uc_frame[mask == label_index]),
                "mean_intensity_converted": np.mean(conv_frame[mask == label_index]),
            }
        )

df = pd.DataFrame(results)

sns.set_theme(style="ticks")
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "axes.labelweight": "bold",  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)

palette = {
    3: "#1F80DB",  # Activated granule
    2: "#DB7B1F",
    1: "#E7A76A",
    4: "#F3D3B4",
}


ncols = 5
col_width = 1.0

nrows = 2
row_height = 1.5


fig = plt.figure(figsize=(col_width * ncols, nrows * row_height))
gs = GridSpec(nrows, ncols, figure=fig)  # , hspace=0.3, wspace=0.3)

ax_plot_activated = fig.add_subplot(gs[0, 3:])
ax_plot_recovered = fig.add_subplot(gs[1, 3:], sharex=ax_plot_activated)

sns.lineplot(
    df.query("label_index == 3"),
    x="time",
    y="mean_intensity_converted",
    hue="label_index",
    ax=ax_plot_activated,
    palette=palette,
)
ax_plot_activated.axvspan(bleach_start, bleach_end, color="#AE00FF", alpha=0.3)
ax_plot_activated.legend_.remove()
ax_plot_activated.set_title("Photoconverted granule")
ax_plot_activated.set_xlabel("Time (s)")
ax_plot_activated.set_ylabel("Mean intensity (a.u.)")
ax_plot_activated.set_xlabel("")
ax_plot_activated.tick_params(labelbottom=False)


sns.lineplot(
    df.query("label_index != 3"),
    x="time",
    y="mean_intensity_converted",
    hue="label_index",
    ax=ax_plot_recovered,
    palette=palette,
)
ax_plot_recovered.axvspan(bleach_start, bleach_end, color="#AE00FF", alpha=0.3)
ax_plot_recovered.legend_.remove()
ax_plot_recovered.set_title("Recovered granules")
ax_plot_recovered.set_xlabel("Time (s)")
ax_plot_recovered.set_ylabel("Mean intensity (a.u.)")

gs.tight_layout(fig)
for show_int, show_frame in enumerate([2, 3, -1]):
    current_time = timepoints[show_frame]

    ax_img_uc = fig.add_subplot(gs[0, show_int])
    ax_img_uc.imshow(img_uc[show_frame], cmap=_cmap("cyan"))
    scale_bar(
        ax_img_uc,
        length=5,
        thickness=0.05,
        pixel_size=lsm_pixel_size(lsm_path),
        pad=0.03,
        color="white",
    )
    ax_img_uc.text(
        0,
        0,
        f"t = {current_time:.1f}s",
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=ax_img_uc.transAxes,
        color="white",
        fontsize=6,
    )
    ax_img_uc.set_axis_off()

    ax_img_conv = fig.add_subplot(gs[1, show_int])
    ax_img_conv.imshow(img_conv[show_frame], vmin=0, vmax=5000, cmap="inferno")
    scale_bar(
        ax_img_conv,
        length=5,
        thickness=0.05,
        pixel_size=lsm_pixel_size(lsm_path),
        pad=0.03,
        color="white",
    )

    for lbl_idx, edge_color in palette.items():
        ax_img_conv.contour(
            mask == lbl_idx, colors=edge_color, linestyles="dotted", linewidths=0.5
        )
    ax_img_conv.axis("off")

    ax_img_conv.text(
        0,
        0,
        f"t = {current_time:.1f}s",
        verticalalignment="bottom",
        horizontalalignment="left",
        transform=ax_img_conv.transAxes,
        color="white",
        fontsize=6,
    )
    # ax_img_conv.

sns.despine()
fig.savefig("out.pdf", dpi=1200, transparent=True, bbox_inches="tight")
plt.show()
