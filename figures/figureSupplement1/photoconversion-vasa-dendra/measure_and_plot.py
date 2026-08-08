import tifffile as tiff
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from skimage.measure import regionprops, label
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


COLORS = {
    "cyan": (0, 1, 1),
    "magenta": (1, 0, 1),
    "yellow": (1, 1, 0),
    "grey": (1, 1, 1),
}


def _cmap(name):
    return LinearSegmentedColormap.from_list(name, [(0, 0, 0), COLORS[name]])


img = tiff.imread("data/p11_directHit_bigGranule_same_cell_as_p10_dc_imagej.tif")
img_uc = img[:, 0]
img_conv = img[:, 1]

mask_path = "data/masks/granule_drift_corrected.tif"
mask = label(tiff.imread(mask_path))

timepoints = np.asarray(np.arange(img.shape[0]))
start_bleach = 1
end_bleach = 2.4
timepoints = timepoints - end_bleach


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


fig = plt.figure(figsize=(4.5, 2))

ncols = 5
nrows = 2

gs = GridSpec(nrows, ncols, figure=fig, hspace=0.3, wspace=0.3)

ax_plot_activated = fig.add_subplot(gs[0, 3:])
ax_plot_recovered = fig.add_subplot(gs[1, 3:])

sns.lineplot(
    df.query("label_index == 3"),
    x="time",
    y="mean_intensity_converted",
    hue="label_index",
    ax=ax_plot_activated,
    palette=palette,
)

sns.lineplot(
    df.query("label_index != 3"),
    x="time",
    y="mean_intensity_converted",
    hue="label_index",
    ax=ax_plot_recovered,
    palette=palette,
)

for show_int, show_frame in enumerate([2, 3, -1]):
    ax_img_uc = fig.add_subplot(gs[0, show_int])
    ax_img_uc.imshow(img_uc[show_frame], cmap=_cmap("cyan"))

    ax_img_conv = fig.add_subplot(gs[1, show_int])
    ax_img_conv.imshow(img_conv[show_frame], vmin=0, vmax=5000, cmap="inferno")
    for lbl_idx, edge_color in palette.items():
        ax_img_conv.contour(mask == lbl_idx, colors=edge_color, linestyles="dashed")
    ax_img_uc.set_axis_off()
    ax_img_conv.axis("off")

    # ax_img_conv.

sns.despine()
plt.tight_layout()
fig.savefig("out.pdf")
plt.show()
