import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import measure

from filtering import segment_spots


def detect_and_evaluate(
    pla_image,
    cell_mask,
    granule_mask,
    uid,
    condition,
    sample_id,
    rna_type,
    sigma=1.5,
    k=5.0,
    min_area=4,
    pixel_size=1.0,
):
    """Detect PLA puncta within cell_mask, extract per-cell summary + per-spot records."""
    labels, _, _ = segment_spots(pla_image, sigma=sigma, k=k, min_area=min_area)
    labels = measure.label(labels * (cell_mask > 0))

    # granule_distance_map = ndi.distance_transform_edt(granule_mask == 0) * pixel_size
    outside_distance = ndi.distance_transform_edt(granule_mask == 0) * pixel_size
    inside_distance = ndi.distance_transform_edt(granule_mask > 0) * pixel_size

    granule_distance_map = outside_distance.copy()
    granule_distance_map[outside_distance == 0] = inside_distance[outside_distance == 0]

    base = {
        "uid": uid,
        "condition": condition,
        "sample_id": sample_id,
        "rna_type": rna_type,
    }

    summary = base.copy()
    summary.update(
        {
            "number_of_spots": labels.max(),
            "sum_intensity": np.sum(pla_image[labels > 0]),
            "mean_intensity": np.mean(pla_image[labels > 0])
            if np.sum(labels > 0) > 0
            else 0,
            "area_of_spots": np.sum(labels > 0),
        }
    )

    details = []
    for region in measure.regionprops(labels, intensity_image=pla_image):
        centroid = tuple(np.round(region.centroid).astype(int))
        detail = base.copy()
        detail.update(
            {
                "spot_label": region.label,
                "centroid": region.centroid,
                "area": region.area,
                "intensity_sum": region.image_intensity.sum(),
                "intensity_mean": region.image_intensity.mean(),
                "distance_to_granule": granule_distance_map[centroid],
            }
        )
        details.append(detail)

    return labels, summary, details


def _strip_point(ax, data, y, color_palette):
    sns.stripplot(
        data=data,
        y=y,
        x="condition",
        hue="condition",
        palette=color_palette,
        alpha=0.6,
        size=3,
        ax=ax,
    )
    sns.pointplot(
        data=data,
        y=y,
        x="condition",
        hue="condition",
        errorbar="sd",
        estimator="mean",
        capsize=0.075,
        linestyle="none",
        markersize=10,
        marker="_",
        err_kws=dict(linewidth=0.4, color="black"),
        markeredgewidth=1,
        palette="dark:black",
        zorder=5,
        ax=ax,
    )


def plot_spot_summary(df_cells, color_palette):
    """2-panel strip+point figure: number_of_spots and sum_intensity by condition."""
    fig, axes = plt.subplots(1, 2, figsize=(3, 2))

    _strip_point(axes[0], df_cells, "number_of_spots", color_palette)
    axes[0].set_title("PLA spots per cell")
    axes[0].set_ylabel("Number of spots")
    axes[0].set_xlabel("Condition")

    _strip_point(axes[1], df_cells, "mean_intensity", color_palette)
    axes[1].set_title("PLA spot intensity")
    axes[1].set_ylabel("Summed intensity in spots (a.u.)")
    axes[1].set_xlabel("Condition")

    sns.despine()
    fig.tight_layout()
    return fig


def plot_distance_kde(df_spots, color_palette, control_condition="DMSO"):
    """KDE of spot->nearest-granule distance for control-condition spots, split by rna_type."""
    control = df_spots[df_spots["condition"] == control_condition]
    rna_types = ["actin", "vasa"]

    fig, axes = plt.subplots(len(rna_types), 1, figsize=(3, 3))

    for ax, rna_type in zip(axes, rna_types):
        subset = control[control["rna_type"] == rna_type]
        vals = subset["distance_to_granule"].dropna()

        sns.kdeplot(
            data=subset,
            x="distance_to_granule",
            hue="condition",
            ax=ax,
            fill=True,
            alpha=0.5,
            palette=color_palette,
        )
        ax.axvline(
            vals.mean(), color="k", ls="--", lw=1.2, label=f"Mean = {vals.mean():.2f}"
        )
        ax.axvline(
            vals.median(),
            color="k",
            ls=":",
            lw=1.2,
            label=f"Median = {vals.median():.2f}",
        )

        ax.set_xlabel("Distance to nearest granule ($\\alpha m)")
        ax.set_title(f"{rna_type} ({control_condition})")
        ax.legend()

    sns.despine()
    fig.tight_layout()
    return fig


def plot_distance_vs_intensity_hex(df_spots, control_condition="DMSO"):
    """Hexbin joint plot of spot distance-to-granule vs. summed intensity (control only)."""
    control = df_spots[df_spots["condition"] == control_condition]

    g = sns.jointplot(
        x="distance_to_granule",
        y="intensity_mean",
        data=control,
        kind="hex",
        height=4,
    )
    g.set_axis_labels(
        "Distance to nearest granule ($\\micro m)", "mean spot intensity (a.u.)"
    )
    return g
