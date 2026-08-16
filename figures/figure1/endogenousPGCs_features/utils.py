import numpy as np
from scipy import stats
import pandas as pd
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from skimage import measure


def parametric(data_set1, data_set2):
    set1_normal = stats.normaltest(data_set1).pvalue > 0.05
    set2_normal = stats.normaltest(data_set2).pvalue > 0.05

    return min(set1_normal, set2_normal)


def statistical_analysis(data_set1: list, data_set2: list):
    if parametric(data_set1, data_set2):
        t_statistic, p_value = stats.ttest_ind(data_set1, data_set2)
        test_type = "t-test"
    else:
        t_statistic, p_value = stats.mannwhitneyu(data_set1, data_set2)
        test_type = "Mann-Whitney U test"
    return test_type, t_statistic, p_value


def get_stars(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def print_nested(d, indent=0):
    for key, value in d.items():
        prefix = "  " * indent
        if isinstance(value, dict):
            print(f"{prefix}{key}:")
            print_nested(value, indent + 1)
        elif (
            isinstance(value, list)
            and value
            and all(isinstance(i, dict) for i in value)
        ):
            print(f"{prefix}{key}: list ({len(value)})")
            print_nested(value[0], indent + 1)
        elif isinstance(value, np.ndarray):
            print(f"{prefix}{key}: ndarray {value.shape}")
        else:
            print(f"{prefix}{key}: {type(value).__name__}")


def export_representative_cells(
    df: pd.DataFrame,
    feature: str,
    export_path,
    percentile_low: float = 40,
    percentile_high: float = 60,
):
    """
    Export uids whose per-cell feature values are representative for each stage,
    condition, and marker.

    The feature is aggregated per (uid, stage, condition, marker) before
    computing percentiles. A uid is only exported for a given stage/condition
    when it is representative for every marker present in that group.
    """
    required_columns = {"uid", "stage", "condition", "marker", feature}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required columns: {missing}")

    if percentile_low > percentile_high:
        raise ValueError("percentile_low must be <= percentile_high")

    export_dir = Path(export_path)
    if export_dir.suffix:
        export_dir = export_dir.parent
    export_dir.mkdir(parents=True, exist_ok=True)
    export_file = export_dir / f"Representative_{feature}.txt"

    working = df.loc[:, ["uid", "stage", "condition", "marker", feature]].copy()
    working[feature] = pd.to_numeric(working[feature], errors="coerce")
    working = working.dropna(subset=["uid", "stage", "condition", "marker", feature])

    if working.empty:
        export_file.write_text(
            "\n".join(
                [
                    f"Representative cells for feature: {feature}",
                    f"Percentile range: {percentile_low}-{percentile_high}",
                    "No valid rows were available.",
                    "",
                ]
            )
        )
        return export_file

    per_cell = (
        working.groupby(["uid", "stage", "condition", "marker"], as_index=False)[
            feature
        ]
        .mean()
        .rename(columns={feature: "feature_value"})
    )

    lines = [
        f"Representative cells for feature: {feature}",
        f"Percentile range: {percentile_low}-{percentile_high}",
        "Representative means are computed per uid, stage, condition, and marker.",
        "A uid is selected only if it falls within the percentile window for all markers in the group.",
        "",
    ]

    for (stage, condition), group in per_cell.groupby(
        ["stage", "condition"], sort=True
    ):
        lines.append(f"Stage: {stage}")
        lines.append(f"Condition: {condition}")

        marker_frames = []
        marker_summaries = []
        marker_order = list(dict.fromkeys(group["marker"].tolist()))

        for marker in marker_order:
            marker_group = group[group["marker"] == marker].copy()
            low = marker_group["feature_value"].quantile(percentile_low / 100.0)
            high = marker_group["feature_value"].quantile(percentile_high / 100.0)
            representative = marker_group[
                (marker_group["feature_value"] >= low)
                & (marker_group["feature_value"] <= high)
            ].copy()
            marker_frames.append(
                representative[["uid", "feature_value"]].rename(
                    columns={"feature_value": f"{marker}_value"}
                )
            )
            marker_summaries.append(
                f"  marker={marker} q{percentile_low}={low:.6g} q{percentile_high}={high:.6g} selected={len(representative)}"
            )

        common_uids = set(marker_frames[0]["uid"]) if marker_frames else set()
        for marker_frame in marker_frames[1:]:
            common_uids &= set(marker_frame["uid"])

        lines.extend(marker_summaries)
        lines.append(f"  shared_representative_uids={len(common_uids)}")

        if common_uids:
            marker_value_tables = {
                marker: frame.set_index("uid")[f"{marker}_value"]
                for marker, frame in zip(marker_order, marker_frames)
            }
            for uid in sorted(common_uids, key=str):
                uid_values = ", ".join(
                    f"{marker}={marker_value_tables[marker].loc[uid]:.6g}"
                    for marker in marker_order
                )
                lines.append(f"    uid={uid} {uid_values}")

        lines.append("")

    export_file.write_text("\n".join(lines).rstrip() + "\n")
    return export_file


def _normalize(img):
    """Scale an intensity image to the [0, 1] range."""
    img = img.astype(float)
    mn, mx = np.nanmin(img), np.nanmax(img)
    if mx > mn:
        return (img - mn) / (mx - mn)
    return np.zeros_like(img)


def _plot_contours(ax, seg, color="white", lw=0.8):
    """Draw the outline of every labelled object in `seg`."""
    for label in np.unique(seg):
        if label == 0:  # skip background
            continue
        for contour in measure.find_contours(seg == label, 0.5):
            ax.plot(contour[:, 1], contour[:, 0], color=color, linewidth=lw)


def _plot_scalebar(
    ax,
    img_shape,
    bar_length_um=5,
    sensor_pixel_size=6.5,
    binning=1,
    magnification=63,
    color="white",
    height_frac=0.015,
    pad_frac=0.04,
):
    """
    Draw a scale bar on `ax`.

    pixel_size = (sensor_pixel_size * binning) / magnification   [µm/pixel]
    """
    pixel_size = (sensor_pixel_size * binning) / magnification  # µm per pixel
    bar_length_px = bar_length_um / pixel_size  # length in pixels

    n_rows, n_cols = img_shape[:2]

    # geometry of the bar
    bar_height = n_rows * height_frac
    pad_x = n_cols * pad_frac
    pad_y = n_rows * pad_frac

    x0 = n_cols - pad_x - bar_length_px
    y0 = n_rows - pad_y - bar_height

    ax.add_patch(
        Rectangle((x0, y0), bar_length_px, bar_height, color=color, edgecolor=None)
    )
    ax.text(
        x0 + bar_length_px / 2,
        y0 - bar_height,
        f"{bar_length_um} µm",
        color=color,
        ha="center",
        va="bottom",
        fontsize=8,
    )


def pub_images(ax, thecell, **kwargs):

    dnd_img = thecell.markers["dnd1"].raw_image
    gra_img = thecell.markers["gra"].raw_image
    nls_img = thecell.markers["cell"].raw_image
    cell_seg = thecell.markers["cell"].segmentation
    name = thecell.uid
    condition = thecell.condition
    stage = thecell.stage

    # --- linear black->color colormaps -------------------------------------
    cyan_cmap = LinearSegmentedColormap.from_list("cyan", [(0, 0, 0), (0, 1, 1)])
    yellow_cmap = LinearSegmentedColormap.from_list("yellow", [(0, 0, 0), (1, 1, 0)])

    # --- normalize each channel to [0, 1] ----------------------------------
    gra_n = _normalize(gra_img)  # cyan
    dnd_n = _normalize(dnd_img)  # yellow
    nls_n = _normalize(nls_img)  # magenta

    # --- panel 1: GRA in cyan ----------------------------------------------
    ax[0].imshow(gra_n, cmap=cyan_cmap, vmin=0, vmax=1)
    ax[0].set_title(f"{name} – {condition} - {stage}")

    # --- panel 2: DND in yellow --------------------------------------------
    ax[1].imshow(dnd_n, cmap=yellow_cmap, vmin=0, vmax=1)

    # --- panel 3: color-accurate additive merge ----------------------------
    # cyan    = (0, 1, 1)  -> contributes to G, B
    # yellow  = (1, 1, 0)  -> contributes to R, G
    # magenta = (1, 0, 1)  -> contributes to R, B
    cyan, yellow, magenta = gra_n, dnd_n, nls_n
    R = np.clip(yellow + magenta, 0, 1)
    G = np.clip(cyan + yellow, 0, 1)
    B = np.clip(cyan + magenta, 0, 1)
    merge = np.dstack([R, G, B])
    ax[2].imshow(merge)
    ax[2].set_title(f"{name} – merge")

    # --- overlay segmentation contours on every panel ----------------------
    for a in ax:
        _plot_contours(a, cell_seg)
        _plot_scalebar(a, cell_seg.shape, **kwargs)
        a.axis("off")
        # a.set_xticks([])
        # a.set_yticks([])

    return ax
