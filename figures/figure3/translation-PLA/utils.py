from __future__ import annotations

import numpy as np
from skimage import draw, measure
from scipy import stats


def _cross(origin: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (
        b[0] - origin[0]
    )


def _graham_scan(points: np.ndarray) -> np.ndarray:
    if len(points) <= 1:
        return points

    points = np.array(sorted(points.tolist()))

    lower = []
    for point in points:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)

    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)

    return np.array(lower[:-1] + upper[:-1])


def get_pseudo_nucleus_mask(label_image, image, cell_mask):
    """
    Build a convex hull around all labeled centroids and subtract it from the cell mask.
    """
    if label_image.shape != image.shape or cell_mask.shape != image.shape:
        raise ValueError("label_image, image, and cell_mask must share the same shape")

    centroids = np.array(
        [
            (prop.centroid[1], prop.centroid[0])
            for prop in measure.regionprops(label_image)
        ]
    )

    if len(centroids) == 0:
        return cell_mask.astype(bool)

    hull = _graham_scan(centroids)
    hull_mask = np.zeros(cell_mask.shape, dtype=bool)

    if len(hull) == 1:
        col = int(np.clip(np.round(hull[0][0]), 0, cell_mask.shape[1] - 1))
        row = int(np.clip(np.round(hull[0][1]), 0, cell_mask.shape[0] - 1))
        hull_mask[row, col] = True
    elif len(hull) == 2:
        rr, cc = draw.line(
            int(np.clip(np.round(hull[0][1]), 0, cell_mask.shape[0] - 1)),
            int(np.clip(np.round(hull[0][0]), 0, cell_mask.shape[1] - 1)),
            int(np.clip(np.round(hull[1][1]), 0, cell_mask.shape[0] - 1)),
            int(np.clip(np.round(hull[1][0]), 0, cell_mask.shape[1] - 1)),
        )
        hull_mask[rr, cc] = True
    else:
        rr, cc = draw.polygon(hull[:, 1], hull[:, 0], shape=cell_mask.shape)
        hull_mask[rr, cc] = True

    return hull_mask


def get_measurements(func, image, segmentation):
    return func(image[segmentation > 0])


def calculate_area_ratio(df, image_name, seg1, seg2, measurement="sum"):
    """Calculate ratio of PLA signal between two segmentation areas"""
    subset = df[
        (df["image_name"] == image_name) & (df["measurement_name"] == measurement)
    ].copy()

    pivot = subset.pivot_table(
        index=["uid", "condition", "sample_id", "rna_type"],
        columns="segmentation_name",
        values="value",
    )

    pivot["ratio"] = pivot[seg1] / pivot[seg2]
    return pivot


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


def statistics(df_pivot, rna_type, ratio_desc):
    """
    Perform statistical analysis on pivoted ratio data.

    Args:
        df_pivot: Pivoted dataframe with ratio column
        rna_type: RNA type being analyzed (e.g., 'vasa' or 'actin')
        ratio_desc: Description of ratio (e.g., 'granule / granule_periphery')

    Returns:
        List of formatted strings with statistical results
    """
    df_pivot = df_pivot.copy()  # Avoid modifying the original dataframe
    df_pivot.reset_index(inplace=True)  # Ensure 'condition' is a column for filtering
    results = []
    results.append(f"{rna_type} mRNA probe:")
    results.append(f"    Ratio: {ratio_desc}")

    for drug in ["CHX", "PatA"]:
        ctrl_vals = df_pivot[df_pivot["condition"] == "DMSO"]["ratio"].to_list()
        drug_vals = df_pivot[df_pivot["condition"] == drug]["ratio"].to_list()

        if len(ctrl_vals) == 0 or len(drug_vals) == 0:
            continue

        ctrl_mean = np.mean(ctrl_vals)
        drug_mean = np.mean(drug_vals)
        ctrl_count = len(ctrl_vals)
        drug_count = len(drug_vals)

        test_name, _, p_value = statistical_analysis(ctrl_vals, drug_vals)
        stars = get_stars(p_value)

        results.append(f"    DMSO-{drug}:")
        results.append(f"        DMSO mean:    {ctrl_mean:.5f}")
        results.append(f"        {drug} mean:    {drug_mean:.5f}")
        results.append(f"        DMSO count:   {ctrl_count}")
        results.append(f"        {drug} count:   {drug_count}")
        results.append(f"        Test name:    {test_name}")
        results.append(f"        p value:      {p_value:.7f}")
        results.append(f"        stars:        {stars}")

    return results
