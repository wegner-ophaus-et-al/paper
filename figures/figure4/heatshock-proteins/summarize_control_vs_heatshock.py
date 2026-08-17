from pathlib import Path

import numpy as np
import pandas as pd

from utils.utils import get_stars, statistical_analysis

root = Path("/Volumes/Kur/paper_data_sorted/heatshock_stress-granule_markers/figures")

df = pd.read_csv(root / "granule_segmentation_results.csv")

# --- Counts per RNA type (construct) and condition -------------------------
counts = df.groupby(["construct_id", "condition"]).size().unstack(fill_value=0)
counts["total"] = counts.sum(axis=1)
counts = counts.reset_index().rename_axis(None, axis=1)
counts.to_csv(root / "cell_counts_per_construct_condition.csv", index=False)

# --- Control vs heatshock significance for every numeric feature -----------
non_feature_columns = ["uid", "sample_name", "construct_id", "condition"]
feature_columns = [
    col for col in df.select_dtypes(include=np.number).columns if col not in non_feature_columns
]

all_stats_results = []
for construct_id in df["construct_id"].unique():
    df_construct = df[df["construct_id"] == construct_id]
    for feature in feature_columns:
        ctrl_vals = df_construct[df_construct["condition"] == "control"][feature].to_list()
        hs_vals = df_construct[df_construct["condition"] == "heatshock"][feature].to_list()

        test_type, _, p_val = statistical_analysis(ctrl_vals, hs_vals)
        all_stats_results.append(
            {
                "Construct": construct_id,
                "Feature": feature,
                "P-value": p_val,
                "Test type": test_type,
                "significance": get_stars(p_val),
                "Mean Control": np.mean(ctrl_vals),
                "Mean Heatshock": np.mean(hs_vals),
                "Percentage Change": (
                    (np.mean(hs_vals) - np.mean(ctrl_vals)) / np.mean(ctrl_vals) * 100
                ),
                "Count Control": len(ctrl_vals),
                "Count Heatshock": len(hs_vals),
            }
        )

df_stat_results = pd.DataFrame(all_stats_results)
df_stat_results.to_csv(
    root / "granule_segmentation_statistical_results_all_features.csv", index=False
)

print("Counts per construct/condition:")
print(counts)
print("\nWrote:")
print(f"  {root / 'cell_counts_per_construct_condition.csv'}")
print(f"  {root / 'granule_segmentation_statistical_results_all_features.csv'}")
