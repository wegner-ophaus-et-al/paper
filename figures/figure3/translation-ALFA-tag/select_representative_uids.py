import pandas as pd
import numpy as np
from pathlib import Path

output_dir = Path(__file__).parent / "output"
df = pd.read_csv(output_dir / "cell_features.csv")

all_features = [
    "nucleus_target_expression",
    "mean_alfa_cytoplasm",
    "sum_alfa_cytoplasm",
    "cv_alfa_cytoplasm",
    "mean_alfa_granules",
    "sum_alfa_granules",
    "cv_alfa_granules",
    "mean_alfa_granule_periphery",
    "sum_alfa_granule_periphery",
    "cv_alfa_granule_periphery",
    "number_of_spots",
    "area_of_spots",
    "sum_intensity_alfa_spots",
    "mean_distance_to_granules",
    "MeanRNA_Cytoplasm_Granules_Ratio",
    "SumRNA_Cytoplasm_Granules_Ratio",
]

conditions = ["Control", "PatA"]

lines = []
for feature in all_features:
    lines.append(f"# {feature}")
    for condition in conditions:
        subset = df.loc[df["condition"] == condition, ["uid", feature]].dropna()
        values = subset[feature].to_numpy()
        low, high = np.percentile(values, [40, 60])
        representative_uids = subset.loc[
            (subset[feature] >= low) & (subset[feature] <= high), "uid"
        ].tolist()

        lines.append(f"## {condition} (40th-60th pct: {low:.4g} - {high:.4g})")
        lines.extend(str(uid) for uid in representative_uids)

        print(
            f"{feature} | {condition}: {len(representative_uids)}/{len(subset)} "
            f"cells in [{low:.4g}, {high:.4g}]"
        )
    lines.append("")

with open(output_dir / "representative_uids.txt", "w") as f:
    f.write("\n".join(lines).rstrip() + "\n")
