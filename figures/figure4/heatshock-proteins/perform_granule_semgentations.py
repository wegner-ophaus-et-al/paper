from segmenter import SegmentationModel
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from processor import process_sample
from tqdm import tqdm
from multiprocessing import Pool
from utils.utils import statistical_analysis, get_stars, get_representative_ids
import numpy as np
import matplotlib as mpl


sm_g = SegmentationModel(
    "/Users/julian/local_files/msam_models/sam_granules_refined_up3_34917658",
    model_type="vit_l_lm",
)

root = Path("/Volumes/Kur/paper_data_sorted/heatshock_stress-granule_markers")

sample_folders = []
for construct_id in ["F617", "F618"]:
    per_construct_folder_list = [
        sf
        for sf in (root / construct_id).glob("*")
        if sf.is_dir() and not sf.name.startswith(".")
    ]
    sample_folders.extend(per_construct_folder_list)
print(f"There are a total of {len(sample_folders)} sample folders to process.")


data_records = []
fig_scaler = 1.5
fig_summary, ax_summary = plt.subplots(
    ncols=3,
    nrows=len(sample_folders),
    figsize=(3 * fig_scaler, fig_scaler * len(sample_folders)),
)
axes_indices = list(range(len(sample_folders)))
for sample_folder, ax_idx in tqdm(
    zip(sample_folders, axes_indices),
    total=len(sample_folders),
    desc="Processing samples",
):
    data_records.append(process_sample(sample_folder, sm_g, ax=ax_summary[ax_idx, :]))

for axs in ax_summary.flatten():
    axs.axis("off")
fig_summary.tight_layout()
plt.savefig(root / "granule_segmentation_summary.pdf", dpi=300, bbox_inches="tight")

df = pd.DataFrame(data_records)
df["mean_int_ratio_stress_cytoplasm_to_granule"] = (
    df["mean_int_stress_cytoplasm"] / df["mean_int_stress_granules"]
)

df.to_csv(root / "granule_segmentation_results.csv", index=False)

color_palette = {"control": "#888888", "heatshock": "#b04a46"}

sns.set_style("ticks")
# fig, ax = plt.subplots(figsize=(4, 6))
# sns.violinplot(
#     data=df,
#     y="mean_int_ratio_stress_cytoplasm_to_granule",
#     x="construct_id",
#     hue="condition",
#     hue_order=["control", "heatshock"],
#     palette=color_palette,
#     ax=ax,
#     split=True,
#     inner="quart",
# )
# sns.swarmplot(
#     data=df,
#     y="mean_int_ratio_stress_cytoplasm_to_granule",
#     x="construct_id",
#     hue="condition",
#     hue_order=["control", "heatshock"],
#     palette=["black"],
#     dodge=True,
#     ax=ax,
# )
#

plt.rcParams.update(
    {
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)


mpl.rcParams["font.family"] = "Arial"

fig_mean_cyto_gran_newstyle, ax_mean_cyto_gran_newstyle = plt.subplots(
    1, 1, figsize=(1.63, 2.5)
)

ax = ax_mean_cyto_gran_newstyle

sns.stripplot(
    data=df,
    y="mean_int_ratio_stress_cytoplasm_to_granule",
    x="construct_id",
    order=["F618", "F617"],
    hue="condition",
    hue_order=["control", "heatshock"],
    palette=color_palette,
    dodge=True,
    alpha=0.4,
    size=3,
    ax=ax,
    jitter=0.3,
)
sns.pointplot(
    data=df,
    y="mean_int_ratio_stress_cytoplasm_to_granule",
    x="construct_id",
    order=["F618", "F617"],
    hue="condition",
    hue_order=["control", "heatshock"],
    dodge=0.4,
    errorbar="sd",
    estimator="median",
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


sns.despine()
plt.tight_layout()
plt.savefig(root / "cytoplasm_granule_stress_pointplot.pdf", bbox_inches="tight")
plt.close()

feature_to_get_stats_from = [
    "mean_int_ratio_stress_cytoplasm_to_granule",
    "total_cell_area",
    "total_granules_area",
]
all_stats_results = []
for construct_id in df["construct_id"].unique():
    for feature in feature_to_get_stats_from:
        df_construct = df[df["construct_id"] == construct_id]
        ctrl_vals = df_construct[df_construct["condition"] == "control"][
            feature
        ].to_list()
        hs_vals = df_construct[df_construct["condition"] == "heatshock"][
            feature
        ].to_list()

        test_type, _, p_val = statistical_analysis(ctrl_vals, hs_vals)
        stat_result = {
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
        all_stats_results.append(stat_result)

        dict_representative_ids = get_representative_ids(
            df.query("construct_id == @construct_id"), feature
        )

        with open(root / f"representative_ids_{construct_id}_{feature}.txt", "w") as f:
            for condition, ids in dict_representative_ids.items():
                f.write(f"{condition}:\n")
                for uid in ids:
                    f.write(f"{uid}\n")
                f.write("\n")

df_stat_results = pd.DataFrame(all_stats_results)
df_stat_results.to_csv(
    root / "granule_segmentation_statistical_results.csv", index=False
)
