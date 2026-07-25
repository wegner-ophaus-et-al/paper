import json
from pathlib import Path
from frappy import FrapExperiment
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

root = Path(
    "/Volumes/HELHEIM/analyzed_data/diffusivity/202207_FRAP_germ-granule_components/"
)

export_path = Path("/Users/julian/local_files/frap_test")

experiment_dict = {
    "gra_10hpf": root / "FRAP_gra" / "10hpf",
    "vasa(AA1-164)_10hpf": root / "FRAP_hypergerm" / "10hpf",
    "vasa_10hpf": root / "FRAP_full_vasa" / "10hpf",
    "dnd_10hpf": root / "FRAP_dnd" / "10hpf",
    "nanos_10hpf": root / "FRAP_nanos" / "10hpf",
}


molecular_weights = {
    "gra": 30.0,  # kDa
    "vasa": 65.0,  # kDa
    "vasa(AA1-164)": 18.0,  # kDa
    "dnd": 65.0,  # kDa
    "nanos": 38.0,  # kDa
}

palette = {
    "gra": "#2A5F6E",  # muted dark teal
    "vasa": "#4E9AA6",  # muted teal
    "vasa(AA1-164)": "#7FB3C4",  # muted light blue
    "dnd": "#C9922B",  # muted gold
    "nanos": "#8A6D1E",  # muted bronze
}


all_timepoints = []
fit_data = []
fit_summary = {}
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,  # default text
        "axes.titlesize": 8,  # subplot titles
        "axes.labelsize": 6,  # x/y axis labels
        "xtick.labelsize": 6,  # x tick labels
        "ytick.labelsize": 6,  # y tick labels
        "legend.fontsize": 6,  # legend
        "figure.titlesize": 8,  # suptitle
    }
)
fig, ax = plt.subplots(1, 2, figsize=(9, 2.4))
for fused_name, exp_path in experiment_dict.items():
    print(f"Processing {fused_name}...")

    exp_protein_name, exp_stage = fused_name.split("_", 2)
    exp = FrapExperiment(exp_path)
    exp.generate_report(ax=ax[0])

    tau = exp.fit_params["tau"]
    fit_summary[fused_name] = {
        "t_half": tau * np.log(2),
        "tau": tau,
        "Iinf": exp.fit_params["Iinf"],
        "r_squared": exp.fit_params["r_squared"],
    }

    for timepoint_data in exp.timepoint_dicts:
        timepoint_data.update(
            {
                "protein_name": exp_protein_name,
                "stage": exp_stage,
            }
        )
    all_timepoints.extend(exp.timepoint_dicts)


with open(export_path / "fit_summary.txt", "w") as f:
    json.dump(fit_summary, f, indent=4)

print(all_timepoints[0])
df = pd.DataFrame(all_timepoints)
print(df["protein_name"].unique())
sns.lineplot(
    data=df,
    x="index",
    y="normalized_intensity",
    hue="protein_name",
    ax=ax[1],
    errorbar="se",
    palette=palette,
    lw=0.75,
)
ax[1].set_xlim(-5, 180)
ax[1].set_ylim(-0.1, 1.25)
ax[1].set_xlabel("Frame", fontweight="bold")
ax[1].set_ylabel("Normalized intensity", fontweight="bold")
ax[1].set_title("Diffsitivity of germ line proteins", fontweight="bold")
plt.tight_layout()
sns.despine()
fig.savefig(export_path / "frap_report.pdf")
# plt.show()

statistics = (
    df.groupby(by=["protein_name", "index"]).count()["normalized_intensity"].copy()
)
statistics.to_csv(export_path / "statistics.csv")
print(statistics.groupby(by=["protein_name"]).max())
