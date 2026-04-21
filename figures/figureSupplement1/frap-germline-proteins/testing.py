from pathlib import Path
import tifffile as tiff
from frappy import FrapSample, FrapExperiment
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

root = Path(
    "/Volumes/HELHEIM/analyzed_data/diffusivity/202207_FRAP_germ-granule_components/"
)

export_path = Path("/Users/julian/local_files/frap_test")

experiment_dict = {
    "gra_10hpf": root / "FRAP_gra" / "10hpf",
    # "vasa(AA1-164)_10hpf": root / "FRAP_hypergerm" / "10hpf",
    # "vasa_10hpf": root / "FRAP_full_vasa" / "10hpf",
}

all_timepoints = []

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
for fused_name, exp_path in experiment_dict.items():
    print(f"Processing {fused_name}...")

    exp_protein_name, exp_stage = fused_name.split("_", 2)
    exp = FrapExperiment(exp_path)
    exp.generate_report(ax=ax)
    for timepoint_data in exp.timepoint_dicts:
        timepoint_data.update(
            {
                "protein_name": exp_protein_name,
                "stage": exp_stage,
            }
        )
plt.tight_layout()
plt.show()
