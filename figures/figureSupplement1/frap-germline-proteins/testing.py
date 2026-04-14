from pathlib import Path
import tifffile as tiff
from frappy import FrapSample, FrapExperiment
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

root = Path(
    "/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_gra/10hpf"
)

list_of_sample_paths = [
    p for p in list(root.glob("*")) if not p.name.startswith(".") and p.is_dir()
]

fe = FrapExperiment(root)
df = fe.process_sample_data()  # Process all samples in the experiment
fe.generate_report()  # Generate a report for the entire experiment

raise ValueError("Stop here for testing")
number_of_samples = len(list_of_sample_paths)


fig = plt.figure(figsize=(21, 3 * number_of_samples))
gs = gridspec.GridSpec(
    number_of_samples, 7, figure=fig, height_ratios=[1] * number_of_samples
)
for i, sample_path in enumerate(list_of_sample_paths):
    try:
        fs = FrapSample(sample_path)
        fs.process_data()
        print(fs.metadata["time"].index)
        print(type(fs.metadata["time"].index))
        break
        ax0 = fig.add_subplot(gs[i, :2])
        ax1 = fig.add_subplot(gs[i, 2:4])
        ax2 = fig.add_subplot(gs[i, 4])
        ax3 = fig.add_subplot(gs[i, 5])
        ax4 = fig.add_subplot(gs[i, 6])

        fs.generate_report(axes=[ax0, ax1, ax2, ax3, ax4])
    except Exception as e:
        print(f"Error processing sample at {sample_path}: {e}")

sns.despine()
# plt.tight_layout()
# plt.savefig(root / "overview.pdf")
