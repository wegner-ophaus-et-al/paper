import sys

sys.path.insert(0, "/Users/icb_remote/Documents/JW/py/packages/")

from gamgee.segmenter import SegmentationModel

from pathlib import Path
from cell import Cell
import matplotlib.pyplot as plt
from tqdm import tqdm

data_root = Path(__file__).parent / "data"

cells = []
for lsm_path in data_root.rglob("*.lsm"):
    if not lsm_path.name.startswith("."):
        cells.append(Cell(lsm_path.parent.parent))

sm_granule = SegmentationModel(
    "/Users/icb_remote/Documents/JW/py/packages/gamgee/models/msam/sam_granules_refined_up3_34917658",
    model_type="vit_l_lm",
)


ncols = 3
nrows = len(cells)
width_size = 1
height_size = 1
fig, axes = plt.subplots(
    nrows=len(cells), ncols=3, figsize=(width_size * ncols, height_size * nrows)
)

with tqdm(total=len(cells), desc="Processing cells") as pbar:
    for ax, cell in zip(axes, cells):
        cell.segment(sm_granule)  # Uncomment this line to perform segmentation
        cell.write_images()  # Uncomment this line if segmentations are performed
        cell.plot_images(ax)
        pbar.update(1)


plt.tight_layout()
fig.savefig(data_root / "segmentation_results.pdf")
