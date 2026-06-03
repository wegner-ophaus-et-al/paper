from thecell import TheCell
import sys

sys.path.append("/Users/icb_remote/Documents/JW/py/packages/")

from gamgee.instance import ModelHandler
from pathlib import Path
from multiprocessing import Pool
import pickle
from tqdm import tqdm

model_handler = ModelHandler()


sample_dir = Path(
    "/Volumes/icb_remote/Documents/JW/py/data/endogenous_PGCs_size_characterization/8hpf"
)
cell_objects = []

# Initialize cell objects
for test_cell in sample_dir.iterdir():
    if test_cell.is_dir():
        print(f"Processing cell at {test_cell}")
        cell_obj = TheCell(test_cell, conditions=["ctrl", "kd"])
        cell_objects.append(cell_obj)

# Segment marks in paralel
with Pool(processes=3) as pool:
    for cell in tqdm(cell_objects, desc="Segmenting Cells", total=len(cell_objects)):
        cell.read_segmentations()
        for marker in cell.markers.values():
            # Remove this line when the model works more clear
            if marker.name == "cell" or marker.name == "nucleus":
                continue
            future = pool.apply_async(marker.segment)
        cell.write_segmentations()

# Save cells as objects
for cell in cell_objects:
    with open(cell.path / "thecell.pkl", "wb") as f:
        pickle.dump(cell, f)
