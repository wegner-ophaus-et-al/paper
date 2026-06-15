from thecell import TheCell
import sys

sys.path.append("/Users/icb_remote/Documents/JW/py/packages/")

from gamgee.instance import ModelHandler
from pathlib import Path
from multiprocessing import Pool
import pickle
from tqdm import tqdm

model_handler = ModelHandler()
model_handler.add_model(
    model_name="granules",
    cell_compartment="granules",
    model_type="vit_l_lm",
    model_path="/Users/icb_remote/Documents/JW/py/packages/gamgee/models/msam/sam_granules_refined_up3_39637264/",
)


sample_dir = Path(
    "/Users/icb_remote/Documents/JW/py/data/iPGC_24hpf_tdrd7a-modulation/fm_nt_data"
)
cell_objects = []

# Initialize cell objects
for test_cell in sample_dir.iterdir():
    if test_cell.is_dir() and "figures" not in test_cell.name:
        cell_obj = TheCell(
            test_cell,
            conditions=["full_mix", "no_tdrd7"],
            model_handler=model_handler,
            granuleA="dnd",
            granuleB="gra",
        )
        cell_objects.append(cell_obj)

# Segment marks in paralel
for cell in tqdm(cell_objects, desc="Segmenting Cells", total=len(cell_objects)):
    cell.read_segmentations()
    for marker in cell.markers.values():
        if "cell" in marker.name or "nucleus" in marker.name or "nls" in marker.name:
            continue
        # print(f"Segmenting {marker.name} in cell {cell.path.name}")
        # future = pool.apply_async(marker.segment)
        marker.segment()
    cell.write_segmentations()
    with open(cell.path / "thecell.pkl", "wb") as f:
        pickle.dump(cell, f)
