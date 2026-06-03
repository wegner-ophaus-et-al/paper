from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import print_nested
from plotting import plot_individial_granule_profile, contact_sheet


data_dir = Path(__file__).parent / "data"

cell_objects = []
for cell_dir in data_dir.iterdir():
    if not cell_dir.is_dir():
        continue
    with open(cell_dir / "thecell.pkl", "rb") as f:
        cell = pickle.load(f)
        cell_objects.append(cell)


data = []
for cell in cell_objects:
    # cell.clean_segmentations()
    # cell.compute_features()

    with open(cell.path / "thecell.pkl", "wb") as f:
        pickle.dump(cell, f)

    data.extend(cell.get_granule_features())


df = pd.DataFrame(data)
df.to_pickle(data_dir / "granule_features.pkl")
df.to_csv(data_dir / "granule_features.csv")

plot_individial_granule_profile(df, data_dir)
