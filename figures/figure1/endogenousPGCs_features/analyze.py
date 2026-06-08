from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import print_nested
from plotting import (
    plot_individial_granule_profile,
    contact_sheet,
    per_cell_summary,
    ridgeplot_per_marker,
)


data_root = Path(
    "/Users/icb_remote/Documents/JW/py/data/endogenous_PGCs_size_characterization"
)

redo_feature_extraction = False

if redo_feature_extraction:
    cell_objects = []
    for data_dir in [data_root / "8hpf", data_root / "24hpf"]:
        for cell_dir in data_dir.iterdir():
            if not cell_dir.is_dir() or "figures" in cell_dir.name:
                continue
            with open(cell_dir / "thecell.pkl", "rb") as f:
                cell = pickle.load(f)
                cell_objects.append(cell)

    data = []
    for cell_counter, cell in enumerate(cell_objects):
        if not cell.cell_segmentation_exists():
            continue
        print(f"Cell {cell_counter + 1}/{len(cell_objects)}")
        cell.read_segmentations()
        cell.clean_segmentations()
        cell.compute_features()

        with open(cell.path / "thecell.pkl", "wb") as f:
            pickle.dump(cell, f)

        data.extend(cell.get_granule_features())

    df = pd.DataFrame(data)
    df.to_pickle(data_root / "granule_features.pkl")
    df.to_csv(data_root / "granule_features.csv")
    contact_sheet(cell_objects, data_root)
else:
    df = pd.read_pickle(data_root / "granule_features.pkl")

# Calc aspect ratio
df["AspectRatio"] = df["MajorAxisLength"] / df["MinorAxisLength"]

# Acquisition parameters
magnification = 63
pixel_pitch = 6.5
binning = 1
pixel_size = pixel_pitch * binning / magnification  # in microns

for col in df.columns:
    if "fourier" in col.lower():
        continue

# Scale features to um
for lin_feature in [
    "Perimeter",
    "MajorAxisLength",
    "MinorAxisLength",
    "EdgeDistanceToNucleus",
    "EdgeDistanceToNucleusSigned",
    "CentroidDistanceToNucleus",
    "CentroidDistanceToNucleusSigned",
]:
    df[lin_feature] = df[lin_feature] * pixel_size
df["Area"] = df["Area"] * pixel_size**2
df["TouchAreaNucleus"] = df["TouchAreaNucleus"] * pixel_size**2
for vol_feature in [
    "SphericalVolume",
    "EllipsoidVolumeProlate",
    "EllipsoidVolumeOblate",
]:
    df[vol_feature] = df[vol_feature] * pixel_size**3

color_palette = {
    "ctrl": "#878787",
    "kd": "#8a38a6",
}


plot_individial_granule_profile(df, data_root)

per_cell_summary(df, data_root, color_palette=color_palette)
ridgeplot_per_marker(df, data_root)
