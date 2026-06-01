from pathlib import Path
import pickle
import matplotlib.pyplot as plt

data_dir = Path(__file__).parent / "data"

cell_objects = []
for cell_dir in data_dir.iterdir():
    with open(cell_dir / "thecell.pkl", "rb") as f:
        cell = pickle.load(f)
        cell_objects.append(cell)


# Generate a contact sheet for all cells
fig, ax = plt.subplots(len(cell_objects), 1, figsize=(10, 5 * len(cell_objects)))
for ax, cell in zip(ax, cell_objects):
    cell.plot_markers_on_axis(ax)
