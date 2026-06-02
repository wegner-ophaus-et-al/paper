from pathlib import Path
import pickle
import matplotlib.pyplot as plt

data_dir = Path(__file__).parent / "data"

cell_objects = []
for cell_dir in data_dir.iterdir():
    if not cell_dir.is_dir():
        continue
    with open(cell_dir / "thecell.pkl", "rb") as f:
        cell = pickle.load(f)
        cell_objects.append(cell)


def contact_sheet(list_of_cells: list, save_dir: Path):
    """
    Generate a contanct sheet for all cells in the list
    """
    fig, axs = plt.subplots(len(list_of_cells), 3, figsize=(10, 4 * len(list_of_cells)))
    for ax, cell in zip(axs, list_of_cells):
        cell.plot_markers_on_axis(ax, segmentation_cmap="summer")
        for a in ax:
            a.axis("off")
    fig.savefig(save_dir / "contact_sheet.pdf")


for cell in cell_objects:
    for marker in cell.markers.values():
        if marker.compartment.lower() == "granules":
            marker.features = marker.get_features()
            print(marker.features.keys())
