from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import print_nested, export_representative_cells, pub_images
from multiprocessing import Pool
from plotting import (
    plot_foldchange,
    plot_individial_granule_profile,
    contact_sheet,
    per_cell_summary,
    ridgeplot_per_marker,
    cell_features,
)


data_root = Path(
    "/Users/icb_remote/Documents/JW/py/data/iPGC_24hpf_tdrd7a-modulation/fm_nt_data"
)

only_data_processing = False
redo_feature_extraction = True


def process_cell(args):
    cell_counter, cell, total = args
    if not cell.cell_segmentation_exists():
        return [[], {}]
    if redo_feature_extraction:
        print(f"Cell {cell_counter + 1}/{total}")
        cell.read_segmentations()
        cell.clean_segmentations()
        cell.compute_features()
        # cell.write_segmentations()

        with open(cell.path / "thecell.pkl", "wb") as f:
            pickle.dump(cell, f)

    return (cell.get_granule_features(), cell.get_cell_features())


def main():
    if not only_data_processing:
        cell_objects = []
        for data_dir in [data_root]:
            for cell_dir in data_dir.iterdir():
                if not cell_dir.is_dir() or "figures" in cell_dir.name:
                    continue
                with open(cell_dir / "thecell.pkl", "rb") as f:
                    cell = pickle.load(f)
                    cell_objects.append(cell)

        data = []
        data_cell = []
        with Pool() as pool:
            for features in pool.map(
                process_cell,
                [
                    (cell_counter, cell, len(cell_objects))
                    for cell_counter, cell in enumerate(cell_objects)
                ],
            ):
                data.extend(features[0])
                data_cell.append(features[1])

        df = pd.DataFrame(data)
        df.to_pickle(data_root / "granule_features.pkl")
        # df.to_csv(data_root / "granule_features.csv")

        # print(data_cell)
        df_cell = pd.DataFrame(data_cell)
        df_cell.to_pickle(data_root / "cell_features.pkl")

        if redo_feature_extraction:
            contact_sheet(cell_objects, data_root)

        print("---" * 20)
        print_nested(cell_objects[0].features, indent=2)
        print("---" * 20)
    else:
        df = pd.read_pickle(data_root / "granule_features.pkl")
        df_cell = pd.read_pickle(data_root / "cell_features.pkl")

    # Calc aspect ratio
    df["AspectRatio"] = df["MajorAxisLength"] / df["MinorAxisLength"]

    # Acquisition parameters
    magnification = 63
    pixel_pitch = 6.5
    binning = 1
    pixel_size = pixel_pitch * binning / magnification  # in microns

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
        "full_mix": "#878787",
        "no_tdrd7a": "#8a38a6",
    }

    # plot_individial_granule_profile(df, data_root)

    per_cell_summary(df, data_root, color_palette=color_palette)
    plot_foldchange(df, data_root)
    cell_features(df_cell, data_root, color_palette=color_palette)

    export_representative_cells(
        df, "SphericalVolume", data_root / "figures" / "representative_cells"
    )

    if not only_data_processing and cell_objects:
        # Representative cells for SphericalVolume feature, exported by export_representative_cells function
        representative_cells = [
            uid.lower()
            for uid in [
                "02c4431d",
                "65a4e350",
                "6e528d49",
                "758d46a5",
                "7b386751",
                "8081a738",
                "b7c19c49",
                "168d3394",
                "3b2cbed0",
                "6a8e011a",
                "716495ad",
                "91f2ff2a",
                "da3efd6b",
            ]
        ]

        ncols = 3
        nrows = len(representative_cells)
        col_width = 2
        row_height = 2
        fig_representative, ax_representative = plt.subplots(
            nrows=nrows, ncols=ncols, figsize=(col_width * ncols, row_height * nrows)
        )

        for rep_uid, ax in zip(representative_cells, ax_representative):
            for cell in cell_objects:
                if cell.uid.lower() == rep_uid.lower():
                    pub_images(ax, cell)

        fig_representative.tight_layout()
        path_representative = (
            data_root / "figures" / "representative_cells" / "representative_cells.pdf"
        )
        path_representative.parent.mkdir(parents=True, exist_ok=True)
        fig_representative.savefig(path_representative, dpi=600)


if __name__ == "__main__":
    main()
