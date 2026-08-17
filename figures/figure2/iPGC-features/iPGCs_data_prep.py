import marimo

__generated_with = "0.23.13"
app = marimo.App()


@app.cell
def _():
    import sqlite3 as sql
    import pandas as pd
    from pathlib import Path
    import matplotlib.pyplot as plt
    import tifffile as tiff
    import numpy as np
    import pickle
    import sys

    return Path, plt, sys


@app.cell
def _(Path):
    root = Path('/Volumes/HELHEIM/analyzed_data/size_and_periphery/iPGC_24hpf_tdrd7a-modulation')
    return (root,)


@app.cell
def _():
    # db = sql.connect(root / 'iPGC_24hpf_tdrd7a-modulation.db')
    #
    # # Create columns for uid (UNIQUE), condition, and filename and absolute path
    # query = """
    # CREATE TABLE IF NOT EXISTS samples (
    #     uid TEXT PRIMARY KEY,
    #     condition TEXT,
    #     filename TEXT,
    #     absolute_path TEXT
    # );
    # """
    #
    # with db:
    #     db.execute(query)
    return


@app.cell
def _():
    # condition_id_dict = {
    #     "full_mix": ["_full-mix", "_full", "_fm"],
    #     "no_tdrd7a": ["_no-tdrd7a", "_no-tdrd", "_nt"],
    #     "tmd": ["_tmd"],
    # }
    return


@app.cell
def _():
    # for sample in root.iterdir():
    #     if not sample.is_dir() or sample.name.startswith('.'):
    #         continue
    #     print(f'Processing sample folder: {sample.name}')
    #     if not sample.name.count('__') == 1:
    #         print(f"Sample {sample.name} is already preped. Skipping...")
    #         continue
    #     uid, filename = sample.name.split('__')
    #
    #     condition = None
    #     # Set condition based on the identifiers in the folder name
    #     for cond, identifiers in condition_id_dict.items():
    #         if any(identifier in sample.name.lower() for identifier in identifiers):
    #             condition = cond
    #             break
    #     if condition is None:
    #         print(f"Could not determine condition for sample {sample.name}. Skipping...")
    #         continue
    #
    #     # Rename folder and file to {uid}__{condition}__{filename}
    #     file_path = sample / "original_raw" / f"{filename}.tif"
    #     if file_path.exists():
    #         new_folder_name = f"{uid}__{condition}__{filename}"
    #         new_file_name = f"{uid}__{condition}__{filename}.tif"
    #         new_folder_path = root / new_folder_name
    #         new_file_path = new_folder_path / "original_raw"/ new_file_name
    #
    #
    #
    #         # Move and rename the file
    #         sample.rename(new_folder_path)  # Rename the folder first
    #         (new_folder_path / "original_raw"/ f"{filename}.tif").rename(new_file_path)
    #
    #         # Add entry to the database
    #         with db:
    #             db.execute(
    #                 "INSERT OR IGNORE INTO samples (uid, condition, filename, absolute_path) VALUES (?, ?, ?, ?)",
    #                 (uid, condition, filename, str(new_folder_path)),
    #             )
    #
    #     else:
    #         print(f"File {file_path} does not exist. Skipping...")
    return


@app.cell
def _():
    # # Double check if all folder have an entry if not add
    # for sample in root.iterdir():
    #     if not sample.is_dir() or sample.name.startswith('.'):
    #         continue
    #     # print(f'Checking sample folder: {sample.name}')
    #     if not sample.name.count('__') == 2:
    #         # print(f"Sample {sample.name} does not follow the naming convention. Skipping...")
    #         continue
    #     uid, condition, filename = sample.name.split('__')
    #     with db:
    #         result = db.execute("SELECT * FROM samples WHERE uid = ?", (uid,)).fetchone()
    #         if result is None:
    #             print(f"No entry found for {sample.name}. Adding to database...")
    #             db.execute(
    #                 "INSERT OR IGNORE INTO samples (uid, condition, filename, absolute_path) VALUES (?, ?, ?, ?)",
    #                 (uid, condition, filename, str(sample)),
    #             )
    return


@app.cell
def _():
    # # Save contact sheet as pdf
    # fig, axes = None, None
    # channel_counter = 0
    #
    # number_of_samples = sum(1 for sample in root.iterdir() if sample.is_dir() and not sample.name.startswith('.'))
    # print(f"Number of samples: {number_of_samples}")
    # for sample_dir in root.iterdir():
    #     if not sample_dir.is_dir() or sample_dir.name.startswith('.'):
    #         continue
    #     uid, condition, filename = sample_dir.name.split('__')
    #     img_path  = sample_dir / "original_raw" / f"{sample_dir.name}.tif"
    #     if img_path.exists():
    #         img = tiff.imread(img_path)
    #         if img.ndim == 3:
    #             n_channels = img.shape[0]
    #             export = True
    #         elif img.ndim == 4:
    #             # Project z-stacks
    #             img = np.max(img, axis=0)  # Sum across z-stacks if 4D
    #             n_channels = img.shape[0]
    #             export = True
    #         else:
    #             n_channels = 1
    #             export = False
    #
    #         if export:
    #             img_dir = sample_dir / "imgs"
    #             img_dir.mkdir(exist_ok=True, parents=True)
    #             for ch in range(n_channels):
    #                 img[ch] = img[ch] / np.max(img[ch]) * 2**16 - 1  # Normalize to 16-bit range
    #             if "2025-03-13" in sample_dir.name:
    #                 tiff.imwrite(img_dir / f"dnd.tif", img[0], dtype=np.uint16)
    #                 tiff.imwrite(img_dir / f"nls.tif", img[1], dtype=np.uint16)
    #                 tiff.imwrite(img_dir / f"gra.tif", img[2], dtype=np.uint16)
    #             else:
    #                 tiff.imwrite(img_dir / f"nls.tif", img[0], dtype=np.uint16)
    #                 tiff.imwrite(img_dir / f"dnd.tif", img[1], dtype=np.uint16)
    #                 tiff.imwrite(img_dir / f"gra.tif", img[2], dtype=np.uint16)
    #
    #         if fig is None:
    #             fig, axes = plt.subplots(number_of_samples, n_channels, figsize=(5*n_channels, 5*number_of_samples))
    #             channel_counter = 0
    #         else:
    #             channel_counter += 1
    #         if n_channels == 1:
    #             axes.imshow(img, cmap='gray')
    #             axes.set_ylabel(f"{uid} - {condition}")
    #             axes.axis('off')
    #         else:
    #             for ch in range(n_channels):
    #                 axes[channel_counter,ch].imshow(img[ch], cmap='gray')
    #                 # axes[channel_counter, ch].set_title(f"Channel {ch+1}")
    #                 axes[channel_counter, ch].axis('off')
    #             axes[channel_counter, 0].set_title(f"{uid}")
    #             axes[channel_counter, 1].set_title(f"{condition}")
    #
    #     else:
    #         print(f"File {img_path} does not exist. Skipping...")
    # plt.tight_layout()
    # plt.savefig(root / "contact_sheet.pdf")
    # plt.close(fig)
    # # %%
    return


@app.cell
def _(sys):
    sys.path.append('..')
    from gamgee.instance import ModelHandler, TheCell

    return ModelHandler, TheCell


@app.cell
def _(ModelHandler):
    mh = ModelHandler()
    # mh.preconfigureations["granules"] = dict(path="/Users/julian/local_files/microSAM/sam_granules_refined_up3_35558740",
    #                                          model_type="vit_l_lm", friendly_name="35416497",
    #                                          cell_compartment="granules", upsample_factor=3)
    # mh.preconfigureations["cell_nls"] = dict(path="/Users/julian/local_files/microSAM/sam_large_blobs_up1_35569278",
    #                                          model_type="vit_l_lm", friendly_name="35569278",
    #                                          cell_compartment="cell", upsample_factor=1)
    return (mh,)


@app.cell
def _(TheCell, mh, root):
    for sample_dir in root.iterdir():
        if not sample_dir.is_dir() or sample_dir.name.startswith('.'):
            continue
        uid, condition, filename = sample_dir.name.split('__')
        print(f"Processing {uid} - {condition}")
        tc = TheCell(
            uid=uid,
            name=filename,
            root_path=sample_dir,
            model_handler=mh,
        )
        # tc.save_instance()
        tc.save_segmentations()
        tc.pickle_images()
        break  # Process only the first sample for testing
    return (tc,)


@app.cell
def _(mh, sys):
    # Calculate size of tc object
    size_in_bytes = sys.getsizeof(mh.segmentation_models)
    print(f"Size of TheCell object: {size_in_bytes} bytes")
    print(f"That's {size_in_bytes/1000000} MB for a family of four!")
    return


@app.cell
def _(plt, tc):
    plt.imshow(tc.markers["nls"].segmentation)
    return


@app.cell
def _(tc):
    tc.markers["nls"].logs["Segmentation Info"]["checkpoints_name"]
    return


@app.cell
def _(tc):
    tc.logs
    return


if __name__ == "__main__":
    app.run()
