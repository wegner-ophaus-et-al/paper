import sys
sys.path.append('/Users/julian/Documents/General Science/Programming/py/packages')
from gamgee.segmenter import SegmentationModel
from pathlib import Path
import matplotlib.pyplot as plt
import tifffile as tiff


sm_g = SegmentationModel('/Users/julian/Library/Mobile Documents/com~apple~CloudDocs/Documents/General Science/Programming/py/general_analysis/20250722_msam_granule_sizes/gamgee/models/msam/granules/sam_granules_refined_up3_35416497',
                       model_type='vit_l_lm')


root = Path("/Users/julian/local_files/chx_kim/data")

list_of_sample_folders = [path for path in root.iterdir() if path.is_dir() and not path.name.startswith(".")]


fig, ax = plt.subplots(len(list_of_sample_folders), 3, figsize=(2, 2*len(list_of_sample_folders)))

for idx, sample_folder_path in enumerate(list_of_sample_folders):
    print(f"Processing {sample_folder_path.name}...")
    image_file_path = next((sample_folder_path / "original_raw").glob("*.lsm"))

    if not image_file_path:
        print(f"No .lsm file found in {sample_folder_path}. Skipping.")
        continue

    image = tiff.imread(image_file_path)

