import sys
sys.path.append('/Users/julian/Documents/General Science/Programming/py/packages')
from gamgee.segmenter import SegmentationModel, SegmentationInstance
from pathlib import Path
import matplotlib.pyplot as plt
import tifffile as tiff


sm_g = SegmentationModel('/Users/julian/Library/Mobile Documents/com~apple~CloudDocs/Documents/General Science/Programming/py/general_analysis/20250722_msam_granule_sizes/gamgee/models/msam/granules/sam_granules_refined_up3_35416497',
                       model_type='vit_l_lm')


root = Path("/Users/julian/local_files/chx_kim/data")

list_of_sample_folders = [path for path in root.iterdir() if path.is_dir() and not path.name.startswith(".")]


fig, ax = plt.subplots(len(list_of_sample_folders), 3, figsize=(4, 2*len(list_of_sample_folders)))

for i, sample_folder_path in enumerate(list_of_sample_folders):
    print(f"Processing {sample_folder_path.name}...")
    image_file_path = next((sample_folder_path / "original_raw").glob("*.lsm"))

    if not image_file_path:
        print(f"No .lsm file found in {sample_folder_path}. Skipping.")
        continue

    image = tiff.imread(image_file_path)

    img_out_path = sample_folder_path / "imgs"
    seg_out_path = sample_folder_path / "segmentations"

    for out_path in [img_out_path, seg_out_path]:
        out_path.mkdir(exist_ok=True)

    # Vasa channel segmentation and writing
    vasa_seg = SegmentationInstance(image[0], segmentation_model=sm_g)
    tiff.imwrite(img_out_path / "vasa.tif", image[0])
    tiff.imwrite(seg_out_path / "vasa.tif", vasa_seg.segmentation)

    # Nanos channel segmentation and writing
    nanos_seg = SegmentationInstance(image[1], segmentation_model=sm_g)
    tiff.imwrite(img_out_path / "nanos.tif", image[1])
    tiff.imwrite(seg_out_path / "nanos.tif", nanos_seg.segmentation)

    # tdrd7a channel segmentation and writing
    tdrd7a_seg = SegmentationInstance(image[2], segmentation_model=sm_g)
    tiff.imwrite(img_out_path / "tdrd7a.tif", image[2])
    tiff.imwrite(seg_out_path / "tdrd7a.tif", tdrd7a_seg.segmentation)

    if image.ndim == 3:
        ax[i, 0].imshow(image[0], cmap='gray')
        ax[i, 0].contour(vasa_seg.segmentation, colors='r', linewidths=0.5)
        ax[i, 1].imshow(image[1], cmap='gray')
        ax[i, 1].contour(nanos_seg.segmentation, colors='r', linewidths=0.5)
        ax[i, 2].imshow(image[2], cmap='gray')
        ax[i, 2].contour(tdrd7a_seg.segmentation, colors='r', linewidths=0.5)
        ax[i, 0].set_title("Channel 0 - vasa")
        ax[i, 1].set_title("Channel 1 - nanos")
        ax[i, 2].set_title("Channel 2 - tdrd7a")
    else:
        print(f"Unexpected image dimensions {image.shape} in {image_file_path}. Skipping.")


    for a in ax.flatten():
        a.axis('off')
    plt.tight_layout()
    fig.savefig('contactsheed.pdf', dpi=400)


