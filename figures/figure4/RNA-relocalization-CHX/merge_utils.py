import numpy as np
import matplotlib.colors as mcolors



def merge_images_with_cmaps(img1, img2, cmap1, cmap2, vmin1=None, vmax1=None, vmin2=None, vmax2=None):
    """
    Merge two images while preserving their individual colormaps.

    :param img1: First image array
    :param img2: Second image array
    :param cmap1: Colormap for first image (e.g., magenta)
    :param cmap2: Colormap for second image (e.g., yellow)
    :param vmin1, vmax1: Min/max values for first image normalization
    :param vmin2, vmax2: Min/max values for second image normalization
    :return: Merged RGB image
    """
    # Normalize images
    norm1 = mcolors.Normalize(vmin=vmin1 or img1.min(), vmax=vmax1 or img1.max())
    norm2 = mcolors.Normalize(vmin=vmin2 or img2.min(), vmax=vmax2 or img2.max())

    # Convert to RGB using colormaps
    rgb1 = cmap1(norm1(img1))[..., :3]  # Remove alpha channel
    rgb2 = cmap2(norm2(img2))[..., :3]

    # Merge by adding RGB values (clipped to [0, 1])
    merged = np.clip(rgb1 + rgb2, 0, 1)

    return merged
