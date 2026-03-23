from scipy.ndimage import distance_transform_edt
from pathlib import Path
import tifffile as tiff
from scipy import stats

def compute_masks(seg_cell, seg_nucleus, seg_granule, safety_margin=2):
    """
    Compute binary masks for cell, nucleus, granules, and cytoplasm based on segmentations. The cytoplasm mask is defined as the cell mask minus the nucleus and granule masks, with an optional safety margin to avoid overlap.
    Args:
        seg_cell:
        seg_nucleus:
        seg_granule:
        safety_margin:

    Returns: Cell mask, nucleus mask, granule mask, cytoplasm mask

    """
    # Create binary masks
    cell_mask = seg_cell > 0
    nucleus_mask = seg_nucleus > 0
    granule_mask = seg_granule > 0

    # Only consider granules within cells
    granule_mask = granule_mask & cell_mask

    # Only consider nucleus within cells
    nucleus_mask = nucleus_mask & cell_mask

    # Cytoplasm mask is cell minus nucleus and granules with safety margin
    nucleus_granule_mask = nucleus_mask | granule_mask
    distance_map = distance_transform_edt(~nucleus_granule_mask)
    cytoplasm_mask = cell_mask & (distance_map > safety_margin)
    return cell_mask, nucleus_mask, granule_mask, cytoplasm_mask


def get_confocal_pixel_size(sample_path:Path):
    """
    Get the pixel size from a confocal microscopy image. The pixel size is typically stored in the metadata of the image file. This function reads the metadata and extracts the pixel size information.

    Args:
        sample_path (Path): Path a sample that contains a confocal microscopy image file (.lsm).
    Returns:
        tuple: A float that represents the pixel sizein micrometers (µm) for the x and y dimensions.
    """
    # Search the sample directory and its subdirectories for .lsm files
    lsm_dir = list(sample_path.rglob('*.lsm'))
    if not lsm_dir:
        raise FileNotFoundError(f"No .lsm files found in {sample_path} or its subdirectories.")

    with tiff.TiffFile(lsm_dir[0]) as lsm_file:
        # Get the metadate
        metadata = lsm_file.lsm_metadata

        # Check if the x and y pixel sizes are available in the metadata and are equal
        x_px_size = metadata.get('VoxelSizeX', None)
        y_px_size = metadata.get('VoxelSizeY', None)

        if x_px_size == y_px_size and x_px_size is not None:
            return x_px_size
        else:
            raise ValueError("Pixel size information is missing or inconsistent in the metadata.")


def parametric(data_set1, data_set2):
    set1_normal = stats.normaltest(data_set1).pvalue > 0.05
    set2_normal = stats.normaltest(data_set2).pvalue > 0.05

    return min(set1_normal, set2_normal)


def statistical_analysis(data_set1, data_set2):
    if parametric(data_set1, data_set2):
        t_statistic, p_value = stats.ttest_ind(data_set1, data_set2)
        test_type = "t-test"
    else:
        t_statistic, p_value = stats.mannwhitneyu(data_set1, data_set2)
        test_type = "Mann-Whitney U test"
    return test_type, t_statistic, p_value

def get_stat_stars(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'