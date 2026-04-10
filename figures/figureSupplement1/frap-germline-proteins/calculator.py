import numpy as np


def phair_double_normalization(
    image_sequence: np.ndarray,
    irradiated_mask: np.ndarray,
    reference_mask: np.ndarray,
    background_mask: np.ndarray,
    pre_bleach_frames: int = 5,
):
    """
    Perform PHair double normalization on a sequence of images using specified masks.

    Parameters:
    image_sequence (np.ndarray): 3D array of shape (time, x, y) representing the image sequence.
    irradiated_mask (np.ndarray): 2D boolean array indicating the irradiated region.
    correction_mask (np.ndarray): 2D boolean array indicating the correction region.
    background_mask (np.ndarray): 2D boolean array indicating the background region.

    Returns:
    np.ndarray: Normalized image sequence after applying PHair double normalization.
    """

    intensity_irradiated = image_sequence[:, irradiated_mask].mean(axis=1)
    intensity_reference = image_sequence[:, reference_mask].mean(axis=1)
    intensity_background = image_sequence[:, background_mask].mean(axis=1)

    irr_norm = intensity_irradiated - intensity_background
    ref_norm = intensity_reference - intensity_background

    irr_norm_pre_bleach = irr_norm[:pre_bleach_frames].mean()
    ref_norm_pre_bleach = ref_norm[:pre_bleach_frames].mean()

    double_normalized = (ref_norm_pre_bleach / irr_norm_pre_bleach) * (
        irr_norm / ref_norm
    )

    return double_normalized

