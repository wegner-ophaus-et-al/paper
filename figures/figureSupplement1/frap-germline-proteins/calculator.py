import numpy as np


def phair_double_normalization(
    image_sequence: np.ndarray,
    irradiated_mask: np.ndarray,
    reference_mask: np.ndarray,
    background_mask: np.ndarray,
    pre_bleach_frames: int = 5,
    clip_pre_bleach: bool = True,
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

    refernce_mask_xor = np.logical_xor(irradiated_mask > 0, reference_mask > 0)

    intensity_irradiated = image_sequence[:, irradiated_mask > 0].mean(axis=1)
    intensity_reference = image_sequence[:, refernce_mask_xor > 0].mean(axis=1)
    intensity_background = image_sequence[:, background_mask > 0].mean(axis=1)

    irr_norm = intensity_irradiated - intensity_background
    ref_norm = intensity_reference - intensity_background

    irr_norm_pre_bleach = irr_norm[: pre_bleach_frames - 1].mean()
    ref_norm_pre_bleach = ref_norm[: pre_bleach_frames - 1].mean()

    double_normalized = (ref_norm_pre_bleach / irr_norm_pre_bleach) * (
        irr_norm / ref_norm
    )

    if clip_pre_bleach:
        # Normalize pre-bleach values to 1 and first post-bleach value to 0
        tripple_normalized = (
            double_normalized - double_normalized[pre_bleach_frames]
        ) / (
            double_normalized[: pre_bleach_frames - 1].mean()
            - double_normalized[pre_bleach_frames]
        )
    else:
        tripple_normalized = double_normalized
    return {
        "intensity_irradiated_raw": intensity_irradiated,
        "intensity_reference_raw": intensity_reference,
        "intensity_background_raw": intensity_background,
        "intensity_normalized": tripple_normalized,
    }


def subsample(arr, n):
    indices = np.linspace(0, len(arr) - 1, n, dtype=int)
    return arr[indices]
