import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_laplace, label
from skimage.morphology import remove_small_objects


def dog_filter(img, sigma, ratio=1.6):
    """Difference of Gaussians. sigma ~ PSF sigma in pixels (blob radius / sqrt(2))."""
    img = img.astype(np.float32)
    return gaussian_filter(img, sigma) - gaussian_filter(img, sigma * ratio)


def log_filter(img, sigma):
    """Laplacian of Gaussian, sign-flipped so bright blobs give positive response."""
    return -gaussian_laplace(img.astype(np.float32), sigma)


def anscombe(img):
    """Variance-stabilizing transform for shot-noise-dominated data.
    Apply to the raw offset-subtracted image, before filtering."""
    return 2.0 * np.sqrt(np.maximum(img, 0) + 3.0 / 8.0)


def mad_threshold(img, k=5.0):
    """Robust background + k * robust sigma. Returns a scalar threshold."""
    bg = np.median(img)
    sigma = 1.4826 * np.median(np.abs(img - bg))
    return bg + k * sigma


def segment_spots(img, sigma=1.5, k=5.0, min_area=4, offset=0.0, stabilize=False):
    """Returns (labels, n_objects, threshold_used).

    offset    : camera baseline in ADU, from a dark frame
    k         : calibrate on your no-signal condition, then freeze it
    min_area  : drop sub-PSF-sized survivors
    """
    img = img.astype(np.float32) - offset
    if stabilize:
        img = anscombe(img)

    filtered = dog_filter(img, sigma)
    thr = mad_threshold(filtered, k)

    mask = filtered > thr
    if min_area > 1:
        mask = remove_small_objects(mask, min_size=min_area)

    labels, n = label(mask)
    return labels, n, thr


def sweep_k(images, ks=np.arange(2, 10, 0.5), **kwargs):
    """False-positive curve on negative controls: objects detected vs k."""
    return {
        float(k): [segment_spots(im, k=k, **kwargs)[1] for im in images] for k in ks
    }
