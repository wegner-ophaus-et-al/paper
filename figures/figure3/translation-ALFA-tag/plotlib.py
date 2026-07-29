import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle


COLORS = {
    "cyan": (0, 1, 1),
    "magenta": (1, 0, 1),
    "yellow": (1, 1, 0),
    "grey": (1, 1, 1),
}


def _cmap(name):
    return LinearSegmentedColormap.from_list(name, [(0, 0, 0), COLORS[name]])


def plot_image(image, ax, colormap: str = "grey", vmin=None, vmax=None):
    return ax.imshow(image, cmap=_cmap(colormap), vmin=vmin, vmax=vmax)


def plot_merge(images: dict, ax, vmins: dict = None, vmaxs: dict = None):
    vmins, vmaxs = vmins or {}, vmaxs or {}
    rgb = 0
    for name, img in images.items():
        norm = Normalize(vmins.get(name), vmaxs.get(name), clip=True)(img)
        rgb = rgb + np.asarray(norm)[..., None] * np.array(COLORS[name])
    return ax.imshow(np.clip(rgb, 0, 1))


def scale_bar(ax, length, thickness, pixel_size, pad=0.03, color="white"):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    sx, sy = np.sign(x1 - x0), np.sign(y1 - y0)
    px = length / pixel_size
    x = x1 - sx * (pad * abs(x1 - x0) + px)
    y = y0 + sy * pad * abs(y1 - y0)
    ax.add_patch(Rectangle((x, y), sx * px, sy * thickness, color=color, zorder=10))
