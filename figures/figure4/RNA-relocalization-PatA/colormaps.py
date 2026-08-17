import matplotlib as mpl

# Create custom colormaps from black to each color
cyan = mpl.colors.LinearSegmentedColormap.from_list("cyan", [(0, 0, 0), (0, 1, 1)], N=256)
magenta = mpl.colors.LinearSegmentedColormap.from_list("magenta", [(0, 0, 0), (1, 0, 1)], N=256)
yellow = mpl.colors.LinearSegmentedColormap.from_list("yellow", [(0, 0, 0), (1, 1, 0)], N=256)
green = mpl.colors.LinearSegmentedColormap.from_list("green", [(0, 0, 0), (0, 1, 0)], N=256)
red = mpl.colors.LinearSegmentedColormap.from_list("red", [(0, 0, 0), (1, 0, 0)], N=256)