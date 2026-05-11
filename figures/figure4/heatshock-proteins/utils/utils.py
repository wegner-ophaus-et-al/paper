import numpy as np
from scipy import stats
from pathlib import Path
from PIL import Image
import tifffile as tiff


def upsampling(image, scale_factor):
    """
    Fast vectorized upsampling of a 2D image using NumPy.

    Parameters:
    -----------
    image : numpy.ndarray
        Input 2D image array
    scale_factor : int
        Factor by which to scale the image

    Returns:
    --------
    numpy.ndarray
        Upsampled image
    """
    if image.ndim != 2:
        raise ValueError("Input image must be 2D")

    # Repeat each row scale_factor times
    upsampled = np.repeat(image, scale_factor, axis=0)
    # Repeat each column scale_factor times
    upsampled = np.repeat(upsampled, scale_factor, axis=1)

    return upsampled


def preprocess(image_data, upsample_factor=3, final_size=450, verbose=False):
    """
    Process the image data by upsampling and saving it.

    :param final_size:
    :param image_data: The image data to process.
    :param upsample_factor: The factor by which to upscale the image.
    :return: None
    """
    if not isinstance(image_data, np.ndarray):
        raise ValueError("Image data must be a numpy array.")

    if image_data.ndim != 2:
        raise ValueError("Image data must be 2D.")

    if image_data.shape[0] == image_data.shape[1]:
        if image_data.shape[0] < final_size:
            # If the image is smaller than the final size, pad it
            padding_total = final_size - image_data.shape[0]
            padding_x = padding_total // 2
            padding_y = padding_x if padding_x * 2 == padding_total else padding_x + 1

            image_padded = np.pad(image_data, (padding_x, padding_y), mode="edge")

            if image_padded.shape[0] < final_size:
                # If padding still results in a size smaller than final_size, raise an error
                raise ValueError(
                    f"Image data is too small after padding: {image_padded.shape}."
                )
        elif image_data.shape[0] > final_size:
            # If the image is larger than or equal to the final size, crop it
            start = (image_data.shape[0] - final_size) // 2
            end = start + final_size
            image_padded = image_data[start:end, start:end]

            if verbose:
                print(
                    f"Cropping image to final size: {final_size}x{final_size} from original size: {image_data.shape}"
                )
        elif image_data.shape[0] == final_size:
            # If the image is already the correct size, no padding or cropping needed
            image_padded = image_data
            if verbose:
                print(
                    f"Image is already at final size: {final_size}x{final_size}. No processing needed."
                )
        else:
            raise ValueError(f"Image≠≠")
    else:
        raise ValueError(f"Image data must be square but has shape {image_data.shape}.")

    # Upsample the image
    upsampled_image = upsampling(image_padded, upsample_factor)

    # Set range to [0, 255] for uint8 images
    if upsampled_image.dtype != np.uint8:
        upsampled_image = (
            (upsampled_image - upsampled_image.min())
            / (upsampled_image.max() - upsampled_image.min())
            * 255
        )
        upsampled_image = upsampled_image.astype(np.uint8)

    return upsampled_image


def imread(path) -> np.ndarray:
    """Read an image from a given path.
    Args:
        path (Path): Path to the image file.

    Returns:
        np.ndarray: Image data as a numpy array.
    """
    if not isinstance(path, Path):
        path = Path(path)
    if not path.is_file():
        raise ValueError("Path is not a file.")

    if path.suffix.lower() == ".tif" or path.suffix.lower() == ".tiff":
        return tiff.imread(path)
    elif path.suffix.lower() in [".png", ".jpg", ".jpeg", ".bmp", ".gif"]:
        with Image.open(path) as img:
            return np.array(img)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")


def normalize(img, percent_saturation=0.002):
    int_cut = np.percentile(img, 100 - percent_saturation)
    normalized = img / int_cut
    img_norm = np.clip(normalized, 0, 1)
    img_norm_bit = img_norm * 65535
    return img_norm_bit.astype(np.uint16)


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


def get_stars(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def get_representative_ids(df, feature, rep_min=40, rep_max=60):
    """
    For a given feature returns the IDs of samples that are representative all conditions, meaning they are around the median of the distribution for that feature. This can be used to select samples for visualizations that are representative of the overall trends in the data.
    params:
    df: DataFrame containing results of an expeiment
    feature: Feature, meaning column, of the dataframe in the
    rep_min: lower percentile threshold for selecting representative samples
    rep_max: upper percentile threshold for selecting representative samples
    """
    condition_dict = {}
    for condition in df["condition"].unique():
        representative_ids = []
        for rna_type in df["construct_id"].unique():
            sub_df = df.query("condition == @condition")
            min_val = np.percentile(sub_df[feature].values, rep_min)
            max_val = np.percentile(sub_df[feature].values, rep_max)
            if len(representative_ids) == 0:
                # print(sub_df.head(10))
                # print(sub_df.columns)
                representative_ids = sub_df.query(f"@min_val <= {feature} <= @max_val")[
                    "uid"
                ].values
            else:
                representative_ids = np.intersect1d(
                    representative_ids,
                    sub_df.query(f"@min_val <= {feature} <= @max_val")["uid"].values,
                )

            condition_dict[condition] = representative_ids

    return condition_dict
