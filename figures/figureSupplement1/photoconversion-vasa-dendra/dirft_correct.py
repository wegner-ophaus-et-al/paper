from pystackreg import StackReg
import numpy as np
import tifffile as tiff


sr = StackReg(StackReg.RIGID_BODY)
rel_path = "data/p11_directHit_bigGranule_same_cell_as_p10.lsm"
img = tiff.imread(rel_path)

img_uc = img[:, 0]
img_conv = img[:, 1]

reg = sr.register_stack(img_uc)
img_uc_dc = sr.transform_stack(img_uc)
img_conv_dc = sr.transform_stack(img_conv)

tiff.imwrite(
    "data/p11_directHit_bigGranule_same_cell_as_p10_dc.tif",
    np.array([img_uc_dc, img_conv_dc]).astype("float32"),
    imagej=True,
)
