# FRAP Reanalysis Package

This package processes FRAP timelapse datasets, computes normalized recovery curves, fits recovery models, and exports experiment results to HDF5.

## Main modules

- `frappy.py`: high-level workflow classes:
  - `FrapSample`: loads one sample (image, masks, metadata), computes intensities and fit.
  - `FrapExperiment`: loads multiple samples from one experiment directory, computes composed/averaged data and experiment-level fit.
- `ios.py`: input/output helpers for reading microscopy files and writing HDF5 groups.
- `calculator.py`: normalization routines (`phair_double_normalization`).
- `fitting.py`: curve models and fitting (`recovery_fit`).

## Expected input folder layout

For one sample folder:

- `*.lsm` (source metadata)
- `*.tif` (timelapse image sequence)
- `masks/` (or another subfolder containing `"mask"` in the name), with:
  - `*irrad*.tif`
  - `*correction*.tif`
  - `*background*.tif`

For one experiment folder:

- one subfolder per sample

## Quick usage

```python
from pathlib import Path
from frappy import FrapExperiment

experiment_root = Path("/path/to/experiment_folder")
fe = FrapExperiment(experiment_root)

# Build composed tables and experiment-level fit
fe.process_sample_data()

# Optional visualization
fe.generate_report()
```

## HDF5 export

Export one experiment into a shared HDF5 file:

```python
from pathlib import Path
from frappy import FrapExperiment

fe = FrapExperiment(Path("/path/to/experiment_folder"))
fe.export_to_hdf(
    hdf_path=Path("/path/to/all_frap_experiments.h5"),
    experiment_name="MyExperimentName",
)
```

If `experiment_name` is omitted, the experiment folder name is used.

### HDF5 layout

```text
/<experiment_name>
  /meta_data
  /FrapSamples
    /<sample_name>
      /images
        image_sequence
        mask_irradiated
        mask_background
        mask_reference
      /meta_data
      /intensity_data
        time
        index
        normalized_intensity
        intensity_irradiated_raw
        intensity_reference_raw
        intensity_background_raw
      /fit_data
        fit_type
        /fit_params
        /fit_values
  /ComposedData
    /Data
      /df
      /df_averaged
```

`df` stores one row per timepoint per sample with columns:

- `sample_name`
- `time`
- `index`
- `normalized_intensity`
- `intensity_irradiated_raw`
- `intensity_reference_raw`
- `intensity_background_raw`

`df_averaged` stores the experiment-averaged numeric time series grouped by `index`.

## Dependencies

Core dependencies used by this codebase:

- `numpy`
- `pandas`
- `tifffile`
- `h5py`
- `scipy`
- `scikit-learn`
- `matplotlib`
- `seaborn`
