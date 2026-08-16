# Photoconversion of germ granules

Analysis of a photoconversion  experiment in a primordial
germ cell of an early (10 hpf) zebrafish embryo. A region containing a germ
granule is photoconverted, and the fluorescence intensity of the converted and
unconverted channels is tracked over time in the photoconverted granule and in
neighboring, non-converted granules to assess exchange/recovery behavior.

## Input data

- `data/p11_directHit_bigGranule_same_cell_as_p10.lsm` — raw two-channel
  (unconverted / converted) Zeiss LSM confocal time-lapse, including the
  photoconversion (bleach) event metadata used to locate the conversion
  window in time.
- `data/masks/granule_drift_corrected.tif` — a labeled mask identifying
  individual granules (by pixel label) in the drift-corrected image; origin
  of the segmentation is not otherwise documented in this repo.
- All of the above, plus the derived/registered stacks
  (`*_dc.tif`, `*_dc_imagej.tif`, `*.corrected.tif`), are committed directly
  under `data/`.

## Pipeline

The primary pipeline is two scripts run in order:

1. `dirft_correct.py` — registers the two-channel stack with `pystackreg`
   (rigid body) against the unconverted channel, and writes the registered
   stack to `data/p11_directHit_bigGranule_same_cell_as_p10_dc.tif`.
   (`data/p11_directHit_bigGranule_same_cell_as_p10_dc_imagej.tif`, consumed
   by the next step, was produced the same way / is the ImageJ-compatible
   version of this output.)
2. `measure_and_plot.py` — loads the registered stack and the granule mask,
   locates the photoconversion window from the LSM metadata
   (`utils.bleach_times`), measures per-granule mean/summed intensity in
   both channels over time, and produces the figure panel: example frames
   (unconverted and converted channel, with granule outlines and a scale
   bar) alongside intensity-over-time traces for the photoconverted granule
   and the other (recovering) granules. Saves `out.pdf`.

`utils.py` holds shared helpers (scale bar drawing, pixel size and bleach
timing from LSM metadata).

`process.py` is a separate, legacy pipeline that does its own
drift correction (`phase_cross_correlation`) and *automated* granule
detection/tracking (DoG blob detection + Hungarian-algorithm linking across
frames) instead of using the manual mask. It writes per-spot measurements to
`results.csv` (the version currently committed in this repo) and shows an
exploratory plot. It is not the pipeline used to produce the figure.

## Outputs

- `data/p11_directHit_bigGranule_same_cell_as_p10_dc.tif` — registered stack
  (from `dirft_correct.py`).
- `out.pdf` — the figure panel (example frames + intensity traces),
  produced by `measure_and_plot.py`.

***
This README was generated using Claude Code (Opus 5) and proofread by JW
