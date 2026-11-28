# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Documentation
- State the measured size of the gray-conversion effect (small) and of scale mismatches (large); link the write-up
  on why fusion metrics differ between implementations. The Chinese README overstated the gray-conversion effect.

## [1.2.0] — 2026-11-27

### Added
- Eleven metrics from the remaining `.m` files, ported line by line and verified against them in GNU Octave
  (7 synthetic and 194 real cases): `Q_TE` (Cvejic Tsallis entropy), `Q_NABF` (fusion artifacts), `Q_P` (Zhao
  phase congruency), `Q_S` / `Q_Sw` / `Q_Se` (Piella, `Q_S.m` with sw = 1, 2, 3), `Q_C` (Cvejic), `Q_Y` (Yang),
  `Q_H` (Hossny), `Q_CV` (Chen–Varshney) and `Q_CB` (Chen–Blum). Every `.m` file in `matlab/` now has a port.
- Presets `perception` (`Q_CV`, `Q_CB`) and `liu2012` (the twelve metrics compared by Liu et al., TPAMI 2012);
  `information`, `feature` and `structural` include the new metrics (`Q_H` in no preset, see Notes).
- MATLAB helpers `filter2_same`, `imfilter_corr_symmetric`, `blkproc`, `freqspace_meshgrid`, `im2double`.
- Octave shims `blkproc.m` and `freqspace.m` for verification, and `verify/octave_exactflat/` (see below).

### Notes
- `Q_H` is lower-is-better: `Q_H.m` returns mean(1 − E) where E is an importance-weighted SSIM (a similarity). The
  value is kept as in the reference file and `Q_H` is not included in the `structural` preset.
- `Q_C` and `Q_Y` evaluate windows in which an image is constant in exact arithmetic: with C = 2e-16 the reference
  files return round-off-dependent values there (MATLAB, Octave and NumPy disagree, by up to ~1e-2 on real images).
- `Q_CB` is reproducible across FFT libraries only to a median 2e-5 (max 5e-3) on real images; see README.

## [1.1.0] — 2026-11-26

### Added
- Six metrics ported line by line from their `.m` files and verified against them (≤ 1.4e-13 on synthetic and on
  194 real cases): `Q_CI` (contrast index = SD), `Q_EI` (edge intensity), `Q_CE` (cross entropy), `Q_NMI`
  (normalized mutual information), `MS_SSIM` (SSIM averaged over the two sources) and `Q_UIQI` (universal image
  quality index). With them the package covers SD and SSIM, completing the metric table most infrared–visible
  fusion papers report (EN, SD, SF, MI, SCD, VIF, Qabf, SSIM).
- Presets `structural` (`MS_SSIM`, `Q_UIQI`) and `ivif` (`Q_En`, `Q_CI`, `SF`, `Q_MI`, `Q_SCD`, `VIF_Sheikh`,
  `Q_ABF`, `MS_SSIM`, i.e. EN, SD, SF, MI, SCD, VIF, Qabf, SSIM);
  `image` and `information` now include the new no-reference and information metrics.
- MATLAB helpers `imfilter_conv_replicate` (`imfilter(..., 'replicate', 'conv')`) and `mtimes_elem` (MATLAB
  class rules for element-wise products, used by `Q_UIQI` on `uint8` input).
- Self-test: the new metrics, plus `Q_UIQI` and `Q_CE` on `float64` and [0, 1] inputs.
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue and pull-request templates.
- API reference with worked examples (`docs/API.md`).

## [1.0.1] — 2026-11-26

### Changed
- First release archived on Zenodo (DOI [10.5281/zenodo.24536175](https://doi.org/10.5281/zenodo.24536175)); no code changes.

## [1.0.0] — 2026-11-25

### Added
- 18 image-fusion quality metrics ported line by line from reference MATLAB code and verified against it
  (≤ 1e-12): `Q_MLI`, `Q_AG`, `Q_En`, `SF`, `Q_SF`, `Q_PSNR`, `Q_CC`, `Q_SCD`, `Q_MI`, `Q_ABF`, `Q_VIFF`, `Q_M`,
  `Q_NCIE`, `Q_TE_Nava`, `Q_FSIM`, `Q_FMI`, `AG_VIFB`, `VIF_Sheikh`.
- Presets `image`, `information`, `feature`, `similarity`, `all`; `evaluate` / `evaluate_files`.
- `fusion-metrics` command-line tool (single image or whole folders → CSV/JSON).
- Self-test against Octave reference values; CI on Python 3.10–3.13; PyPI release via Trusted Publishing.

[Unreleased]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/VuVietAnh2309/fusion-metrics/releases/tag/v1.0.0
