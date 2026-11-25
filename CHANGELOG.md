# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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

[Unreleased]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/VuVietAnh2309/fusion-metrics/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/VuVietAnh2309/fusion-metrics/releases/tag/v1.0.0
