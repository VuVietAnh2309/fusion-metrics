# Contributing to fusion-metrics

Thanks for helping make image-fusion evaluation reproducible! This project has one rule above all others:

> **A metric is only accepted when it reproduces its reference MATLAB implementation value by value.**

## Getting help

- **Questions / usage help** — open a [Discussion](https://github.com/VuVietAnh2309/fusion-metrics/discussions)
  (Q&A category). Please include your `fusion-metrics` version (`fusion-metrics --version`), Python version and a
  minimal example.
- **Bug reports** — open an [issue](https://github.com/VuVietAnh2309/fusion-metrics/issues/new/choose) with the
  *Bug report* template. If a value differs from MATLAB, attach (or link) the three images and the MATLAB value.
- **Security or private matters** — email vietanhresearcher@gmail.com.

## Development setup

```bash
git clone https://github.com/VuVietAnh2309/fusion-metrics.git
cd fusion-metrics
python -m venv .venv && source .venv/bin/activate
pip install -e .
python verify/selftest.py        # must print "ALL OK"
```

## Porting a new metric

The [roadmap issues](https://github.com/VuVietAnh2309/fusion-metrics/issues?q=is%3Aissue+is%3Aopen+label%3Ametric-port)
list the `.m` files in `matlab/` that still need a Python port. Comment on one before you start so work isn't duplicated.

1. **Port line by line.** Keep the `.m` file's function name and argument order (`Q_X(A, B, F)`, or `Q_X(F)` for
   no-reference metrics). Reproduce MATLAB semantics with the helpers in `fusion_metrics/_matlab.py`
   (`conv2_same`, `filter2_valid`, `fspecial_gaussian`, `mround`, `to_uint8`, …) rather than "equivalent" NumPy/SciPy
   calls — they often differ at borders or for even-sized kernels.
2. **Register it** in `METRICS` (and, if it fits, a `PRESETS` group) in `fusion_metrics/__init__.py`, and add a row to
   the metric table in `README.md` with the original paper.
3. **Add reference values.** With GNU Octave (see *Verification* in the README):

   ```bash
   cd verify
   python make_reference_cases.py && octave-cli run_reference.m \
     && python selftest.py --export && python selftest.py
   ```

   Add the metric to the `names` list and the computation loop in `run_reference.m`; `selftest.py` reads the names
   from the reference file and calls `METRICS[name]`, so no change is needed there. The default tolerance is
   `1e-9 × max(1, |reference|)`; any relaxed tolerance must be justified in a comment (see `RELAXED` in `selftest.py`).
4. **Document MATLAB-only behaviour.** If the `.m` file relies on something Octave does not reproduce (a toolbox
   function, a mex file), add a shim in `verify/octave_shims/` and explain it in the PR.
5. **Licensing.** If the `.m` file carries its own copyright notice, keep it, and add the port to `NOTICE`.

## Pull requests

- One metric (or one fix) per PR; reference the issue it closes.
- CI must be green (`verify/selftest.py` on Python 3.10–3.13).
- Add a line to `CHANGELOG.md` under *Unreleased*.
- Code style: follow the surrounding code; comments should explain *why* a line mirrors MATLAB, not *what* it does.

## Releases (maintainer)

Bump `version` in `pyproject.toml`, `__version__` in `fusion_metrics/__init__.py` and `version`/`date-released` in
`CITATION.cff`, move *Unreleased* in `CHANGELOG.md` under the new version, then publish a GitHub Release `vX.Y.Z`.
PyPI (Trusted Publishing) and Zenodo (DOI) update automatically.

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
