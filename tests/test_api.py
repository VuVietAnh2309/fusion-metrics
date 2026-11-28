"""Behavioural tests of the public API (numerical agreement with MATLAB is in verify/selftest.py)."""
import warnings

import numpy as np
import pytest

import fusion_metrics as fm


def _images(seed=0, n=64):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    A = (127 + 100 * np.sin(xx / 5.0) + rng.normal(0, 5, (n, n))).clip(0, 255).astype(np.uint8)
    B = (127 + 100 * np.cos(yy / 7.0) + rng.normal(0, 5, (n, n))).clip(0, 255).astype(np.uint8)
    F = ((A.astype(float) + B) / 2).round().astype(np.uint8)
    return A, B, F


def test_presets_reference_known_metrics():
    for name, metrics in fm.PRESETS.items():
        assert metrics, name
        assert set(metrics) <= set(fm.METRICS), name


def test_ivif_preset_is_finite():
    A, B, F = _images()
    with warnings.catch_warnings():
        warnings.simplefilter("error", fm.ScaleMismatchWarning)
        out = fm.evaluate(A, B, F, "ivif")
    assert list(out) == fm.PRESETS["ivif"]
    assert all(np.isfinite(v) for v in out.values())


def test_mixed_range_warns_and_collapses_qabf():
    A, B, F = _images()
    good = fm.Q_ABF(A, B, F)
    with pytest.warns(fm.ScaleMismatchWarning, match="F"):
        bad = fm.evaluate(A, B, F / 255.0, ["Q_ABF"])["Q_ABF"]
    assert bad < 0.1 * good          # the silent failure the warning is there for


def test_all_normalised_warns():
    A, B, F = (x / 255.0 for x in _images())
    with pytest.warns(fm.ScaleMismatchWarning, match="all images"):
        fm.evaluate(A, B, F, ["Q_En"])


def test_check_can_be_disabled():
    A, B, F = _images()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fm.evaluate(A, B, F / 255.0, ["Q_En"], check=False)


def test_check_range_ok():
    assert fm.check_range(*_images())


def test_fused_equal_to_source():
    A, B, _ = _images()
    # identical inputs: SSIM-type metrics are 1, entropy is unchanged
    assert fm.MS_SSIM(A, A, A) == pytest.approx(1.0)
    assert fm.Q_En(A) == pytest.approx(fm.Q_En(A.astype(np.float64)))


def test_unknown_metric_raises():
    A, B, F = _images()
    with pytest.raises(KeyError):
        fm.evaluate(A, B, F, ["NOT_A_METRIC"])
