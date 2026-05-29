"""Spot deconvolution invariants."""

from __future__ import annotations

import numpy as np

from spatialniche import deconv, synth


def test_proportions_are_valid_simplex() -> None:
    d = synth.generate(side=12, seed=1)
    est = deconv.deconvolve(d.expr, d.reference)
    assert est.shape == d.proportions.shape
    assert (est >= -1e-9).all()
    assert np.allclose(est.sum(axis=1), 1.0)


def test_deconvolution_recovers_truth() -> None:
    d = synth.generate(side=20, seed=0)
    est = deconv.deconvolve(d.expr, d.reference)
    cors = deconv.per_celltype_correlation(d.proportions, est)
    mae = deconv.mean_absolute_error(d.proportions, est)
    # estimated proportions track the ground truth strongly, low absolute error
    assert cors.mean() > 0.8, cors
    assert mae < 0.12, mae


def test_deconvolution_is_deterministic() -> None:
    d = synth.generate(side=10, seed=4)
    a = deconv.deconvolve(d.expr, d.reference)
    b = deconv.deconvolve(d.expr, d.reference)
    assert np.allclose(a, b)
