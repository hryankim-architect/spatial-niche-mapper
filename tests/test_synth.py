"""Synthetic spatial generator invariants."""

from __future__ import annotations

import numpy as np

from spatialniche import synth


def test_generate_is_deterministic() -> None:
    a = synth.generate(side=10, seed=5)
    b = synth.generate(side=10, seed=5)
    assert np.array_equal(a.expr, b.expr)
    assert np.array_equal(a.proportions, b.proportions)


def test_shapes_and_proportions_sum_to_one() -> None:
    d = synth.generate(side=12, seed=1)
    assert d.n_spots == 144
    assert d.proportions.shape[1] == len(synth.CELL_TYPES)
    assert np.allclose(d.proportions.sum(axis=1), 1.0)
    assert d.reference.shape[0] == len(synth.CELL_TYPES)


def test_three_ground_truth_niches_present() -> None:
    d = synth.generate(side=20, seed=2)
    assert set(d.niche_labels) == set(synth.NICHES)


def test_tumor_core_is_tumor_enriched() -> None:
    d = synth.generate(side=20, seed=2)
    ti = synth.CELL_TYPES.index("Tumor")
    core = d.proportions[d.niche_labels == "tumor_core", ti].mean()
    stroma = d.proportions[d.niche_labels == "stroma", ti].mean()
    assert core > stroma
