"""Spatial niche-detection invariants."""

from __future__ import annotations

from spatialniche import niche, synth


def test_neighborhood_composition_shape() -> None:
    d = synth.generate(side=12, seed=1)
    comp = niche.neighborhood_composition(d.proportions, d.coords, k=6)
    assert comp.shape == d.proportions.shape


def test_niches_recover_ground_truth_regions() -> None:
    d = synth.generate(side=22, seed=0)
    pred = niche.detect_niches(d.proportions, d.coords, n_niches=3)
    ari = niche.niche_recovery_ari(d.niche_labels, pred)
    # detected niches align with the ground-truth concentric regions
    assert ari > 0.7, ari


def test_niche_detection_is_deterministic() -> None:
    d = synth.generate(side=14, seed=2)
    a = niche.detect_niches(d.proportions, d.coords, n_niches=3, seed=0)
    b = niche.detect_niches(d.proportions, d.coords, n_niches=3, seed=0)
    assert niche.niche_recovery_ari(a, b) == 1.0
