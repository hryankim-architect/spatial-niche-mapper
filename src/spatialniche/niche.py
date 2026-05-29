"""Spatial niche detection from per-spot cell-type composition.

A tissue *niche* is a recurrent local cell-type neighborhood. Given per-spot
cell-type proportions and spot coordinates, we:

1. build a spatial k-nearest-neighbour graph over the coordinates,
2. summarize each spot by the *mean cell-type composition of its neighbourhood*
   (the spot plus its neighbours),
3. cluster those neighbourhood vectors with KMeans into ``n_niches`` niches.

This is the standard neighborhood-enrichment / niche-clustering pattern used in
spatial analysis (squidpy-style), implemented transparently with scikit-learn.
Recovery is scored against the synthetic ground-truth regions with the adjusted
Rand index.
"""

from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


def neighborhood_composition(
    proportions: np.ndarray, coords: np.ndarray, *, k: int = 6
) -> np.ndarray:
    """Mean cell-type composition over each spot's spatial k-NN neighbourhood."""
    n = coords.shape[0]
    k_eff = min(k + 1, n)  # +1 because the spot itself is its nearest neighbour
    nn = NearestNeighbors(n_neighbors=k_eff).fit(coords)
    _, idx = nn.kneighbors(coords)
    out = np.zeros_like(proportions)
    for i in range(n):
        out[i] = proportions[idx[i]].mean(axis=0)
    return out


def detect_niches(
    proportions: np.ndarray,
    coords: np.ndarray,
    *,
    n_niches: int = 3,
    k: int = 6,
    seed: int = 0,
) -> np.ndarray:
    """Cluster spots into niches by neighbourhood composition. Returns labels."""
    comp = neighborhood_composition(proportions, coords, k=k)
    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=10)
    return km.fit_predict(comp)


def niche_recovery_ari(true_labels: np.ndarray, pred_labels: np.ndarray) -> float:
    """Adjusted Rand index between ground-truth regions and detected niches."""
    return float(adjusted_rand_score(np.asarray(true_labels).astype(str), pred_labels))
