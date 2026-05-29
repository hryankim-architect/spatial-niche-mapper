"""Spot deconvolution by non-negative least squares.

Each spatial spot is a mixture of cell types. Given a single-cell reference
signature matrix (cell type x genes), deconvolution estimates the per-spot
cell-type proportions. We solve, per spot, the non-negative least-squares
problem

    minimize || R^T p - x ||_2   subject to  p >= 0,

then renormalize ``p`` to sum to 1. This is the legible public method behind
reference-based spatial deconvolution (RCTD / SpatialDWLS / cell2location use
richer probabilistic variants; NNLS is the transparent baseline).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import nnls


def deconvolve(expr: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Estimate per-spot cell-type proportions.

    Parameters
    ----------
    expr:
        spots x genes expression.
    reference:
        cell_types x genes signature matrix.

    Returns
    -------
    proportions : np.ndarray
        spots x cell_types, each row non-negative and summing to 1.
    """
    a = reference.T  # genes x cell_types
    n_spots = expr.shape[0]
    k = reference.shape[0]
    out = np.zeros((n_spots, k))
    for i in range(n_spots):
        coef, _ = nnls(a, expr[i])
        total = coef.sum()
        out[i] = coef / total if total > 0 else np.full(k, 1.0 / k)
    return out


def per_celltype_correlation(true_prop: np.ndarray, est_prop: np.ndarray) -> np.ndarray:
    """Pearson correlation between true and estimated proportions, per cell type."""
    k = true_prop.shape[1]
    cors = np.zeros(k)
    for j in range(k):
        t, e = true_prop[:, j], est_prop[:, j]
        if t.std() < 1e-12 or e.std() < 1e-12:
            cors[j] = 0.0
        else:
            cors[j] = float(np.corrcoef(t, e)[0, 1])
    return cors


def mean_absolute_error(true_prop: np.ndarray, est_prop: np.ndarray) -> float:
    """Mean absolute error across all spot x cell-type proportion entries."""
    return float(np.abs(true_prop - est_prop).mean())
