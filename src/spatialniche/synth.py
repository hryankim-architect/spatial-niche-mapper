"""Deterministic synthetic spatial-transcriptomics generator.

Fabricates a small Visium-style tissue section with a *known* ground truth so
deconvolution and niche detection can be scored offline and reproducibly. No
real spatial data is used.

Design
------
- A square grid of **spots**, each at an (x, y) coordinate.
- Three spatial **regions / niches** laid out across the grid — a tumor core, an
  immune-infiltrated margin, and surrounding stroma — each with a characteristic
  *cell-type composition* (the ground-truth proportions per spot).
- A single-cell **reference signature** matrix (cell type x genes): each cell
  type has a marker-driven expression profile.
- Each spot's expression is the proportion-weighted mix of reference signatures
  plus Poisson-like noise — exactly what a deconvolution method must invert.

Cell types: Tumor, T_cell, B_cell, Myeloid, Fibroblast, Endothelial.
Everything derives from one integer seed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

CELL_TYPES = ["Tumor", "T_cell", "B_cell", "Myeloid", "Fibroblast", "Endothelial"]
NICHES = ["tumor_core", "immune_margin", "stroma"]

# Ground-truth mean cell-type composition per niche (rows sum to 1).
NICHE_COMPOSITION = {
    "tumor_core":    {"Tumor": 0.70, "Myeloid": 0.10, "T_cell": 0.05, "B_cell": 0.02, "Fibroblast": 0.08, "Endothelial": 0.05},
    "immune_margin": {"Tumor": 0.20, "Myeloid": 0.20, "T_cell": 0.35, "B_cell": 0.15, "Fibroblast": 0.05, "Endothelial": 0.05},
    "stroma":        {"Tumor": 0.05, "Myeloid": 0.08, "T_cell": 0.07, "B_cell": 0.05, "Fibroblast": 0.55, "Endothelial": 0.20},
}


@dataclass
class SpatialData:
    expr: np.ndarray            # spots x genes
    coords: np.ndarray          # spots x 2 (x, y)
    proportions: np.ndarray     # spots x cell_types (ground-truth)
    niche_labels: np.ndarray    # spots, ground-truth niche per spot
    reference: np.ndarray       # cell_types x genes signature matrix
    genes: list[str]
    meta: dict = field(default_factory=dict)

    @property
    def n_spots(self) -> int:
        return self.expr.shape[0]


def _build_reference(rng, n_genes: int) -> np.ndarray:
    """Cell-type x gene signature matrix; each type elevates a marker block."""
    k = len(CELL_TYPES)
    ref = rng.gamma(shape=1.0, scale=1.0, size=(k, n_genes)) + 0.1
    block = n_genes // k
    for i in range(k):
        ref[i, i * block:(i + 1) * block] += 6.0  # marker block for this type
    # normalize each signature to a expression profile (sums to 1 across genes)
    ref = ref / ref.sum(axis=1, keepdims=True)
    return ref


def _assign_niche(x: int, y: int, side: int) -> str:
    """Concentric layout: core in the centre, margin around it, stroma outside."""
    cx = cy = (side - 1) / 2.0
    r = np.hypot(x - cx, y - cy) / (side / 2.0)
    if r < 0.35:
        return "tumor_core"
    if r < 0.65:
        return "immune_margin"
    return "stroma"


def generate(
    side: int = 24,
    n_genes: int = 120,
    *,
    seed: int = 0,
    counts_per_spot: float = 400.0,
) -> SpatialData:
    """Generate a `side` x `side` spot grid with ground-truth niches + mixtures."""
    rng = np.random.default_rng(seed)
    genes = [f"g{j:03d}" for j in range(n_genes)]
    reference = _build_reference(rng, n_genes)

    coords = []
    niche_labels = []
    proportions = []
    for x in range(side):
        for y in range(side):
            niche = _assign_niche(x, y, side)
            base = np.array([NICHE_COMPOSITION[niche][c] for c in CELL_TYPES])
            # Dirichlet jitter around the niche mean so spots vary
            prop = rng.dirichlet(base * 40.0 + 0.05)
            coords.append((x, y))
            niche_labels.append(niche)
            proportions.append(prop)

    coords = np.array(coords, dtype=float)
    proportions = np.array(proportions)
    niche_labels = np.array(niche_labels, dtype=object)

    # spot expression = proportion-weighted mixture of reference signatures,
    # scaled to a per-spot count budget, with Poisson noise (log1p-normalized).
    mix = proportions @ reference            # spots x genes (expected fractions)
    lam = mix * counts_per_spot
    counts = rng.poisson(lam).astype(np.float64)
    expr = np.log1p(counts)

    return SpatialData(
        expr=expr, coords=coords, proportions=proportions,
        niche_labels=niche_labels, reference=reference, genes=genes,
        meta={"seed": seed, "side": side, "n_genes": n_genes},
    )


def reference_frame(data: SpatialData) -> pd.DataFrame:
    """Reference signatures as a labeled DataFrame (cell types x genes)."""
    return pd.DataFrame(data.reference, index=CELL_TYPES, columns=data.genes)
