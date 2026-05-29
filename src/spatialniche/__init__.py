"""spatialniche: a spatial-transcriptomics niche-mapping capability portrait.

A clean-room demonstration of the spatial-transcriptomics analysis pattern:

- **spot deconvolution** — estimate per-spot cell-type proportions from a
  Visium-style spot x gene matrix against a single-cell reference signature,
  by non-negative least squares (the legible public method behind RCTD /
  SpatialDWLS / cell2location-style deconvolution);
- **spatial niche detection** — build a neighborhood graph over spot
  coordinates, summarize each spot's local cell-type composition, and cluster
  recurrent tissue niches (e.g. tumor-core vs immune-infiltrated vs stroma).

Everything runs on synthetic, deterministically-generated data with a known
ground truth. No patient data and no proprietary code or parameters are present.
See the README honest-scope preamble and ``docs/what-is-out-of-scope.md``.
"""

__version__ = "0.1.0"
