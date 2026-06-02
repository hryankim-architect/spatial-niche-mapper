# Architecture

One Python process, three method modules, three substrate hooks, the same
house style as the rest of the capability-portrait portfolio.

## Control flow

```
                make run / scripts/run_lab.sh
                          │
                          ▼
              spatialniche.pipeline.run_pipeline
                          │
        ┌─────────────────┼──────────────────────────────┐
        ▼                 ▼                               ▼
  audit.emit         tracking.run                       body
 (NDJSON +         (MLflow active run,        synth.generate → deconv.deconvolve
  optional POST)    no-op if unset)            → niche.detect_niches → score
                          │
                          ▼
              artifacts/<name>.json  (deconv accuracy + niche ARI)
```

## Method modules

| Module | Responsibility |
|---|---|
| `synth.py` | Deterministic Visium-style section: a spot grid with concentric ground-truth niches (tumor core / immune margin / stroma), per-spot ground-truth cell-type proportions, a single-cell reference signature matrix, and proportion-mixed + Poisson-noised spot expression. |
| `deconv.py` | Per-spot deconvolution by non-negative least squares against the reference, renormalized to a proportion simplex; accuracy scored by per-cell-type correlation and MAE vs the ground truth. |
| `niche.py` | Spatial k-NN graph over coordinates → per-spot neighbourhood cell-type composition → KMeans niches; recovery scored by adjusted Rand index vs the ground-truth regions. |

## Why NNLS deconvolution

Each spot is a mixture of cell types, so its expression is approximately the
proportion-weighted sum of the cell-type reference signatures. Recovering the
proportions is therefore a non-negative least-squares problem, the transparent
public baseline behind RCTD / SpatialDWLS / cell2location (which add
probabilistic noise models and regularization). NNLS keeps the demo legible,
deterministic, and dependency-light while demonstrating the method.

## Why neighbourhood composition for niches

A niche is defined by what is *around* a spot, not just the spot itself.
Summarizing each spot by its spatial-neighbourhood mean composition and
clustering those vectors recovers recurrent tissue niches, the squidpy-style
neighborhood pattern, implemented with scikit-learn. Recovery is scored against
the synthetic ground-truth regions (ARI), so the claim is measurable.

## Substrate integration

| Channel | Module | Env var | Behaviour when unset |
|---|---|---|---|
| Audit | `audit` | `AUDIT_HOST` | local NDJSON only (source of truth) |
| MLflow | `tracking` | `MLFLOW_TRACKING_URI` | no-op |
| Canary | `canary` | `SPATIALNICHE_CANARY_FIXTURE` | uses the bundled fixture |

The canary asserts the core invariant (deconvolution recovers the ground-truth
proportions above a floor correlation) in well under a second.

## What this architecture intentionally avoids

No deep generative model, no GPU, no real spatial file formats (h5ad / Zarr),
and no probabilistic deconvolution. The point is the method and its measurable
invariants, runnable on a laptop.
