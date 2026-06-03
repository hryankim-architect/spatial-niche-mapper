# What is out of scope

Explicit boundaries for this repo. Small scope kept intentionally.

## Hard boundaries

- **No proprietary material and no patient data.** The tissue, reference
  signatures, and niches are synthetic and deterministically generated.
- **A clean-room implementation of public methods.** NNLS deconvolution and
  neighbourhood niche clustering are the transparent baselines behind RCTD,
  SpatialDWLS, cell2location, and squidpy; this repo implements them from
  scratch and does not vendor those packages.

## Default out-of-scope items

- **Probabilistic deconvolution** (negative-binomial noise models, variational
  inference à la cell2location). NNLS is the legible baseline used here.
- **Real spatial file formats** (AnnData h5ad, Zarr, 10x SpaceRanger output) and
  the `[spatial]` stack (scanpy / squidpy), documented as the real-data path,
  not exercised in the synthetic demo.
- **Image / H&E registration.** Spot coordinates are used; histology image
  alignment is out of scope.
- **Statistical-power / benchmark claims.** The reported correlation, MAE, and
  ARI describe the synthetic demo only; they are an existence proof that the
  method runs end to end, not a real-cohort benchmark.
- **Single-cell reference construction.** The reference is synthetic; on real
  data it would come from a matched annotated scRNA-seq atlas
  (`sc-tumor-annotator` is the sibling repo for that side).

## How to add an item

File an issue first. Then open a PR that appends the item to the relevant
section above, states the reason in one sentence, and links back to that issue.
Scope changes go through review the same way code changes do.
