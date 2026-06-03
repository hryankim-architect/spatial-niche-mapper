# Architecture

The pipeline lives in a single Python process. Calling `make run` or
`scripts/run_lab.sh` invokes `spatialniche.pipeline.run_pipeline`, which is the
only public entry point. Everything else — data synthesis, deconvolution, niche
detection, scoring, and side-channel hooks — is called from there in a fixed
order.

## Method modules

`synth.py` builds a deterministic Visium-style tissue section. It lays out a
spot grid with three concentric ground-truth regions (tumor core, immune margin,
stroma), assigns per-spot cell-type proportions drawn from those regions, builds
a single-cell reference signature matrix, and produces proportion-mixed,
Poisson-noised spot expression. Nothing is random: given the same seed, the
section is identical.

`deconv.py` runs per-spot deconvolution by non-negative least squares against
the reference matrix. Each spot's expression vector is decomposed into cell-type
proportions, renormalized to a simplex. Accuracy is reported as per-cell-type
Pearson correlation and mean absolute error against the ground-truth proportions.

`niche.py` builds a spatial k-NN graph over the spot coordinates, summarizes
each spot by the mean cell-type composition of its neighbourhood, and clusters
those neighbourhood vectors with KMeans. The recovered clusters are compared to
the ground-truth concentric regions using adjusted Rand index.

## Why NNLS deconvolution

Each Visium spot aggregates signal from a mixture of cell types. Its expression
profile is therefore approximately the proportion-weighted sum of the cell-type
reference signatures. Recovering those proportions is a non-negative
least-squares problem. RCTD, SpatialDWLS, and cell2location all build on the
same decomposition but add probabilistic noise models and regularization that
are not needed here. NNLS makes the method legible and deterministic without
hiding the core idea behind a generative model.

## Why neighbourhood composition for niches

A spatial niche is defined by what surrounds a spot, not by the spot in
isolation. By replacing each spot's own composition with the mean composition of
its k nearest neighbours and then clustering the resulting vectors, the pipeline
recovers recurrent tissue environments. This mirrors the neighbourhood-pattern
approach in squidpy but uses only scikit-learn, so there is no additional spatial
stack to install. Because the ground-truth regions are known, the ARI score is
a direct, reproducible measure of recovery quality.

## Substrate hooks

Three optional channels attach to the pipeline without altering its logic.

The audit channel calls `audit.emit` on every pipeline run. Each call appends a
record to a local NDJSON file. Every record's `prev_hash` field holds the
SHA-256 digest of the previous record's canonical JSON serialization, forming a
hash chain. The local NDJSON file is the durable source of truth regardless of
whether the optional remote sink is reachable. If `AUDIT_HOST` is unset the
remote POST is skipped and only the local file is written. Write throughput for
the hash-chain bookkeeping is around 6.19 µs per entry.

The tracking channel calls `tracking.run` to open an MLflow active run and log
the accuracy metrics. If `MLFLOW_TRACKING_URI` is not set the call becomes a
no-op and the pipeline continues without recording anything to MLflow.

The canary channel, controlled by `SPATIALNICHE_CANARY_FIXTURE`, runs a fast
invariant check that confirms deconvolution recovers the ground-truth proportions
above a minimum correlation floor. When the env var is unset the bundled fixture
is used. The check completes in well under a second.

## What this architecture intentionally avoids

There is no deep generative model, no GPU requirement, no real spatial file
formats such as h5ad or Zarr, and no probabilistic deconvolution. The design
goal is a method that runs end to end on a laptop and produces scores that are
directly interpretable against a known ground truth.
