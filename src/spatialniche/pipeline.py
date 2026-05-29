"""End-to-end pipeline: synthetic tissue -> deconvolution -> niche map.

House-style shape: audit_start -> tracking_start -> body -> tracking_end ->
audit_end. The body generates a synthetic Visium-style section, deconvolves each
spot against the single-cell reference, detects spatial niches, and scores both
against the known ground truth. Deterministic given the seed.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import click

from spatialniche import audit, deconv, niche, synth, tracking
from spatialniche.synth import CELL_TYPES


def _run_id(name: str) -> str:
    return f"{name}-{datetime.now(UTC).strftime('%Y%m%d-%H%M%S')}"


def run_pipeline(
    run_name: str,
    out_dir: Path,
    *,
    side: int = 24,
    seed: int = 0,
) -> dict[str, Any]:
    """Generate -> deconvolve -> niche-detect -> score -> artifact."""
    out_dir.mkdir(parents=True, exist_ok=True)
    job_id = _run_id(run_name)

    audit.emit(action="pipeline_start", job_id=job_id, fields={"side": side, "seed": seed})

    metrics: dict[str, float] = {}
    with tracking.run(name=job_id, experiment="spatialniche"):
        data = synth.generate(side=side, seed=seed)

        est = deconv.deconvolve(data.expr, data.reference)
        cors = deconv.per_celltype_correlation(data.proportions, est)
        mae = deconv.mean_absolute_error(data.proportions, est)

        pred_niches = niche.detect_niches(est, data.coords, n_niches=len(set(data.niche_labels)))
        ari = niche.niche_recovery_ari(data.niche_labels, pred_niches)

        metrics["deconv_mean_celltype_corr"] = float(cors.mean())
        metrics["deconv_mae"] = mae
        metrics["niche_recovery_ari"] = ari
        tracking.log_metrics(metrics)

    artifact = {
        "run_name": run_name,
        "job_id": job_id,
        "n_spots": int(data.n_spots),
        "n_genes": len(data.genes),
        "deconvolution": {
            "per_celltype_corr": {c: float(v) for c, v in zip(CELL_TYPES, cors, strict=False)},
            "mean_celltype_corr": float(cors.mean()),
            "mae": mae,
        },
        "niche_recovery_ari": ari,
    }
    artifact_path = out_dir / f"{run_name}.json"
    with artifact_path.open("w", encoding="utf-8") as fh:
        json.dump(artifact, fh, indent=2, sort_keys=True)

    audit.emit(action="pipeline_end", job_id=job_id, fields={
        "artifact_path": str(artifact_path),
        "mean_celltype_corr": metrics["deconv_mean_celltype_corr"],
        "niche_recovery_ari": ari,
    })
    return {"job_id": job_id, "artifact_path": str(artifact_path), "metrics": metrics}


@click.group()
def cli() -> None:
    """spatialniche spatial-transcriptomics niche-mapping pipeline (synthetic POC)."""


@cli.command()
@click.option("--manifest", type=click.Path(path_type=Path), default=Path("data/manifest.yaml"))
@click.option("--out", type=click.Path(file_okay=False, path_type=Path), default=Path("data"))
def fetch(manifest: Path, out: Path) -> None:
    """No-op for the synthetic demo: tissue is generated, not downloaded."""
    click.echo(json.dumps(
        {"status": "synthetic-demo", "note": "spatial tissue is generated deterministically; "
         "see data/manifest.yaml for the public datasets the method targets",
         "manifest": str(manifest), "out": str(out)}, indent=2))


@cli.command()
@click.option("--name", default="demo")
@click.option("--out", type=click.Path(file_okay=False, path_type=Path), default=Path("artifacts"))
@click.option("--side", default=24, type=int)
@click.option("--seed", default=0, type=int)
def run(name: str, out: Path, side: int, seed: int) -> None:
    """Run the end-to-end pipeline."""
    result = run_pipeline(name, out, side=side, seed=seed)
    click.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    cli()
