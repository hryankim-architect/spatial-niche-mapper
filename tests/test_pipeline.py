"""End-to-end pipeline + audit-chain smoke tests."""

from __future__ import annotations

import json
from pathlib import Path

from spatialniche import audit, pipeline


def test_pipeline_runs_and_scores(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    monkeypatch.delenv("MLFLOW_TRACKING_URI", raising=False)

    result = pipeline.run_pipeline("smoke", tmp_path / "artifacts", side=18, seed=0)
    m = result["metrics"]
    assert m["deconv_mean_celltype_corr"] > 0.8
    assert m["niche_recovery_ari"] > 0.7

    payload = json.loads(Path(result["artifact_path"]).read_text())
    assert "deconvolution" in payload and "niche_recovery_ari" in payload


def test_audit_chain_valid_after_pipeline(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    pipeline.run_pipeline("smoke", tmp_path / "artifacts", side=14, seed=0)
    ok, n, first_bad = audit.verify()
    assert ok, f"audit chain invalid at {first_bad}"
    assert n >= 2
