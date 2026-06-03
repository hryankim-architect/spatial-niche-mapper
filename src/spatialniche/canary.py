"""Deterministic canary smoke test.

Probed daily by the ``lab_semantic_check.py`` runner. Contract:
completes in well under 30 s, deterministic, exit 0 on success / non-zero on any
deviation, no external services required.

The check generates a tiny synthetic section, deconvolves the spots, and asserts
the central invariant: estimated cell-type proportions recover the ground-truth
proportions above a floor correlation.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

from spatialniche import audit, deconv, synth, tracking

DEFAULT_FIXTURE = Path("tests/fixtures/canary.json")
EXPECTED_KEYS = {"name", "tier", "min_corr"}


def _load_fixture(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"canary fixture not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def check() -> dict[str, Any]:
    fixture_path = Path(os.environ.get("SPATIALNICHE_CANARY_FIXTURE", str(DEFAULT_FIXTURE)))
    fixture = _load_fixture(fixture_path)
    missing = EXPECTED_KEYS - set(fixture.keys())
    if missing:
        return {"ok": False, "reason": f"fixture missing keys: {sorted(missing)}"}

    job_id = f"canary-{fixture['name']}"
    audit.emit(action="canary_start", job_id=job_id, fields={"tier": fixture["tier"]})

    data = synth.generate(side=12, seed=3)
    est = deconv.deconvolve(data.expr, data.reference)
    mean_corr = float(deconv.per_celltype_correlation(data.proportions, est).mean())
    ok = mean_corr >= float(fixture["min_corr"])

    with tracking.run(name=job_id, experiment="canary"):
        tracking.log_metric("deconv_mean_corr", mean_corr)

    audit.emit(action="canary_end", job_id=job_id, fields={"ok": ok, "mean_corr": mean_corr})
    return {"ok": ok, "job_id": job_id, "deconv_mean_corr": mean_corr,
            "min_required": float(fixture["min_corr"])}


def main() -> int:
    result = check()
    print(json.dumps(result, indent=2))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
