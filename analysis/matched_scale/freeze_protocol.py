"""Copy fixed inputs and hash a local prospective protocol; no data arithmetic."""
from pathlib import Path
import hashlib
import json
import shutil
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
WORKSPACE = STAGE.parents[3]
SOURCE = STAGE / "source_data" / "single_cell"
FILES = ["library_category_candidates_all25.csv", "library_candidates_all25.csv",
         "category_sampling_counts.csv", "kitagawa_summary_all25.csv",
         "weighting_sensitivity_all25.csv"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

if (HERE / "protocol_freeze.json").exists():
    raise SystemExit("Protocol already frozen; will not replace the original freeze record.")
(HERE / "inputs").mkdir(parents=True, exist_ok=True)
records = []
for name in FILES:
    original = SOURCE / name
    copied = HERE / "inputs" / name
    shutil.copy2(original, copied)
    assert sha(original) == sha(copied)
    records.append({"source": str(original), "copy": str(copied.relative_to(HERE)),
                    "bytes": copied.stat().st_size, "sha256": sha(copied)})
context_paths = [WORKSPACE / "PFHN_投稿针对性审查_2026-10-10" / "PFHN_审查结论与修改清单.md",
                 WORKSPACE / "投稿收尾版_2026-10-10" / "manuscripts" / "manuscript_en.md"]
freeze = {"status": "frozen_before_new_expression_calculation",
          "frozen_utc": datetime.now(timezone.utc).isoformat(), "protocol_version": "1.0",
          "protocol_sha256": sha(HERE / "protocol.md"), "inputs": records,
          "read_only_context": [{"path": str(p), "sha256": sha(p)} for p in context_paths],
          "no_expression_calculation_in_freeze_step": True}
(HERE / "protocol_freeze.json").write_text(json.dumps(freeze, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps({"frozen_utc": freeze["frozen_utc"], "protocol_sha256": freeze["protocol_sha256"],
                  "input_files": len(records)}, indent=2))
