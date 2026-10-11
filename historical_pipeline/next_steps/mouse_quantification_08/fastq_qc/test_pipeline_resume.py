"""Isolated synthetic tests of parent run_once; no production artifacts mutated."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
spec = importlib.util.spec_from_file_location("isolated_common08", PARENT / "common08.py")
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
fixture = Path(tempfile.mkdtemp(prefix="synthetic_resume_", dir=HERE)).resolve()
common.HERE = fixture / "metadata"
common.DATA = fixture / "data"
input_path = fixture / "input.txt"
output_path = fixture / "output.txt"
counter = fixture / "executions.txt"
input_path.write_text("AAAA", encoding="utf-8")
checks = []


def check(name, result, detail=None):
    checks.append({"name": name, "passed": bool(result), "detail": detail})


def refusal(name, action):
    try:
        action()
    except (AssertionError, OSError, RuntimeError) as exc:
        check(name, True, {"type": type(exc).__name__, "message": str(exc)})
    else:
        check(name, False, "Unexpectedly accepted")


def fingerprints():
    return {"input_sha256": common.digest(input_path), "qc_script_sha256": "fixed-fixture-script", "executable_sha256": "fixed-fixture-tool", "index_sha256": "fixed-fixture-index"}


child = "from pathlib import Path; import sys; p,c=map(Path,sys.argv[1:3]); p.write_text('1111'); c.write_text(str(int(c.read_text())+1 if c.exists() else 1))"
argv = [sys.executable, "-c", child, str(output_path), str(counter)]
first = common.run_once("fixture_success", argv, input_fingerprints=fingerprints(), output_paths=[output_path])
check("First actual synthetic command passed", first["status"] == "passed" and counter.read_text() == "1")
check("First record stores input identities", first["input_fingerprints"] == fingerprints())
check("First record stores output digest", first["output_sha256"][str(output_path)] == common.digest(output_path))
second = common.run_once("fixture_success", argv, input_fingerprints=fingerprints(), output_paths=[output_path])
check("Unchanged resume reuses command without child rerun", counter.read_text() == "1" and first == second)
input_path.write_text("BBBB", encoding="utf-8")
refusal("Same-size changed input fingerprint refused", lambda: common.run_once("fixture_success", argv, input_fingerprints=fingerprints(), output_paths=[output_path]))
input_path.write_text("AAAA", encoding="utf-8")
for key in ("qc_script_sha256", "executable_sha256", "index_sha256"):
    altered = fingerprints()
    altered[key] = "changed-fixture-value"
    refusal("Changed " + key + " refused", lambda altered=altered: common.run_once("fixture_success", argv, input_fingerprints=altered, output_paths=[output_path]))
output_path.write_text("2222", encoding="utf-8")
refusal("Same-size changed output refused", lambda: common.run_once("fixture_success", argv, input_fingerprints=fingerprints(), output_paths=[output_path]))
output_path.write_text("1111", encoding="utf-8")
missing_backup = fixture / "output_moved_for_missing_test.txt"
output_path.rename(missing_backup)
refusal("Missing saved output refused", lambda: common.run_once("fixture_success", argv, input_fingerprints=fingerprints(), output_paths=[output_path]))
missing_backup.rename(output_path)
refusal("Changed argv refused", lambda: common.run_once("fixture_success", argv + ["different"], input_fingerprints=fingerprints(), output_paths=[output_path]))
check("Refused reuse never reran child", counter.read_text() == "1")

# An isolated legacy record deliberately lacks provenance; it must not be trusted.
legacy = dict(first)
legacy["key"] = "fixture_legacy"
legacy.pop("input_fingerprints")
legacy.pop("output_sha256")
common.write_json(common.HERE / "execution/fixture_legacy.json", legacy)
refusal("Absent input provenance refused", lambda: common.run_once("fixture_legacy", argv, input_fingerprints=fingerprints(), output_paths=[output_path]))

failed_argv = [sys.executable, "-c", "raise SystemExit(7)"]
refusal("Actual failed child reported", lambda: common.run_once("fixture_failure", failed_argv, input_fingerprints=fingerprints(), output_paths=[]))
failed_record = json.loads((common.HERE / "execution/fixture_failure.json").read_text(encoding="utf-8"))
check("Failed child exit retained", failed_record["status"] == "failed" and failed_record["exit_code"] == 7)
refusal("Failed command cannot be reused as passing", lambda: common.run_once("fixture_failure", failed_argv, input_fingerprints=fingerprints(), output_paths=[]))

report = {
    "scope": "Isolated synthetic tests of run_once fingerprint-based reuse; not study-data execution",
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "fixture_directory": str(fixture),
    "parent_globals_overridden_only_in_test_process": ["HERE", "DATA"],
    "interpreter": sys.executable,
    "reviewed_sha256": {name: hashlib.sha256((PARENT / name).read_bytes()).hexdigest() for name in ("common08.py", "run_study.py")},
    "passed": all(x["passed"] for x in checks),
    "total_checks": len(checks),
    "checks": checks,
}
(HERE / "pipeline_resume_test_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"passed": report["passed"], "checks": len(checks), "fingerprints": report["reviewed_sha256"]}))
raise SystemExit(0 if report["passed"] else 1)
