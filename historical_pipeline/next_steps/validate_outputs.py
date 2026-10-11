"""Check delivered files, cached source integrity and explicit sample mappings."""
import ast
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(name):
    with (HERE / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    files = list(HERE.glob("*.json")) + list((HERE / "pilot_descriptive").glob("*.json"))
    for path in files:
        load(path)
    scripts = list(HERE.glob("*.py"))
    for path in scripts:
        ast.parse(path.read_text(encoding="utf-8-sig"))
    records = load(HERE / "additional_input_audit.json")["downloads"]
    records += load(HERE / "organoid_metadata_audit.json")["downloads"]
    records += [load(HERE / "platform_download.json")]
    for r in records:
        assert r.get("status", "read") == "read"
        path = HERE / r["local_file"]
        assert path.stat().st_size == r["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == r["sha256"]
    for r in load(HERE.parent / "processed_table_audit.json")["records"]:
        assert hashlib.sha256((HERE.parent / r["local_file"]).read_bytes()).hexdigest() == r["sha256"]
    human = read_csv("human_snRNA_sample_manifest.csv")
    organoid = read_csv("organoid_sample_manifest.csv")
    butyrate = read_csv("butyrate_sample_manifest.csv")
    assert len(human) == len({r["donor_id"] for r in human}) == 48
    assert len(organoid) == len({r["gsm"] for r in organoid}) == 4
    assert len(butyrate) == len({r["gsm"] for r in butyrate}) == 15
    additional = load(HERE / "additional_input_audit.json")
    assert [r["matrix_column_exact"] for r in butyrate] == additional["butyrate"]["count_columns"]
    assert len(additional["vitamin_E_vs_D"]["exact_duplicate_pairs"]) == 3
    assert all(r["exactly_equal_features"] == r["matching_features"] == 19978 and r["rmse"] == 0 for r in additional["vitamin_E_vs_D"]["exact_duplicate_pairs"])
    pilot = load(HERE / "pilot_descriptive/summary.json")
    assert pilot["human_transcript_clusters"] + pilot["excluded_nonhuman_or_control_rows"] == 19978
    assert pilot["unique_symbol_features"] + pilot["ambiguous_features"] + pilot["unmapped_features"] == 19911
    assert len(read_csv("pilot_descriptive/all_feature_descriptives.csv")) == 19911
    assert len(read_csv("pilot_descriptive/literature_panel_samples.csv")) == 98
    assert pilot["script_sha256"] == hashlib.sha256((HERE / "run_descriptive_pilot.py").read_bytes()).hexdigest()
    assert pilot["protocol_sha256"] == hashlib.sha256((HERE / "pilot_protocol.json").read_bytes()).hexdigest()
    links = []
    documents = list(HERE.glob("*.md")) + [HERE.parent / "project_plan.md", HERE.parent.parent / "research_brief.md"]
    for path in documents:
        for link in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if not re.match(r"https?://", link):
                target = (path.parent / link.split("#")[0]).resolve()
                assert target.exists(), str(target)
                links.append(str(target))
    decisions = load(HERE.parent.parent / "decisions.json")
    assert decisions["confirmed_answers"]["Q19"].startswith("B:")
    assert decisions["confirmed_answers"]["Q20"].startswith("A:")
    assert load(HERE / "ai_evaluation_protocol.json")["actual_model_calls"] == 0
    out = {"validated_at_utc": datetime.now(timezone.utc).isoformat(), "checks_passed": True,
           "scope": "Artifact integrity, counts, explicit mappings and local links; not proof of biological validity",
           "json_files_parsed": len(files), "python_files_syntax_checked": len(scripts), "new_source_hashes_matched": len(records), "prior_source_hashes_unchanged": 3,
           "human_donors": 48, "organoid_libraries": 4, "butyrate_samples": 15, "resolved_local_links": len(links),
           "review": {"kind": "Bounded Codex internal mathematical and mapping check", "agent": "/root/hepatic_intervention_data", "reported_max_absolute_recalculation_error": 7.11e-15, "reported_scope": "All 19911 probe means/SDs and four contrasts/SEs; 98 fixed-panel values; no blocking issue", "claude_called": False, "expert_validation": False},
           "execution_issues": ["Platform text first exceeded 12 MB cap and family gzip exceeded 25 MB cap; retried platform-only text under a 100 MB cap and cached successfully", "Additional-input script first used st_size as a function; corrected to the integer attribute and reran successfully using cached source", "One alternative-search subagent did not complete a final report; root re-fetched primary metadata and verified the two relevant files independently"],
           "not_executed": ["Dependency installation pending user response", "Inferential differential-expression analysis", "Single-cell expression linkage and biological validation", "AI model evaluation", "Novel biological discovery"]}
    (HERE / "validation.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out))


if __name__ == "__main__":
    main()
