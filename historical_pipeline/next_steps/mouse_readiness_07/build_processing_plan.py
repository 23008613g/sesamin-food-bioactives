"""Prepare command arguments and a study plan; never install, download reads or run quantification."""
import hashlib
import json
from pathlib import Path, PureWindowsPath

BASE = Path(__file__).resolve().parent


def load(name):
    return json.loads((BASE / name).read_text(encoding="utf-8"))


def sha(name):
    return hashlib.sha256((BASE / name).read_bytes()).hexdigest()


def main():
    fastq = load("fastq_manifest.json")
    mapping = load("orthology/orthology.json")
    headers = load("reference_headers.json")
    tool = load("tool_package_plan.json")
    environment = load("environment_preflight.json")
    eligible = [r for r in mapping["candidates"] if r["primary_comparison_eligible"]]
    assert len(eligible) == 15 and len(mapping["candidates"]) == 25
    assert len(fastq["entries"]) == 18 and fastq["total_bytes"] == 31230644073
    assert all(h["status"] == 200 and h["bytes"] for h in headers)
    assert tool["release_tag"] == "v0.51.1" and tool["archive_inspection_ok"]
    root = PureWindowsPath("D:/ResearchData/sesamin_mouse_2026-10-08")
    executable = str(root / "tools/kallisto-0.51.1/kallisto/kallisto.exe")
    index = str(root / "reference/mouse_ensembl116.idx")
    combined_fasta = str(root / "reference/mouse_ensembl116_combined.fa")
    references = [{**h, "destination": str(root / "reference" / h["url"].rsplit("/", 1)[-1]), "expected_sha256": None, "content_validated": False} for h in headers]
    commands = [{"step": "tool_version_after_authorized_extraction", "argv": [executable, "version"], "executed": False},
                {"step": "build_index_after_reference_validation", "argv": [executable, "index", "-i", index, "-k", "31", combined_fasta], "executed": False}]
    runs = sorted({r["run"] for r in fastq["entries"]})
    for run in runs:
        files = sorted([r for r in fastq["entries"] if r["run"] == run], key=lambda r: r["mate"])
        assert [r["mate"] for r in files] == [1, 2]
        commands.append({"step": "pilot_quantification" if run == runs[0] else "subsequent_quantification", "run": run, "argv": [executable, "quant", "-i", index, "-o", str(root / "quant" / run), "-t", "8", "--plaintext", str(root / "fastq" / files[0]["filename"]), str(root / "fastq" / files[1]["filename"])], "executed": False})
    plan = {
        "user_facing_date": "2026-10-08", "status": "Prepared only; installation and raw-data processing not executed",
        "data_root_proposed": str(root), "data_root_created": root.as_posix() and Path(str(root)).exists(),
        "workspace_policy": "Keep existing OneDrive project reports/code in place. Put only new tool/reference/FASTQ/quant products in the proposed D: project folder; no system PATH or base conda changes.",
        "permission_needed": "User-provided AGENTS instructions prohibit dependency installation without authorization. Current accepted step covers preparation; seek authorization for project-local portable kallisto configuration and the pilot-to-all9 execution before installing.",
        "tool": {"name": "kallisto", "version": "0.51.1", "platform": "official native Windows x64 archive", "archive_sha256": tool["sha256"], "archive_bytes": tool["archive_bytes"], "upstream_digest_available": tool["github_asset_digest"] is not None, "package_read_only_inspected": True, "runtime_tested": False, "reason_for_pinned_version": "Latest checked release0.52.0 lists Linux/macOS assets;0.51.1 has an official Windows executable. This is a native route requiring no WSL setup, not a claim that0.51.1 is newest.", "global_install": False},
        "input_route": "ENA paired FASTQ.gz directly; all18 HEAD200 and exact metadata size. No SRA conversion in primary route; ENA/NGDC/SRA are mirrors of one study.",
        "fastq_bytes": fastq["total_bytes"], "reference_download_bytes": sum(int(h["bytes"]) for h in headers),
        "proposed_free_space_gate_bytes": 100_000_000_000,
        "free_space_gate_interpretation": "Conservative engineering reserve for compressed reads, references/index, quant products and logs; not a measured peak or a promise. Read gzip directly without retaining expanded FASTQ. Recheck live free space before execution.",
        "cpu_threads_per_run": 8, "concurrent_quant_runs": 1,
        "runtime_duration_estimate": None, "pilot_purpose": "Measure actual resource use and processing validity before all9; no wall-clock estimate presented as measured.",
        "pilot_run": runs[0], "pilot_rule": "Lexicographically smallest accession; selected before mouse expression results, not best mapping or desired response.",
        "pilot_fastq_bytes": sum(r["bytes"] for r in fastq["entries"] if r["run"] == runs[0]),
        "reference": {"release": 116, "species": "Mus musculus", "assembly": "GRCm39", "assembly_accession": "GCA_000001635.9", "files": references, "combination_rule": "Read cDNA and ncRNA FASTA headers, retain original transcript versions, verify one transcript maps to exactly one stable gene in same-release GTF, detect duplicates before combining. An identical transcript ID/sequence may be de-duplicated with a log; conflicting duplicates must stop. Preserve all biotypes, not only candidate genes.", "post_download_required": ["gzip integrity", "file SHA256 and exact release/assembly", "FASTA transcript uniqueness", "GTF transcript-to-gene agreement", "all primary mouse gene coverage", "index build log"]},
        "quantification": {"method": "kallisto paired-end transcript abundance estimation", "k": 31, "strandedness": "Unknown in metadata; unstranded primary, no inferred biological pairing", "bootstrap_samples": 0, "bias_option": "Not requested; no obsolete --bias flag assumed for this binary", "note": "Index/help/version and bundled toy test must pass after installation. Planned command compatibility is not an executed runtime test."},
        "outputs_planned": ["transcript_est_counts and TPM for all transcripts in9 libraries", "gene_est_counts=sum transcript estimated counts per gene", "gene_TPM=sum transcript TPM per gene", "effective length and tx2gene tables retained", "per-run QC/version/log/hash records", "all25 human-candidate coverage with15 primary one-to-one mouse pairs"],
        "statistical_boundary": "Estimated counts may be fractional and must not be rounded or called raw integer counts. No DESeq2/edgeR count model, biological pvalues, confidence intervals or successful independent validation claim in this stage.",
        "primary_mouse_gene_ids": [r["mouse_orthologues"][0]["mouse_id"] for r in eligible],
        "primary_human_gene_ids": [r["original_human_id"] for r in eligible],
        "mapping_limits": {"CXCL8": "No mouse orthologue returned in release116; no Cxcl1/Cxcl2 substitution", "ACAT2": "Two orthologues retained as secondary only; do not select one by expression", "no_orthologue_count": 9},
        "exploratory_comparisons_after_QC": {"groups": "ctrl/HFD/sess remain inferred, no matching by suffix", "primary_summary": "Within provisional group arithmetic mean gene TPM; report all9 library values", "contrasts": ["HFD vs ctrl", "sess vs HFD"], "primary_display": "log2((mean_A_TPM+0.1)/(mean_B_TPM+0.1))", "sensitivity_pseudocounts_TPM": [0, 0.5, 1], "zero_denominator_rule": "When pseudocount0 makes either group mean zero, show undefined/NA and raw means rather than inventing a finite log ratio", "low_expression_flag": "Flag gene TPM<1 in all libraries; retain row and do not declare absence or remove from coverage denominator", "cross_species_comparison": "Direction only and complete candidate accounting; human normalized counts and mouse TPM magnitudes not pooled or compared directly", "same_direction_not_causality": True},
        "qc_and_stop_rules": ["Verify full-file MD5 before reading FASTQ; HTTP HEAD is not integrity", "Keep all9 accessions; do not relabel based on expression", "Paired read counts/names and gzip integrity must pass", "Record read lengths, read quality, adapter indicators, processed/pseudoaligned fragments and abundance sums", "Zero processed or pseudoaligned reads, conflicting reference mappings or nonfinite/negative abundance must stop", "Unusually low mapping or inconsistent read structure is a diagnostic flag requiring explanation, not automatic sample deletion", "Tool failure/missing DLL stops; do not install system runtime or modify OS without separate authorization", "Retain input and failure logs; no cleanup/delete commands in this prepared stage"],
        "input_snapshot_sha256": {name: sha(name) for name in ["fastq_manifest.json", "orthology/orthology.json", "design/sample_design.json", "reference_headers.json", "tool_package_plan.json"]},
        "new_external_fees": False, "no_system_python_or_base_conda_modification": True,
    }
    (BASE / "processing_plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (BASE / "planned_commands.json").write_text(json.dumps(commands, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    pilot = [r for r in fastq["entries"] if r["run"] == runs[0]]
    lines = []
    for r in pilot:
        lines.extend([r["url"], "  out=" + r["filename"], "  checksum=md5=" + r["expected_md5"]])
    (BASE / "pilot_downloads.aria2").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = {"runs_planned": len(runs), "primary_gene_pairs": len(eligible), "pilot_fastq_bytes": plan["pilot_fastq_bytes"], "all_fastq_bytes": plan["fastq_bytes"], "reference_bytes": plan["reference_download_bytes"], "data_root_created": plan["data_root_created"], "commands_executed": 0}
    (BASE / "plan_preflight.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
