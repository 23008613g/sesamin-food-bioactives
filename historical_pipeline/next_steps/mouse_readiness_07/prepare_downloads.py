"""Build deterministic ENA download manifest; probe headers only. Inspect an official tool archive without installing it."""
import csv
import hashlib
import json
import re
import tarfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
PREVIOUS = BASE.parent / "evidence_closure_06/intervention_access"


def head(row):
    record = {"url": row["url"], "filename": row["filename"], "checked_utc": datetime.now(timezone.utc).isoformat(), "method": "HEAD", "payload_bytes_downloaded": 0}
    try:
        req = urllib.request.Request(row["url"], method="HEAD", headers={"User-Agent": "PublicResearchMetadata/1.0"})
        with urllib.request.urlopen(req, timeout=30) as response:
            record.update({"status": response.status, "content_length": response.headers.get("Content-Length"), "last_modified": response.headers.get("Last-Modified")})
        record["matches_declared_size"] = record["content_length"] is not None and int(record["content_length"]) == row["bytes"]
    except Exception as error:
        record["error"] = str(error)
    return record


def main():
    with (PREVIOUS / "raw_read_manifest.csv").open(encoding="utf-8-sig", newline="") as f:
        old = {r["mirror_run"]: r for r in csv.DictReader(f)}
    with (BASE / "sources/ena_files.tsv").open(encoding="utf-8-sig", newline="") as f:
        ena = list(csv.DictReader(f, delimiter="\t"))
    assert len(ena) == 9 and {r["run_accession"] for r in ena} == set(old)
    entries = []
    for e in sorted(ena, key=lambda r: r["run_accession"]):
        original = old[e["run_accession"]]
        assert e["sample_accession"] == original["biosample_mirror"]
        assert e["study_accession"] == original["mirror_project"] and e["library_layout"] == "PAIRED"
        paths, md5s, sizes = (e[k].split(";") for k in ["fastq_ftp", "fastq_md5", "fastq_bytes"])
        assert len(paths) == len(md5s) == len(sizes) == 2
        for mate, (path, checksum, size) in enumerate(zip(paths, md5s, sizes), 1):
            name = path.rsplit("/", 1)[-1]
            assert name == e["run_accession"] + f"_{mate}.fastq.gz"
            assert re.fullmatch(r"[0-9a-f]{32}", checksum)
            entries.append({"run": e["run_accession"], "sample": e["sample_accession"], "original_label": original["sample_label"], "provisional_group": original["group_label"], "condition_mapping": "inferred from label + study abstract; not independent biology established", "mate": mate, "url": "https://" + path, "filename": name, "bytes": int(size), "expected_md5": checksum})
    headers = list(ThreadPoolExecutor(max_workers=3).map(head, entries))
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "source": "sources/ena_files.tsv", "study": "PRJDB40955 mirrors PRJCA021826", "total_bytes": sum(e["bytes"] for e in entries), "files": len(entries), "runs": 9, "large_payload_downloaded": False, "md5_content_verified": False, "head_checks": "fastq_head_checks.json", "entries": entries}
    (BASE / "fastq_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (BASE / "fastq_head_checks.json").write_text(json.dumps(headers, indent=2) + "\n", encoding="utf-8")
    lines = []
    for e in entries:
        lines.extend([e["url"], "  out=" + e["filename"], "  checksum=md5=" + e["expected_md5"]])
    (BASE / "ena_downloads.aria2").write_text("\n".join(lines) + "\n", encoding="utf-8")

    release = json.loads((BASE / "sources/kallisto_v0.51.1.json").read_text())
    asset = next(a for a in release["assets"] if a["name"] == "kallisto_windows-v0.51.1.tar.gz")
    package = BASE / "sources" / asset["name"]
    info = {"url": asset["browser_download_url"], "declared_size": asset["size"], "github_asset_digest": asset.get("digest"), "release_tag": release["tag_name"], "download_only": True, "extracted": False, "executable_run": False}
    try:
        if not package.exists():
            with urllib.request.urlopen(urllib.request.Request(info["url"], headers={"User-Agent": "PublicResearchMetadata/1.0"}), timeout=40) as response:
                body = response.read(asset["size"] + 1)
            assert len(body) == asset["size"]
            package.write_bytes(body)
        body = package.read_bytes()
        assert len(body) == asset["size"]
        info.update({"archive_bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()})
        if info["github_asset_digest"]:
            assert info["github_asset_digest"] == "sha256:" + info["sha256"]
        info["archive_members"] = []
        with tarfile.open(package, "r:gz") as archive:
            for member in archive.getmembers():
                assert not member.name.startswith(("/", "\\")) and ".." not in Path(member.name).parts
                info["archive_members"].append({"name": member.name, "bytes": member.size, "file": member.isfile()})
        info["archive_inspection_ok"] = True
    except Exception as error:
        info["archive_inspection_ok"] = False
        info["error"] = str(error)
    (BASE / "tool_package_plan.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"runs": 9, "fastq_files": 18, "fastq_bytes": manifest["total_bytes"], "head_200": sum(h.get("status") == 200 for h in headers), "head_size_match": sum(h.get("matches_declared_size", False) for h in headers), "package": info}))


if __name__ == "__main__":
    main()
