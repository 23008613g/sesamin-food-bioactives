"""Budgeted public-data retrieval for organoid expression linkage (stdlib).

Run preflight first. Downloaded compressed sources are retained without
overwrite. Gzip metadata is not treated as reliable uncompressed size.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

BASE = Path(__file__).resolve().parent
SOURCE = BASE / "sources"
DOWNLOAD_LIMIT = 1_000_000_000
UNCOMPRESSED_LIMIT = 4_000_000_000
OS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207889/suppl/GSE207889_os.h5ad.gz"


def sources():
    manifest = BASE.parent / "organoid_sample_manifest.csv"
    rows = list(csv.DictReader(manifest.open(encoding="utf-8-sig", newline="")))
    return [OS_URL] + [row["matrix_url"] for row in rows]


def head(url):
    result = {"url": url}
    try:
        with urlopen(Request(url, method="HEAD"), timeout=40) as response:
            result.update(status="reachable", http_status=response.status,
                          bytes=int(response.headers["Content-Length"]),
                          last_modified=response.headers.get("Last-Modified"),
                          etag=response.headers.get("ETag"),
                          content_type=response.headers.get("Content-Type"))
    except Exception as exc:
        result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["preflight", "download"])
    args = parser.parse_args()
    BASE.mkdir(parents=True, exist_ok=True)
    if args.mode == "preflight":
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
            records = list(pool.map(head, sources()))
        total = sum(x.get("bytes", 0) for x in records)
        allowed = all(x["status"] == "reachable" for x in records) and total <= DOWNLOAD_LIMIT
        report = {"checked_utc": datetime.now(timezone.utc).isoformat(), "sources": records,
                  "compressed_bytes_total": total, "download_limit_bytes": DOWNLOAD_LIMIT,
                  "uncompressed_limit_bytes": UNCOMPRESSED_LIMIT,
                  "download_within_budget": allowed,
                  "decompression_status": "unknown; enforce byte cap while decompressing only annotation object",
                  "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "command": sys.argv, "interpreter": sys.executable}
        (BASE / "preflight.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
        print(json.dumps(report, ensure_ascii=True))
        if not allowed:
            raise SystemExit(2)
    else:
        preflight = json.loads((BASE / "preflight.json").read_text(encoding="utf-8"))
        if not preflight["download_within_budget"]:
            raise SystemExit("Budget denied by preflight; require explicit revised limit")
        SOURCE.mkdir(exist_ok=True)
        records = []
        for item in preflight["sources"]:
            url = item["url"]
            path = SOURCE / url.rsplit("/", 1)[1]
            digest = hashlib.sha256()
            started = time.monotonic()
            if path.exists():
                with path.open("rb") as handle:
                    for block in iter(lambda: handle.read(1024*1024), b""):
                        digest.update(block)
                received = path.stat().st_size
                status = "existing_file_verified_without_overwrite"
            else:
                temp = path.with_suffix(path.suffix+".partial")
                if temp.exists():
                    raise RuntimeError(f"Partial file exists; inspect before resuming: {temp}")
                received = 0
                with urlopen(Request(url), timeout=90) as response, temp.open("xb") as handle:
                    declared = response.headers.get("Content-Length")
                    if declared and int(declared) != item["bytes"]:
                        raise ValueError("Download Content-Length changed since preflight")
                    while True:
                        block = response.read(1024*1024)
                        if not block:
                            break
                        received += len(block)
                        if received > item["bytes"]:
                            raise ValueError("Per-file preflight byte allowance exceeded")
                        digest.update(block)
                        handle.write(block)
                temp.rename(path)
                status = "downloaded"
            if received != item["bytes"]:
                raise ValueError(f"Incomplete size: {path.name}")
            record = {**item, "local_file": str(path.relative_to(BASE)), "retrieval": status,
                      "sha256": digest.hexdigest(), "received_bytes": received,
                      "elapsed_seconds": round(time.monotonic()-started, 3)}
            records.append(record)
            log = {"retrieved_utc": datetime.now(timezone.utc).isoformat(), "records": records,
                   "command": sys.argv, "interpreter": sys.executable,
                   "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "expression_analysis": False}
            (BASE / "download_log.json").write_text(json.dumps(log, indent=2)+"\n", encoding="utf-8")
            print(json.dumps({"file": path.name, "status": status, "bytes": received}), flush=True)


if __name__ == "__main__":
    main()
