"""Retrieve the four separate small MTX archives while the OS archive downloads.

Each worker owns one unique file. The original sequential downloader will
subsequently see completed files and verify their size/hash without refetching.
"""
import concurrent.futures
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

BASE = Path(__file__).resolve().parent
SOURCE = BASE / "sources"


def download(item):
    url = item["url"]
    path = SOURCE / url.rsplit("/", 1)[1]
    if path.exists() or path.with_suffix(path.suffix+".partial").exists():
        raise RuntimeError(f"MTX destination already exists; do not race another writer: {path}")
    temp = path.with_suffix(path.suffix+".parallel_partial")
    if temp.exists():
        raise RuntimeError(f"Existing parallel partial requires inspection: {temp}")
    digest = hashlib.sha256()
    count = 0
    started = time.monotonic()
    with urlopen(url, timeout=90) as response, temp.open("xb") as handle:
        if int(response.headers["Content-Length"]) != item["bytes"]:
            raise ValueError("Changed download size")
        while True:
            block = response.read(1024*1024)
            if not block:
                break
            count += len(block)
            if count > item["bytes"]:
                raise ValueError("Preflight size exceeded")
            handle.write(block)
            digest.update(block)
    if count != item["bytes"]:
        raise ValueError("Incomplete MTX")
    temp.rename(path)
    record = {"url": url, "file": str(path.relative_to(BASE)), "bytes": count,
              "sha256": digest.hexdigest(), "status": "downloaded",
              "elapsed_seconds": round(time.monotonic()-started, 2)}
    print(json.dumps(record), flush=True)
    return record


def main():
    preflight = json.loads((BASE / "preflight.json").read_text())
    if not preflight["download_within_budget"]:
        raise SystemExit("Preflight did not approve budget")
    targets = [x for x in preflight["sources"] if x["url"].endswith("matrix.mtx.gz")]
    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(download, item) for item in targets]
        for future in concurrent.futures.as_completed(futures):
            records.append(future.result())
            output = {"timestamp_utc": datetime.now(timezone.utc).isoformat(),
                      "command": sys.argv, "interpreter": sys.executable,
                      "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "records": records}
            (BASE / "parallel_mtx_download_log.json").write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
