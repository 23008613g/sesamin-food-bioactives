"""Resume a stopped contiguous OS gzip download using validated HTTP ranges.

Four workers write disjoint byte intervals in the existing partial file.
No second compressed copy is created. Never run beside fetch_inputs.py.
"""
import concurrent.futures
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

BASE = Path(__file__).resolve().parent
SOURCE = BASE / "sources"
NAME = "GSE207889_os.h5ad.gz"


def digest_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    preflight = json.loads((BASE / "preflight.json").read_text())
    item = next(x for x in preflight["sources"] if x["url"].endswith(NAME))
    partial = SOURCE / (NAME+".partial")
    final = SOURCE / NAME
    plan_file = BASE / "range_resume_log.json"
    if final.exists() or plan_file.exists():
        raise RuntimeError("Final file or range plan already exists; inspect state instead of overwriting")
    prefix = partial.stat().st_size
    total = item["bytes"]
    if not 0 < prefix < total:
        raise ValueError("Expected incomplete contiguous prefix")
    plan = {"started_utc": datetime.now(timezone.utc).isoformat(), "url": item["url"],
            "command": sys.argv, "interpreter": sys.executable,
            "script_sha256": digest_file(Path(__file__)),
            "retained_prefix_bytes": prefix, "retained_prefix_sha256": digest_file(partial),
            "expected_total_bytes": total,
            "reason": "Original single HTTP response stalled; its verified workspace-specific processes were stopped before resuming",
            "completed_ranges": []}
    plan_file.write_text(json.dumps(plan, indent=2)+"\n", encoding="utf-8")
    # Preallocation is within the already approved compressed size allowance.
    with partial.open("r+b") as handle:
        handle.truncate(total)
    chunk = math.ceil((total-prefix)/4)
    intervals = [(start, min(start+chunk-1, total-1)) for start in range(prefix, total, chunk)]

    def worker(interval):
        start, end = interval
        req = Request(item["url"], headers={"Range": f"bytes={start}-{end}"})
        read = 0
        digest = hashlib.sha256()
        started = time.monotonic()
        with urlopen(req, timeout=90) as response, partial.open("r+b") as handle:
            expected_range = f"bytes {start}-{end}/{total}"
            if response.status != 206 or response.headers.get("Content-Range") != expected_range:
                raise ValueError("Server did not honor exact byte interval")
            handle.seek(start)
            while True:
                block = response.read(1024*1024)
                if not block:
                    break
                read += len(block)
                if read > end-start+1:
                    raise ValueError("Range exceeded declared length")
                handle.write(block)
                digest.update(block)
        if read != end-start+1:
            raise ValueError("Range was incomplete")
        record = {"start": start, "end": end, "bytes": read,
                  "sha256": digest.hexdigest(), "elapsed_seconds": round(time.monotonic()-started, 3)}
        print(json.dumps(record), flush=True)
        return record

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(worker, interval) for interval in intervals]
        for future in concurrent.futures.as_completed(futures):
            plan["completed_ranges"].append(future.result())
            plan_file.write_text(json.dumps(plan, indent=2)+"\n", encoding="utf-8")
    if sum(x["bytes"] for x in plan["completed_ranges"])+prefix != total:
        raise ValueError("Range coverage does not match complete file")
    partial.rename(final)
    plan.update(status="completed", final_sha256=digest_file(final),
                finished_utc=datetime.now(timezone.utc).isoformat())
    plan_file.write_text(json.dumps(plan, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "completed", "bytes": total, "sha256": plan["final_sha256"]}), flush=True)


if __name__ == "__main__":
    main()
