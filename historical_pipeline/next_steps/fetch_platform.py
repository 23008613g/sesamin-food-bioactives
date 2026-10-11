"""Download the public platform annotation once; never overwrite source bytes."""
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

BASE = Path(__file__).resolve().parent
URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL31262&targ=self&form=text&view=full"
DEST = BASE / "GPL31262_platform.soft.txt.gz"
if not DEST.exists():
    request = urllib.request.Request(URL, headers={"User-Agent": "Public research data audit"})
    with urllib.request.urlopen(request, timeout=90) as response:
        data = response.read(100_000_001)
    if len(data) > 100_000_000:
        raise ValueError("Platform-only annotation exceeds the 100 MB limit")
    if b"!platform_table_begin" not in data or b"!platform_table_end" not in data:
        raise ValueError("Not the expected GEO annotation")
    with DEST.open("xb") as stream:
        stream.write(gzip.compress(data, mtime=0))
raw = DEST.read_bytes()
record = {"url": URL, "local_file": DEST.name, "bytes": len(raw),
          "sha256": hashlib.sha256(raw).hexdigest(),
          "checked_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "prior_attempts": ["Platform-only text exceeded an initial 12 MB cap; no file saved", "Family SOFT includes sample data and exceeded a 25 MB compressed cap; no file saved. Retried platform-only text with an explicit 100 MB cap."]}
(BASE / "platform_download.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
print(json.dumps(record))
