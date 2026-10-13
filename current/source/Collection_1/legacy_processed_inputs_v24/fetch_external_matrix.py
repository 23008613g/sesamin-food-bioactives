"""Restore the sole externally hosted v24 input without changing its bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.request
import zipfile

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Use an already downloaded original v24 ZIP; no network access")
    parser.add_argument("--target-root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    meta = json.loads(Path(__file__).with_name("EXTERNAL_INPUT.json").read_text(encoding="utf-8"))
    target = args.target_root / meta["matrix_member"]
    if target.exists():
        if target.stat().st_size == meta["matrix_bytes"] and sha(target) == meta["matrix_sha256"]:
            print("Verified existing external matrix:", target)
            return
        raise RuntimeError("Destination already exists with different bytes; it was not overwritten")
    with tempfile.TemporaryDirectory(prefix="sesamin_external_") as temporary:
        archive = args.archive
        if archive is None:
            archive = Path(temporary) / "research_reproduction_v24.zip"
            with urllib.request.urlopen(meta["archive_url"], timeout=120) as response, archive.open("wb") as output:
                shutil.copyfileobj(response, output)
        if archive.stat().st_size != meta["archive_bytes"] or sha(archive) != meta["archive_sha256"]:
            raise RuntimeError("Archive size or SHA-256 does not match the pinned research release")
        checked = Path(temporary) / "verified_matrix.gz"
        with zipfile.ZipFile(archive) as z, z.open(meta["matrix_member"]) as source, checked.open("wb") as output:
            shutil.copyfileobj(source, output)
        if checked.stat().st_size != meta["matrix_bytes"] or sha(checked) != meta["matrix_sha256"]:
            raise RuntimeError("Matrix size or SHA-256 mismatch; destination was not written")
        target.parent.mkdir(parents=True, exist_ok=True)
        with checked.open("rb") as source, target.open("xb") as output:
            shutil.copyfileobj(source, output)
    print("Restored and SHA-256 verified external matrix:", target)

if __name__ == "__main__":
    main()
