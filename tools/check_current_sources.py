#!/usr/bin/env python3
"""Check current source files against their recorded size and SHA-256; no analysis."""

from __future__ import annotations

import hashlib
import json
import ntpath
from pathlib import Path, PurePosixPath
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]


def check_sources(root: Path, records: list[dict]) -> list[str]:
    """Return all manifest/source problems without modifying any files."""
    root = root.resolve()
    errors = []
    seen = set()
    for index, record in enumerate(records):
        try:
            relative = record["path"]
            if not isinstance(relative, str) or not relative:
                raise ValueError("path must be a nonempty string")
            parts = relative.split("/")
            if (
                "\\" in relative
                or PurePosixPath(relative).is_absolute()
                or ntpath.splitdrive(relative)[0]
                or any(part in ("", ".", "..") or ":" in part for part in parts)
            ):
                raise ValueError("path must be a safe repository-relative path")
            key = relative.casefold()
            if key in seen:
                raise ValueError("duplicate manifest path")
            seen.add(key)
            if not isinstance(record["archive_member"], str) or not record["archive_member"]:
                raise ValueError("archive_member is missing")
            expected_size = record["bytes"]
            expected_hash = record["sha256"]
            if type(expected_size) is not int or expected_size < 0:
                raise ValueError("bytes must be a nonnegative integer")
            if not isinstance(expected_hash, str) or len(expected_hash) != 64:
                raise ValueError("sha256 must contain 64 hexadecimal characters")
            int(expected_hash, 16)
            path = root.joinpath(*parts)
            if not path.resolve().is_relative_to(root):
                raise ValueError("path resolves outside the repository")
            if not path.is_file():
                errors.append(f"Missing file: {relative}")
                continue
            if path.stat().st_size != expected_size:
                errors.append(f"Size mismatch: {relative}")
                continue
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest().lower() != expected_hash.lower():
                errors.append(f"SHA-256 mismatch: {relative}")
        except (KeyError, TypeError, ValueError, OSError) as exc:
            errors.append(f"Invalid manifest record {index}: {exc}")
    return errors


def main() -> int:
    manifest = REPO_ROOT / "current" / "source_manifest.json"
    try:
        records = json.loads(manifest.read_text(encoding="utf-8-sig"))
        if not isinstance(records, list) or not records:
            raise ValueError("source_manifest.json must be a nonempty list")
        errors = check_sources(REPO_ROOT, records)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"PASS: all {len(records)} listed files match their recorded byte count and SHA-256.")
    print("This verifies source integrity only; no scientific model was rerun.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
