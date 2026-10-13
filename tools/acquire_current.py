#!/usr/bin/env python3
"""Verify and unpack the current research supplement; no implicit network access."""

from __future__ import annotations

import argparse
import hashlib
import json
import ntpath
import os
from pathlib import Path, PurePosixPath
import stat
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "current" / "archive.json"
DEFAULT_DEST = REPO_ROOT / ".local" / "current_data"


class AcquisitionError(ValueError):
    """An archive or destination did not satisfy the acquisition contract."""


def read_metadata(path: Path) -> dict:
    metadata = json.loads(path.read_text(encoding="utf-8-sig"))
    for field in ("filename", "sha256", "bytes", "release_tag", "download_url"):
        if field not in metadata:
            raise AcquisitionError(f"Archive metadata is missing {field!r}.")
    digest = metadata["sha256"]
    if not isinstance(digest, str) or len(digest) != 64:
        raise AcquisitionError("The expected SHA-256 must contain 64 hexadecimal characters.")
    try:
        int(digest, 16)
    except ValueError as exc:
        raise AcquisitionError("The expected SHA-256 is not hexadecimal.") from exc
    if type(metadata["bytes"]) is not int or metadata["bytes"] <= 0:
        raise AcquisitionError("The expected archive size must be a positive integer.")
    return metadata


def verify_archive(archive: Path, expected_bytes: int, expected_sha256: str) -> None:
    if not archive.is_file():
        raise AcquisitionError(f"Archive is not a regular file: {archive}")
    size = archive.stat().st_size
    if size != expected_bytes:
        raise AcquisitionError(f"Archive size mismatch: expected {expected_bytes}, found {size} bytes.")
    digest = hashlib.sha256()
    with archive.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest().lower() != expected_sha256.lower():
        raise AcquisitionError("Archive SHA-256 mismatch; nothing was extracted.")


def _member_parts(info: zipfile.ZipInfo) -> tuple[str, ...]:
    # ZipInfo normalizes backslashes on Windows and truncates NULs. Inspect the
    # original central-directory name so those transformations cannot hide an
    # unsafe archive path.
    name = info.orig_filename
    if not name or "\\" in name or "\x00" in name:
        raise AcquisitionError(f"Unsafe archive member: {name!r}")
    path = PurePosixPath(name)
    raw_parts = name.rstrip("/").split("/")
    if path.is_absolute() or ntpath.splitdrive(name)[0]:
        raise AcquisitionError(f"Absolute or drive-qualified archive member: {name!r}")
    reserved = {"CON", "PRN", "AUX", "NUL"} | {
        f"{prefix}{number}" for prefix in ("COM", "LPT") for number in range(1, 10)
    }
    for part in raw_parts:
        if (
            part in ("", ".", "..")
            or ":" in part
            or part.endswith((".", " "))
            or any(ord(char) < 32 or char in '<>"|?*' for char in part)
            or part.split(".")[0].upper() in reserved
        ):
            raise AcquisitionError(f"Unsafe archive member: {name!r}")
    mode = (info.external_attr >> 16) & 0xFFFF
    kind = stat.S_IFMT(mode)
    if stat.S_ISLNK(mode) or kind not in (0, stat.S_IFREG, stat.S_IFDIR):
        raise AcquisitionError(f"Links or special files are not allowed: {name!r}")
    if kind == stat.S_IFDIR and not info.is_dir():
        raise AcquisitionError(f"Inconsistent directory entry: {name!r}")
    if info.flag_bits & 1:
        raise AcquisitionError(f"Encrypted members are not supported: {name!r}")
    return tuple(raw_parts)


def safe_extract(archive: Path, destination: Path) -> int:
    """Extract into a new directory after checking every member; never overwrite.

    The caller must verify the expected whole-archive hash before using this
    function. If a read/CRC error occurs after extraction starts, the new partial
    directory is left for inspection and will not be reused automatically.
    """
    destination = Path(os.path.abspath(destination))
    if os.path.lexists(destination):
        raise AcquisitionError(f"Destination already exists; choose a NEW --dest: {destination}")
    with zipfile.ZipFile(archive) as bundle:
        members = []
        known: dict[tuple[str, ...], bool] = {}
        for info in bundle.infolist():
            parts = _member_parts(info)
            key = tuple(part.casefold() for part in parts)
            if key in known:
                raise AcquisitionError(f"Duplicate or case-colliding member: {info.filename!r}")
            known[key] = info.is_dir()
            members.append((info, parts))
        for key in known:
            for length in range(1, len(key)):
                if key[:length] in known and not known[key[:length]]:
                    raise AcquisitionError("An archive file is also used as a parent directory.")
        destination.mkdir(parents=True, exist_ok=False)
        count = 0
        for info, parts in members:
            output = destination.joinpath(*parts)
            if info.is_dir():
                output.mkdir(parents=True, exist_ok=True)
                continue
            output.parent.mkdir(parents=True, exist_ok=True)
            with bundle.open(info) as source, output.open("xb") as target:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    target.write(block)
            count += 1
    return count


def acquire(archive: Path, destination: Path, metadata: dict) -> int:
    if os.path.lexists(destination):
        raise AcquisitionError(f"Destination already exists; choose a NEW --dest: {destination}")
    verify_archive(archive, metadata["bytes"], metadata["sha256"])
    return safe_extract(archive, destination)


def download_archive(metadata: dict, output: Path) -> None:
    url = metadata["download_url"]
    if not isinstance(url, str) or not url.startswith("https://"):
        raise AcquisitionError("The pinned download URL must use HTTPS.")
    request = urllib.request.Request(url, headers={"User-Agent": "sesamin-archive-acquisition/1"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, output.open("wb") as target:
            size = 0
            for block in iter(lambda: response.read(1024 * 1024), b""):
                size += len(block)
                if size > metadata["bytes"]:
                    raise AcquisitionError("Download exceeds the expected archive size.")
                target.write(block)
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 403, 404):
            raise AcquisitionError(
                f"GitHub returned HTTP {exc.code}. The release may be private or unavailable. "
                "Sign in to GitHub in your browser, download the pinned release asset, "
                "then use --archive PATH. This tool does not read login credentials."
            ) from exc
        raise AcquisitionError(f"Download failed: HTTP {exc.code}.") from exc
    except urllib.error.URLError as exc:
        raise AcquisitionError(f"Download failed: {exc.reason}. You can use --archive PATH instead.") from exc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--archive", type=Path, help="Use an already downloaded local archive.")
    source.add_argument("--download", action="store_true", help="Explicitly download the pinned archive.")
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST, help="NEW extraction directory.")
    args = parser.parse_args(argv)
    temporary: Path | None = None
    try:
        metadata = read_metadata(DEFAULT_MANIFEST)
        if not args.archive and not args.download:
            print(f"Archive: {metadata['filename']} ({metadata['bytes']} bytes)")
            print(f"SHA-256: {metadata['sha256']}")
            print(f"Release: {metadata['release_tag']}")
            print(f"Pinned URL: {metadata['download_url']}")
            status = "exists (not checked or modified)" if os.path.lexists(args.dest) else "does not exist"
            print(f"Destination: {args.dest} [{status}]")
            print("No download or extraction performed. Use --archive PATH or --download explicitly.")
            print("A private release requires a browser download after GitHub sign-in, then --archive PATH.")
            return 0
        if os.path.lexists(args.dest):
            raise AcquisitionError(f"Destination already exists; choose a NEW --dest: {args.dest}")
        archive = args.archive
        if args.download:
            handle, filename = tempfile.mkstemp(prefix="sesamin-current-", suffix=".zip")
            os.close(handle)
            temporary = Path(filename)
            print(f"Downloading the pinned release: {metadata['download_url']}")
            download_archive(metadata, temporary)
            archive = temporary
        count = acquire(archive, args.dest, metadata)
        print(f"Verified archive size and SHA-256; extracted {count} files to {args.dest}.")
        print("No alias script, dependency installation, or scientific analysis was executed.")
        return 0
    except (AcquisitionError, OSError, ValueError, zipfile.BadZipFile, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
