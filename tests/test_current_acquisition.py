"""Small offline safety tests; no real archive or network access required."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import stat
import tempfile
import unittest
from unittest.mock import patch
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


acquisition = load_tool("acquire_current")
checker = load_tool("check_current_sources")


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="sesamin-acquisition-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.archive = self.root / "fixture.zip"
        self.destination = self.root / "output"

    def make_zip(self, members):
        with zipfile.ZipFile(self.archive, "w", zipfile.ZIP_DEFLATED) as bundle:
            for name, data in members:
                if isinstance(name, str):
                    raw_name = name
                    name = zipfile.ZipInfo(raw_name)
                    # Preserve intentionally hostile bytes in fixtures even
                    # on Windows, where ZipInfo normally normalizes slashes.
                    name.filename = raw_name
                    name.orig_filename = raw_name
                bundle.writestr(name, data)
        payload = self.archive.read_bytes()
        return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}

    def test_verified_archive_extracts_bytes_and_does_not_execute_code(self):
        metadata = self.make_zip([
            ("Collection_1/results/value.csv", b"gene,value\nA,1\n"),
            ("restore_exact_duplicates.py", b"raise RuntimeError('must not run')\n"),
        ])
        self.assertEqual(acquisition.acquire(self.archive, self.destination, metadata), 2)
        self.assertEqual((self.destination / "Collection_1/results/value.csv").read_bytes(), b"gene,value\nA,1\n")
        self.assertTrue((self.destination / "restore_exact_duplicates.py").is_file())

    def test_hash_mismatch_refuses_before_creating_destination(self):
        metadata = self.make_zip([("value.csv", b"a,b\n")])
        metadata["sha256"] = "0" * 64
        with self.assertRaisesRegex(acquisition.AcquisitionError, "SHA-256 mismatch"):
            acquisition.acquire(self.archive, self.destination, metadata)
        self.assertFalse(self.destination.exists())

    def test_size_mismatch_refuses_before_creating_destination(self):
        metadata = self.make_zip([("value.csv", b"a,b\n")])
        metadata["bytes"] += 1
        with self.assertRaisesRegex(acquisition.AcquisitionError, "size mismatch"):
            acquisition.acquire(self.archive, self.destination, metadata)
        self.assertFalse(self.destination.exists())

    def test_existing_destination_and_its_files_are_untouched(self):
        metadata = self.make_zip([("existing.txt", b"new")])
        self.destination.mkdir()
        retained = self.destination / "existing.txt"
        retained.write_bytes(b"original")
        with self.assertRaisesRegex(acquisition.AcquisitionError, "already exists"):
            acquisition.acquire(self.archive, self.destination, metadata)
        self.assertEqual(retained.read_bytes(), b"original")

    def test_unsafe_paths_refused_before_extracting_any_member(self):
        names = ["../escape.txt", "nested/../../escape.txt", "/absolute.txt", "C:/escape.txt",
                 "C:escape.txt", "folder\\escape.txt", "folder/./escape.txt", "file.txt:stream", "NUL.txt"]
        for name in names:
            with self.subTest(name=name):
                metadata = self.make_zip([("safe.txt", b"safe"), (name, b"unsafe")])
                with self.assertRaises(acquisition.AcquisitionError):
                    acquisition.acquire(self.archive, self.destination, metadata)
                self.assertFalse(self.destination.exists())

    def test_symlink_member_refused(self):
        link = zipfile.ZipInfo("link")
        link.create_system = 3
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        metadata = self.make_zip([(link, b"../target")])
        with self.assertRaisesRegex(acquisition.AcquisitionError, "Links or special"):
            acquisition.acquire(self.archive, self.destination, metadata)
        self.assertFalse(self.destination.exists())

    def test_case_collisions_and_file_parent_collision_refused(self):
        cases = [[("FILE.txt", b"one"), ("file.txt", b"two")],
                 [("parent", b"one"), ("parent/child.txt", b"two")]]
        for entries in cases:
            with self.subTest(entries=entries):
                metadata = self.make_zip(entries)
                with self.assertRaises(acquisition.AcquisitionError):
                    acquisition.acquire(self.archive, self.destination, metadata)
                self.assertFalse(self.destination.exists())

    def test_no_arguments_never_download_or_extract(self):
        metadata = {"filename": "fixture.zip", "bytes": 10, "sha256": "a" * 64,
                    "release_tag": "fixture", "download_url": "https://example.invalid/fixture.zip"}
        with patch.object(acquisition, "read_metadata", return_value=metadata), \
             patch.object(acquisition, "download_archive") as download, \
             patch.object(acquisition, "acquire") as acquire, \
             patch("builtins.print"):
            self.assertEqual(acquisition.main([]), 0)
            download.assert_not_called()
            acquire.assert_not_called()

    def test_source_checker_detects_changed_missing_and_unsafe_paths(self):
        payload = b"print('source')\n"
        source = self.root / "analysis.py"
        source.write_bytes(payload)
        record = {"path": "analysis.py", "archive_member": "Collection_1/analysis.py",
                  "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
        self.assertEqual(checker.check_sources(self.root, [record]), [])
        source.write_bytes(b"X" * len(payload))
        self.assertIn("SHA-256 mismatch", checker.check_sources(self.root, [record])[0])
        source.unlink()
        self.assertIn("Missing file", checker.check_sources(self.root, [record])[0])
        unsafe = dict(record, path="../outside.py")
        self.assertIn("safe repository-relative", checker.check_sources(self.root, [unsafe])[0])


if __name__ == "__main__":
    unittest.main()
