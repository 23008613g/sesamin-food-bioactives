"""Read-only, dependency-free full paired FASTQ.gz structural QC.

Structural checks inspect every record; quality/adaptor summaries inspect the
deterministic first N pairs only and do not imply random sampling or trimming.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import re
import sys
import time
import zlib

IUPAC = b"ACGTRYSWKMBDHVN"
BAD_BASE = re.compile(rb"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]")
BAD_QUALITY = re.compile(rb"[^\x21-\x7e]")
CASAVA = re.compile(rb"^([12]):[YN]:[0-9]+:")
ADAPTOR_MOTIF = b"AGATCGGAAGAG"
BASE_BUFFER_SIZE = 4 * 1024 * 1024


class QCError(Exception):
    pass


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def without_newline(line):
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def name_and_marker(header, expected_mate, pair_number):
    if not header.startswith(b"@"):
        raise QCError(f"Pair {pair_number}, R{expected_mate}: header lacks @")
    fields = header[1:].split(None, 1)
    if not fields or not fields[0]:
        raise QCError(f"Pair {pair_number}, R{expected_mate}: empty read name")
    token = fields[0]
    normalized = token
    style = []
    if token.endswith((b"/1", b"/2")):
        if token[-1] != ord(str(expected_mate)):
            raise QCError(f"Pair {pair_number}, R{expected_mate}: wrong /mate suffix")
        normalized = token[:-2]
        style.append("slash_mate_suffix")
    if len(fields) == 2:
        match = CASAVA.match(fields[1])
        if match:
            if int(match.group(1)) != expected_mate:
                raise QCError(f"Pair {pair_number}, R{expected_mate}: wrong CASAVA mate field")
            style.append("casava_mate_field")
    if not normalized:
        raise QCError(f"Pair {pair_number}, R{expected_mate}: empty normalized read name")
    return normalized, token, "+".join(style) or "no_explicit_mate_marker"


class MateStats:
    def __init__(self):
        self.records = 0
        self.total_bases = 0
        self.length_histogram = Counter()
        self.base_counts = Counter()
        self._bases = bytearray()
        self.header_styles = Counter()
        self.sample_reads = 0
        self.sample_bases = 0
        self.sample_quality_ascii = Counter()
        self.sample_adaptor_reads = 0
        self.sample_adaptor_start_histogram = Counter()

    def flush_bases(self):
        if self._bases:
            seq = self._bases.upper()
            for base in IUPAC:
                self.base_counts[chr(base)] += seq.count(base)
            self._bases.clear()

    def add(self, sequence, quality, style, sampled):
        self.records += 1
        self.total_bases += len(sequence)
        self.length_histogram[len(sequence)] += 1
        self._bases.extend(sequence)
        if len(self._bases) >= BASE_BUFFER_SIZE:
            self.flush_bases()
        self.header_styles[style] += 1
        if sampled:
            self.sample_reads += 1
            self.sample_bases += len(sequence)
            self.sample_quality_ascii.update(quality)
            where = sequence.upper().find(ADAPTOR_MOTIF)
            if where >= 0:
                self.sample_adaptor_reads += 1
                self.sample_adaptor_start_histogram[where + 1] += 1

    def summary(self):
        self.flush_bases()
        acgt = sum(self.base_counts[x] for x in "ACGT")
        q = self.sample_quality_ascii
        quality_n = sum(q.values())
        return {
            "records_validated": self.records,
            "total_bases": self.total_bases,
            "length_min": min(self.length_histogram) if self.length_histogram else None,
            "length_max": max(self.length_histogram) if self.length_histogram else None,
            "length_histogram": dict(sorted(self.length_histogram.items())),
            "base_counts_all_records": dict(sorted(self.base_counts.items())),
            "gc_fraction_among_acgt": (self.base_counts["G"] + self.base_counts["C"]) / acgt if acgt else None,
            "n_fraction_all_bases": self.base_counts["N"] / self.total_bases if self.total_bases else None,
            "other_iupac_fraction_all_bases": sum(self.base_counts[x] for x in "RYSWKMBDHV") / self.total_bases if self.total_bases else None,
            "header_style_counts": dict(self.header_styles),
            "sample_first_pairs_only": {
                "reads": self.sample_reads,
                "bases": self.sample_bases,
                "quality_ascii_histogram": dict(sorted(q.items())),
                "quality_ascii_min": min(q) if q else None,
                "quality_ascii_max": max(q) if q else None,
                "phred_offset_assumption": 33,
                "quality_encoding_empirically_proven": False,
                "mean_phred_assuming_33": sum((k - 33) * v for k, v in q.items()) / quality_n if quality_n else None,
                "q20_base_fraction_assuming_33": sum(v for k, v in q.items() if k >= 53) / quality_n if quality_n else None,
                "q30_base_fraction_assuming_33": sum(v for k, v in q.items() if k >= 63) / quality_n if quality_n else None,
                "exact_adaptor_motif": ADAPTOR_MOTIF.decode("ascii"),
                "reads_containing_motif": self.sample_adaptor_reads,
                "fraction_containing_motif": self.sample_adaptor_reads / self.sample_reads if self.sample_reads else None,
                "first_motif_start_position_1_based_histogram": dict(sorted(self.sample_adaptor_start_histogram.items())),
            },
        }


def check_record(stream, expected_mate, pair_number):
    raw_header = stream.readline()
    if not raw_header:
        return None
    raw_seq, raw_plus, raw_quality = stream.readline(), stream.readline(), stream.readline()
    if not raw_seq or not raw_plus or not raw_quality:
        raise QCError(f"Pair {pair_number}, R{expected_mate}: incomplete four-line record")
    header, seq, plus, quality = map(without_newline, (raw_header, raw_seq, raw_plus, raw_quality))
    normalized, original_token, style = name_and_marker(header, expected_mate, pair_number)
    if not plus.startswith(b"+"):
        raise QCError(f"Pair {pair_number}, R{expected_mate}: third line lacks +")
    if plus[1:].strip():
        plus_name = plus[1:].split(None, 1)[0]
        if plus_name.startswith(b"@"):
            plus_name = plus_name[1:]
        if plus_name != original_token:
            raise QCError(f"Pair {pair_number}, R{expected_mate}: repeated separator name differs from header")
    if not seq:
        raise QCError(f"Pair {pair_number}, R{expected_mate}: empty sequence")
    if len(seq) != len(quality):
        raise QCError(f"Pair {pair_number}, R{expected_mate}: sequence/quality length mismatch {len(seq)} != {len(quality)}")
    if BAD_BASE.search(seq):
        raise QCError(f"Pair {pair_number}, R{expected_mate}: sequence has non-DNA-IUPAC character")
    if BAD_QUALITY.search(quality):
        raise QCError(f"Pair {pair_number}, R{expected_mate}: quality outside printable ASCII 33..126")
    return normalized, seq, quality, style


def run_qc(r1, r2, sample_pairs=100_000, progress_every=1_000_000):
    start = time.perf_counter()
    stats = [MateStats(), MateStats()]
    result = {
        "schema_version": 1,
        "started_utc": utc_now(),
        "python_version": sys.version,
        "input_files": [],
        "status": "failed",
        "passed": False,
        "full_scan_complete": False,
        "gzip_crc_and_isize_validation": "not_established",
        "total_pairs": 0,
        "sampling": {
            "method": "deterministic first N paired records, no random sampling",
            "requested_pairs": sample_pairs,
            "scope": "Quality-score and exact-adaptor-motif summaries only; all structural checks, lengths and base counts scan every validated pair",
        },
        "limitations": [
            "Technical file QC does not verify sample identity, treatment labels or biological independence.",
            "No duplicate read-name scan across different pairs and no sequence-based contamination screen.",
            "Exact 12-base adaptor motif is a limited indicator; no mismatch, other-adaptor or biological-origin discrimination.",
            "First-N summaries are not random or representative-library estimates; no automatic trimming or exclusion occurs.",
            "Quality scores assume Phred+33 for summaries only; ASCII validity does not establish encoding.",
            "Input MD5/SHA256 must be checked separately against the download manifest.",
            "On failure, counts describe validated pairs before the first failure, not the full library.",
        ],
    }
    try:
        for path in (r1, r2):
            meta = path.stat()
            result["input_files"].append({"path": str(path.resolve()), "compressed_bytes": meta.st_size, "mtime_ns": meta.st_mtime_ns})
        if r1.resolve() == r2.resolve():
            raise QCError("R1 and R2 must be different input paths")
        with gzip.open(r1, "rb") as f1, gzip.open(r2, "rb") as f2:
            while True:
                pair_no = result["total_pairs"] + 1
                a = check_record(f1, 1, pair_no)
                b = check_record(f2, 2, pair_no)
                if a is None and b is None:
                    break
                if a is None or b is None:
                    raise QCError(f"Pair {pair_no}: unequal mate record counts (R1 EOF={a is None}, R2 EOF={b is None})")
                if a[0] != b[0]:
                    raise QCError(f"Pair {pair_no}: mate read-name mismatch {a[0][:200]!r} != {b[0][:200]!r}")
                sampled = pair_no <= sample_pairs
                stats[0].add(a[1], a[2], a[3], sampled)
                stats[1].add(b[1], b[2], b[3], sampled)
                result["total_pairs"] = pair_no
                if progress_every and pair_no % progress_every == 0:
                    print(json.dumps({"event": "progress", "pairs_validated": pair_no, "elapsed_seconds": round(time.perf_counter() - start, 3)}), flush=True)
        if result["total_pairs"] == 0:
            raise QCError("Both files contain zero paired records")
        for path, before in zip((r1, r2), result["input_files"]):
            after = path.stat()
            if after.st_size != before["compressed_bytes"] or after.st_mtime_ns != before["mtime_ns"]:
                raise QCError(f"Input changed while reading: {path}")
        result.update(status="passed", passed=True, full_scan_complete=True,
                      gzip_crc_and_isize_validation="passed_end_of_both_streams")
    except (QCError, OSError, EOFError, ValueError, zlib.error) as exc:
        result["error"] = {"type": type(exc).__name__, "message": str(exc)}
    result["mates"] = {"R1": stats[0].summary(), "R2": stats[1].summary()}
    result["completed_utc"] = utc_now()
    result["elapsed_seconds"] = round(time.perf_counter() - start, 6)
    result["pairs_per_second"] = result["total_pairs"] / result["elapsed_seconds"] if result["elapsed_seconds"] else None
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("r1", type=Path)
    parser.add_argument("r2", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-pairs", type=int, default=100_000)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()
    if args.sample_pairs < 0 or args.progress_every < 0:
        parser.error("sample-pairs and progress-every must be nonnegative")
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    if any(path.resolve() in (args.output.resolve(), temporary.resolve()) for path in (args.r1, args.r2)):
        parser.error("Output must not overwrite an input")
    result = run_qc(args.r1, args.r2, args.sample_pairs, args.progress_every)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"event": "complete", "passed": result["passed"], "pairs": result["total_pairs"], "elapsed_seconds": result["elapsed_seconds"], "output": str(args.output)}), flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
