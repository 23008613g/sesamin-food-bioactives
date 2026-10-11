"""Small synthetic fixtures only: not study FASTQ validation or benchmarks."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from fastq_qc import run_qc

HERE = Path(__file__).resolve().parent
FIXTURES = HERE / "synthetic_fixtures"
FIXTURES.mkdir(exist_ok=True)
checks = []


def check(name, condition, detail=None):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def record(name, seq=b"ACGTN", quality=None, plus=b"+"):
    return b"@" + name + b"\n" + seq + b"\n" + plus + b"\n" + (quality if quality is not None else b"I" * len(seq)) + b"\n"


def case(name, r1, r2, expected, message=None, mutate=None, sample_pairs=100_000):
    paths = [FIXTURES / f"{name}_R{i}.fastq.gz" for i in (1, 2)]
    for path, content in zip(paths, (r1, r2)):
        path.write_bytes(gzip.compress(content, mtime=0))
    if mutate:
        mutate(paths[0])
    before = [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
    result = run_qc(*paths, sample_pairs=sample_pairs, progress_every=0)
    (FIXTURES / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    check(name + ": expected pass state", result["passed"] is expected, result.get("error"))
    check(name + ": readonly inputs", before == [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths])
    if message:
        check(name + ": diagnostic", message in result.get("error", {}).get("message", ""), result.get("error"))
    if not expected:
        check(name + ": failure not full validation", not result["full_scan_complete"])
    return result, paths


v, valid_paths = case("valid_slash", record(b"read1/1") + record(b"read2/1", b"acgtryswkmbdhvn"),
                     record(b"read1/2") + record(b"read2/2", b"GGC"), True)
check("Valid pair count", v["total_pairs"] == 2)
check("Full base totals R1", v["mates"]["R1"]["total_bases"] == 20)
check("Full base totals R2", v["mates"]["R2"]["total_bases"] == 8)
check("Full lowercase IUPAC count", v["mates"]["R1"]["base_counts_all_records"]["A"] == 2 and v["mates"]["R1"]["base_counts_all_records"]["N"] == 2)
check("Length histogram", v["mates"]["R1"]["length_histogram"] == {5: 1, 15: 1})
check("Both gzip streams verified", v["gzip_crc_and_isize_validation"] == "passed_end_of_both_streams")
check("Quality mean", v["mates"]["R1"]["sample_first_pairs_only"]["mean_phred_assuming_33"] == 40)
case("valid_casava", record(b"machine:1:2:3 1:N:0:ATCG"), record(b"machine:1:2:3 2:N:0:ATCG"), True)
case("valid_ena_style", record(b"DRR123.1 1 length=5", plus=b"+DRR123.1 1 length=5"),
     record(b"DRR123.1 1 length=5", plus=b"+DRR123.1 1 length=5"), True)
case("valid_crlf", record(b"read/1").replace(b"\n", b"\r\n"), record(b"read/2").replace(b"\n", b"\r\n"), True)
case("valid_no_final_newline", record(b"read/1")[:-1], record(b"read/2")[:-1], True)
s, _ = case("sampling_bound", record(b"a") + record(b"b", b"AGATCGGAAGAG") + record(b"c", b"AGATCGGAAGAG", b"!!!!!!!!!!!?"),
            record(b"a") + record(b"b") + record(b"c"), True, sample_pairs=2)
q = s["mates"]["R1"]["sample_first_pairs_only"]
check("First-N sampling boundary", q["reads"] == 2 and q["bases"] == 17 and q["mean_phred_assuming_33"] == 40)
check("Exact adaptor indicator", q["reads_containing_motif"] == 1 and q["fraction_containing_motif"] == 0.5)
check("All records outside sample still counted", s["total_pairs"] == 3 and s["mates"]["R1"]["total_bases"] == 29)
case("mismatched_name", record(b"a"), record(b"b"), False, "read-name mismatch")
case("unequal_counts", record(b"a") + record(b"b"), record(b"a"), False, "unequal mate record counts")
case("truncated_record", b"@a\nACGT\n+\n", record(b"a"), False, "incomplete four-line record")
case("bad_header", record(b"a").replace(b"@a", b">a"), record(b"a"), False, "header lacks @")
case("bad_plus", record(b"a", plus=b"-"), record(b"a"), False, "third line lacks +")
case("wrong_plus_name", record(b"a", plus=b"+b"), record(b"a"), False, "separator name differs")
case("empty_seq", record(b"a", b""), record(b"a"), False, "empty sequence")
case("length_mismatch", record(b"a", quality=b"II"), record(b"a"), False, "length mismatch")
case("invalid_base", record(b"a", b"AC?TN"), record(b"a"), False, "non-DNA-IUPAC")
case("invalid_after_sample", record(b"a") + record(b"b", b"AC?TN"), record(b"a") + record(b"b"), False, "non-DNA-IUPAC", sample_pairs=1)
case("invalid_quality_low", record(b"a", quality=b"II II"), record(b"a"), False, "quality outside")
case("invalid_quality_high", record(b"a", quality=b"II\x7fII"), record(b"a"), False, "quality outside")
case("wrong_slash_marker", record(b"a/2"), record(b"a/2"), False, "wrong /mate suffix")
case("wrong_casava_marker", record(b"a 2:N:0:ATCG"), record(b"a 2:N:0:ATCG"), False, "wrong CASAVA mate")
case("both_empty", b"", b"", False, "zero paired records")


def corrupt_crc(path):
    data = bytearray(path.read_bytes())
    data[-8] ^= 0xFF
    path.write_bytes(data)


case("corrupt_crc", record(b"a"), record(b"a"), False, "CRC check failed", mutate=corrupt_crc)
case("truncated_gzip", record(b"a"), record(b"a"), False, mutate=lambda p: p.write_bytes(p.read_bytes()[:-6]))
same = run_qc(valid_paths[0], valid_paths[0], progress_every=0)
check("Same input path refused", not same["passed"] and "different input paths" in same["error"]["message"])

cli_output = FIXTURES / "cli_valid.json"
cli = subprocess.run([sys.executable, str(HERE / "fastq_qc.py"), *map(str, valid_paths), "--output", str(cli_output), "--progress-every", "0"], capture_output=True, text=True)
check("CLI valid exit", cli.returncode == 0, cli.stderr)
check("CLI valid output", json.loads(cli_output.read_text(encoding="utf-8"))["total_pairs"] == 2)
reject = subprocess.run([sys.executable, str(HERE / "fastq_qc.py"), *map(str, valid_paths), "--output", str(valid_paths[0])], capture_output=True, text=True)
check("CLI refuses source overwrite", reject.returncode == 2 and "must not overwrite an input" in reject.stderr)

report = {
    "type": "synthetic software tests only",
    "study_fastq_checked": False,
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "interpreter": sys.executable,
    "python_version": sys.version,
    "checks": checks,
    "passed": all(x["passed"] for x in checks),
    "total_checks": len(checks),
    "failures": [x for x in checks if not x["passed"]],
    "script_sha256": hashlib.sha256((HERE / "fastq_qc.py").read_bytes()).hexdigest(),
}
(HERE / "synthetic_test_results.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps({"passed": report["passed"], "checks": len(checks), "failures": report["failures"]}))
raise SystemExit(0 if report["passed"] else 1)
