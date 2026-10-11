"""Read-only Windows resource/tool check and bounded official metadata fetches."""
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
import winreg
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
SOURCES = BASE / "sources"
SOURCES.mkdir(exist_ok=True)


def powershell(code):
    args = ["powershell.exe", "-NoProfile", "-Command", "[Console]::OutputEncoding=[Text.UTF8Encoding]::new(); " + code]
    result = subprocess.run(args, capture_output=True, timeout=30)
    text = result.stdout.decode("utf-8-sig", errors="replace").strip()
    try:
        payload = json.loads(text)
    except ValueError:
        payload = text
    return {"exit_code": result.returncode, "data": payload, "stderr": result.stderr.decode("utf-8", errors="replace")}


def main():
    now = datetime.now(timezone.utc).isoformat()
    record = {"checked_utc": now, "user_facing_date": "2026-10-08", "read_only": True,
              "python": platform.python_version(), "executable": sys.executable,
              "tool_scope": "PATH plus explicit conventional conda/R roots and current-user WSL registration; no whole-machine scan"}
    record["cpu"] = powershell("Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json")
    record["memory_os"] = powershell("Get-CimInstance Win32_OperatingSystem | Select-Object Caption,Version,TotalVisibleMemorySize,FreePhysicalMemory | ConvertTo-Json")
    record["fixed_disks"] = powershell("Get-CimInstance Win32_LogicalDisk -Filter 'DriveType=3' | Select-Object DeviceID,FileSystem,Size,FreeSpace | ConvertTo-Json")
    record["path_tools"] = {name: shutil.which(name) for name in ["wsl", "docker", "salmon", "kallisto", "STAR", "hisat2", "fastp", "prefetch", "fasterq-dump", "Rscript", "conda", "micromamba", "aria2c"]}
    roots = [Path(os.environ["USERPROFILE"]) / x for x in ["miniconda3", "anaconda3", "miniforge3", "mambaforge"]]
    roots += [Path("C:/ProgramData/miniconda3"), Path("C:/ProgramData/Anaconda3"), Path("<LOCAL_PATH>")]
    record["conventional_roots"] = [{"path": str(p), "exists": p.exists(), "child_directories": [q.name for q in p.iterdir() if q.is_dir()] if p.is_dir() else []} for p in roots]
    record["wsl_commands"] = []
    for command in [["wsl.exe", "--list"], ["wsl.exe", "--status"]]:
        r = subprocess.run(command, capture_output=True, timeout=20)
        raw = r.stdout + r.stderr
        encoding = "utf-16-le" if b"\x00" in raw[:100] else "utf-8"
        out = raw.decode(encoding, errors="replace").strip()
        filename = "wsl_" + command[1].strip("-") + ".txt"
        (SOURCES / filename).write_text(out, encoding="utf-8")
        record["wsl_commands"].append({"argv": command, "exit_code": r.returncode, "decoded_as": encoding, "output_file": "sources/" + filename, "output_start": out[:450]})
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r"Software\Microsoft\Windows\CurrentVersion\Lxss") as key:
            distros = []
            i = 0
            while True:
                try:
                    child = winreg.EnumKey(key, i)
                except OSError:
                    break
                with winreg.OpenKey(key, child) as sub:
                    try:
                        name = winreg.QueryValueEx(sub, "DistributionName")[0]
                    except OSError:
                        name = None
                    distros.append({"registry_key": child, "DistributionName": name})
                i += 1
            record["wsl_current_user_registry"] = {"key_exists": True, "distributions": distros}
    except FileNotFoundError:
        record["wsl_current_user_registry"] = {"key_exists": False, "distributions": []}
    record["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (BASE / "environment_preflight.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    targets = [
        ("kallisto_latest.json", "https://api.github.com/repos/pachterlab/kallisto/releases/latest"),
        ("sra_latest.json", "https://api.github.com/repos/ncbi/sra-tools/releases/latest"),
        ("ena_files.tsv", "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRP018055&result=read_run&fields=run_accession,study_accession,sample_accession,library_layout,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp,submitted_bytes&format=tsv"),
        ("sra_fasterq.md", "https://raw.githubusercontent.com/wiki/ncbi/sra-tools/HowTo:-fasterq-dump.md"),
        ("kallisto_manual.html", "https://pachterlab.github.io/kallisto/manual"),
    ]
    downloads = []
    for filename, url in targets:
        item = {"url": url, "file": "sources/" + filename, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "PublicResearchMetadata/1.0", "Accept": "application/json" if filename.endswith(".json") else "*/*"})
            with urllib.request.urlopen(request, timeout=30) as response:
                body = response.read(4 * 1024 * 1024 + 1)
                if len(body) > 4 * 1024 * 1024:
                    raise ValueError("Metadata response exceeds 4 MB limit")
                item.update({"status": response.status, "content_type": response.headers.get("Content-Type"), "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()})
            (SOURCES / filename).write_bytes(body)
            item["ok"] = True
        except Exception as error:
            item.update({"ok": False, "error": str(error)})
        downloads.append(item)
    (BASE / "source_records.json").write_text(json.dumps(downloads, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    small = {"memory": record["memory_os"]["data"], "disks": record["fixed_disks"]["data"], "wsl_registry": record["wsl_current_user_registry"], "wsl": record["wsl_commands"], "conventional_roots": record["conventional_roots"], "fetch_status": [{"file": d["file"], "ok": d["ok"], "bytes": d.get("bytes"), "error": d.get("error")} for d in downloads]}
    for name in ["kallisto", "sra"]:
        f = SOURCES / (name + "_latest.json")
        if f.exists():
            data = json.loads(f.read_text())
            small[name] = {"tag": data.get("tag_name"), "published_at": data.get("published_at"), "assets": [{"name": a["name"], "size": a["size"], "url": a["browser_download_url"], "digest": a.get("digest")} for a in data.get("assets", [])]}
    print(json.dumps(small, ensure_ascii=False))


if __name__ == "__main__":
    main()
