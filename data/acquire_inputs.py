"""Acquire checksummed v24 inputs separately; source terms remain applicable.

Default performs only a presence/hash check. --archive imports inputs from a
local copy of the existing versioned archive. --download explicitly downloads
that public archive first. No archive code is executed, and no file is uploaded.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=Path(__file__).with_name('acquisition_manifest.json')

def sha(path):
    result=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):result.update(block)
    return result.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--archive',type=Path)
    group.add_argument('--download',action='store_true')
    args=parser.parse_args()
    manifest=json.loads(MANIFEST.read_text(encoding='utf-8'))
    target=ROOT/'reproduction'
    archive=args.archive
    with tempfile.TemporaryDirectory(prefix='sesamin-inputs-') as temp:
        if args.download:
            archive=Path(temp)/'research_reproduction_v24.zip'
            print('Downloading the existing public v24 research artifact; consult original data terms.')
            with urllib.request.urlopen(manifest['archive_url']) as response,archive.open('wb') as out:
                shutil.copyfileobj(response,out)
        if archive:
            if sha(archive)!=manifest['archive_sha256']:raise SystemExit('Archive checksum mismatch; no input written.')
            with zipfile.ZipFile(archive) as z:
                for row in manifest['inputs']:
                    destination=(target/row['file']).resolve()
                    destination.relative_to(target.resolve())
                    payload=z.read(manifest['archive_prefix']+row['file'])
                    if hashlib.sha256(payload).hexdigest()!=row['sha256']:
                        raise SystemExit('Input checksum mismatch: '+row['file'])
                    if destination.exists():
                        if sha(destination)!=row['sha256']:raise SystemExit('Refusing to overwrite differing input: '+row['file'])
                    else:
                        destination.parent.mkdir(parents=True,exist_ok=True)
                        destination.write_bytes(payload)
        missing=[]
        for row in manifest['inputs']:
            path=target/row['file']
            if not path.exists():missing.append(row['file'])
            elif sha(path)!=row['sha256']:raise SystemExit('Existing input checksum mismatch: '+row['file'])
        print(json.dumps({'expected_inputs':len(manifest['inputs']),'missing':missing,'source_inputs_tracked_in_git':False},indent=2))
        if missing:raise SystemExit(2)

if __name__=='__main__':main()
