"""Reassemble revision31 figure assets using the recorded revision32 mapping.
This copies existing assets and verifies hashes; it does not run analysis or fit UMAP.
Usage: python promote_figures.py --previous-root PATH_TO_REVISION31 --output NEW_DIR
"""
import argparse, hashlib, json, shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--previous-root',required=True);p.add_argument('--output',required=True);a=p.parse_args()
src=Path(a.previous_root).resolve();dst=Path(a.output).resolve()
assert not dst.exists(), 'Use a new output directory; originals are never overwritten.'
mapping=json.loads(Path(__file__).with_name('figure_asset_map.json').read_text())
for row in mapping:
    f=src/row['old_path'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'],str(f)
dst.mkdir(parents=True)
for row in mapping:
    target=dst/row['new_path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/row['old_path'],target)
    assert hashlib.sha256(target.read_bytes()).hexdigest()==row['sha256']
print('Verified and copied',len(mapping),'unchanged figure assets')
