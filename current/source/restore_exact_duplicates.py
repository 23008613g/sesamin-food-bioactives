"""Restore exact duplicate paths locally; does not download or run analyses."""
from pathlib import Path, PurePosixPath
import argparse,hashlib,json,shutil

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def safe(root,name):
    posix=PurePosixPath(name)
    if posix.is_absolute() or '..' in posix.parts or '\\' in name:
        raise ValueError('Unsafe path in restoration map: '+name)
    path=(root/Path(*posix.parts)).resolve()
    if not path.is_relative_to(root):raise ValueError('Path outside bundle: '+name)
    return path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--restore',action='store_true',help='Copy verified canonical bytes to missing alias paths')
    args=ap.parse_args();root=args.root.resolve()
    entries=json.loads((root/'EXACT_DUPLICATES.json').read_text(encoding='utf-8'))['aliases']
    checked={};missing=[]
    # Complete preflight before writing any aliases.
    for entry in entries:
        src=safe(root,entry['canonical_path']);dst=safe(root,entry['alias_path'])
        if str(src) not in checked:checked[str(src)]=sha(src)
        if checked[str(src)]!=entry['sha256']:raise ValueError('Canonical content hash mismatch: '+str(src))
        if dst.exists():
            if not dst.is_file() or sha(dst)!=entry['sha256']:raise ValueError('Refusing to overwrite different content: '+str(dst))
        else:missing.append((entry,src,dst))
    restored=0
    if args.restore:
        for entry,src,dst in missing:
            dst.parent.mkdir(parents=True,exist_ok=True)
            # Exclusive creation protects files that appeared after preflight.
            with src.open('rb') as inp,dst.open('xb') as out:shutil.copyfileobj(inp,out,1024*1024)
            if sha(dst)!=entry['sha256']:raise ValueError('Restored content hash mismatch: '+str(dst))
            restored+=1
    print(json.dumps({'alias_count':len(entries),'missing_before':len(missing),'restored':restored,
        'mode':'restore' if args.restore else 'check','network_used':False,'analysis_executed':False}))

if __name__=='__main__':main()
