"""Stage08 execution helpers; no dependency installation or shell evaluation."""
from pathlib import Path
import ctypes, datetime, hashlib, json, os, re, subprocess, sys, time

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / 'mouse_readiness_07'
DATA = Path('D:/ResearchData/sesamin_mouse_2026-10-08')
EXE = DATA / 'tools/kallisto-0.51.1/kallisto/kallisto.exe'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def digest(path, algorithm='sha256'):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, algorithm).hexdigest()

def read_kallisto_info(path):
    """0.51.1 Windows writes unescaped call backslashes; keep raw file untouched."""
    raw=Path(path).read_text(encoding='utf-8')
    match=re.search(r'(?m)^(\s*"call":\s*)"(.*)"(,?\s*)$',raw)
    assert match, 'Expected exact kallisto call field'
    fixed=raw[:match.start()]+match.group(1)+json.dumps(match.group(2))+match.group(3)+raw[match.end():]
    data=json.loads(fixed)
    data['_parsing_note']='Raw Windows kallisto call-field backslashes encoded for JSON parsing only; numerical fields unchanged; original file retained.'
    return data

def run_once(key, argv, input_fingerprints=None, output_paths=None, **kwargs):
    recordpath=HERE/'execution'/(key+'.json')
    if recordpath.exists():
        rec=json.loads(recordpath.read_text(encoding='utf-8'))
        assert rec['status']=='passed' and rec['argv']==[str(x) for x in argv]
        if input_fingerprints is not None:
            assert rec.get('input_fingerprints')==input_fingerprints, 'Inputs changed or provenance absent: '+key
            for path in output_paths or []:
                assert rec.get('output_sha256',{}).get(str(path))==digest(path), 'Saved output changed: '+str(path)
        return rec
    rec=run_logged(key,argv,**kwargs)
    if input_fingerprints is not None:
        rec['input_fingerprints']=input_fingerprints
        rec['output_sha256']={str(p):digest(p) for p in output_paths or []}
        write_json(recordpath,rec)
    return rec

class MemoryCounters(ctypes.Structure):
    _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong)] + [(x, ctypes.c_size_t) for x in ['PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage']]

def run_logged(key, argv, cwd=None, allowed_codes=(0,)):
    logdir = DATA / 'logs'
    logdir.mkdir(parents=True, exist_ok=True)
    recordpath = HERE / 'execution' / (key + '.json')
    if recordpath.exists():
        raise RuntimeError(f'Refusing overwrite of executed command record: {recordpath}')
    rec = {'key':key,'argv':[str(x) for x in argv], 'cwd':str(cwd) if cwd else None,'started_utc':now(),'status':'running','stdout':str(logdir/(key+'.stdout.log')),'stderr':str(logdir/(key+'.stderr.log')),'orchestrator_python':sys.executable,'orchestrator_python_version':sys.version,'driver_script':str(Path(sys.argv[0]).resolve()),'driver_script_sha256':digest(sys.argv[0]) if Path(sys.argv[0]).is_file() else None,'common_helper_sha256':digest(__file__)}
    write_json(recordpath, rec)
    begin = time.monotonic()
    peak = 0
    with open(rec['stdout'],'wb') as out, open(rec['stderr'],'wb') as err:
        p = subprocess.Popen(rec['argv'], cwd=cwd, stdout=out, stderr=err)
        rec['pid'] = p.pid
        write_json(recordpath, rec)
        if os.name == 'nt':
            memfun = ctypes.WinDLL('psapi').GetProcessMemoryInfo
            memfun.argtypes = [ctypes.c_void_p, ctypes.POINTER(MemoryCounters), ctypes.c_ulong]
            memfun.restype = ctypes.c_int
        while True:
            if os.name == 'nt':
                counters = MemoryCounters()
                counters.cb = ctypes.sizeof(counters)
                if memfun(int(p._handle), ctypes.byref(counters), counters.cb):
                    peak = max(peak, counters.PeakWorkingSetSize)
            code = p.poll()
            if code is not None:
                break
            time.sleep(1)
    rec.update(exit_code=code, ended_utc=now(), wall_seconds=time.monotonic()-begin, peak_working_set_bytes_observed=peak, memory_measurement='Windows PeakWorkingSetSize sampled every second for direct child only; not total system/tree peak', status='passed' if code in allowed_codes else 'failed')
    write_json(recordpath,rec)
    print(json.dumps({'command':key,'status':rec['status'],'wall_seconds':round(rec['wall_seconds'],2),'peak_working_set_bytes_observed':peak}), flush=True)
    if code not in allowed_codes:
        raise RuntimeError(f'{key} failed exit={code}; see {rec["stderr"]}')
    return rec
