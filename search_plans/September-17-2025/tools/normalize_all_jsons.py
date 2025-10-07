#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSON_ROOT = ROOT / 'json'

WIN_BASE = "C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER"

def build_run_id(stem: str) -> str:
    # Split after experiment and library class, join the rest
    parts = stem.split('__')
    if len(parts) < 3:
        return stem
    var_parts = parts[2:]
    # Map HuberDelta_XXXX -> HuberXXXX
    var_parts = [p.replace('HuberDelta_', 'Huber') for p in var_parts]
    return '_'.join(var_parts)

def normalize_ms_results(fp: Path, data: dict):
    exp = fp.parent.name
    # Skip SingleCell: already handled separately
    if exp == 'SingleCell':
        return False
    changed = False
    paths = data.setdefault('paths', {})
    # ms_data absolute under experiment
    ms_expected = f"{WIN_BASE}\\{exp}\\"
    # Force to expected value
    if paths.get('ms_data') != ms_expected:
        paths['ms_data'] = ms_expected
        changed = True

    # results ParentFolder_RunID under experiment
    rid = build_run_id(fp.stem)
    res_expected = f"{WIN_BASE}\\{exp}\\{exp}_{rid}"
    if paths.get('results') != res_expected:
        paths['results'] = res_expected
        changed = True

    # Library drive colon/backslashes only (do not change which file)
    lib = paths.get('library')
    if isinstance(lib, str):
        val = lib
        if len(val) >= 2 and val[0].upper() == 'C' and val[1] != ':':
            val = 'C:' + val[1:]
        while '\\' in val and '\\\\' in val:
            val = val.replace('\\\\', '\\')
        if val != lib:
            paths['library'] = val
            changed = True

    data['paths'] = paths
    return changed

def main():
    updated = 0
    for fp in sorted(JSON_ROOT.rglob('*.json')):
        try:
            data = json.loads(fp.read_text(encoding='utf-8'))
        except Exception:
            continue
        if normalize_ms_results(fp, data):
            fp.write_text(json.dumps(data, indent=2), encoding='utf-8')
            updated += 1
    print(f'normalized files: {updated}')

if __name__ == '__main__':
    main()
