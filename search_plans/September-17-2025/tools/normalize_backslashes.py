#!/usr/bin/env python3
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSON_ROOT = ROOT / 'json'

WIN_BASE = 'C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER'

def norm_slashes(s: str, ensure_trailing=False):
    if not isinstance(s, str):
        return s
    # Normalize forward to backslashes
    s = s.replace('/', '\\')
    # Collapse multiple backslashes to a single backslash
    s = re.sub(r'\\+', r'\\', s)
    # Ensure C: drive colon if path starts with C and then backslash
    if len(s) >= 2 and s[0].upper() == 'C' and s[1] != ':':
        s = 'C:' + s[1:]
    # Optionally ensure trailing backslash
    if ensure_trailing and not s.endswith('\\'):
        s += '\\'
    return s

def fix_file(fp: Path):
    try:
        data = json.loads(fp.read_text(encoding='utf-8'))
    except Exception:
        return False
    changed = False
    p = data.get('paths', {})
    for key, trailing in (('ms_data', True), ('results', False), ('library', False)):
        v = p.get(key)
        if isinstance(v, str):
            nv = norm_slashes(v, ensure_trailing=trailing if key=='ms_data' else False)
            if nv != v:
                p[key] = nv
                changed = True
    data['paths'] = p
    if changed:
        fp.write_text(json.dumps(data, indent=2), encoding='utf-8')
    return changed

def main():
    total = 0
    for fp in sorted(JSON_ROOT.rglob('*.json')):
        if fix_file(fp):
            total += 1
    print('normalized backslashes in', total, 'files')

if __name__ == '__main__':
    main()

