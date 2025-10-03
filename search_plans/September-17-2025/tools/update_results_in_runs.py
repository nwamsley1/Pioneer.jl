#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'runs.tsv'

def find_json_col(cols):
    for i,c in enumerate(cols):
        if c.endswith('.json'):
            return i
    return None

def main():
    lines = RUNS.read_text(encoding='utf-8').splitlines()
    out = []
    header = lines[0]
    out.append(header)
    for ln in lines[1:]:
        if not ln.strip():
            out.append(ln)
            continue
        cols = ln.split('\t')
        jidx = find_json_col(cols)
        if jidx is None:
            out.append(ln)
            continue
        jrel = cols[jidx]
        jpath = (ROOT / jrel).resolve()
        try:
            data = json.loads(jpath.read_text(encoding='utf-8'))
            res = data.get('paths', {}).get('results', cols[-1])
            cols[-1] = res
        except Exception:
            pass
        out.append('\t'.join(cols))
    RUNS.write_text('\n'.join(out)+"\n", encoding='utf-8')
    print('runs.tsv updated results from JSON')

if __name__ == '__main__':
    main()

