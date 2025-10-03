#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SC_DIR = ROOT / 'json' / 'SingleCell'
RUNS = ROOT / 'runs.tsv'
ALL = ROOT / 'all_json_concatenated.txt'
CMDS = ROOT / 'commands_by_experiment.txt'

def parse_name(stem: str):
    # SingleCell__Human_SCP__<Amount>_<Time>__<Variant>
    parts = stem.split('__')
    amount, time = parts[2].split('_', 1) if '_' in parts[2] else (parts[2], '')
    variant = parts[3]
    return amount, time, variant

def build_row(fp: Path, data: dict):
    amount, time, variant = parse_name(fp.stem)
    mbr = 'false' if variant == 'Standard__MBRoff' else 'true'
    lib_variant = 'Standard-MBRoff' if variant == 'Standard__MBRoff' else variant
    res = data['paths']['results']
    return [
        'SingleCell', 'SingleCell', 'Human', 'SCP', lib_variant,
        'true', '2', 'false', '', mbr,
        str(fp.relative_to(ROOT)).replace('\\','/'), res
    ]

def write_runs(rows):
    lines = RUNS.read_text(encoding='utf-8').splitlines()
    header = lines[0]
    body = [ln for ln in lines[1:] if not ln.startswith('SingleCell\t')]
    body.extend('\t'.join(r) for r in rows)
    RUNS.write_text('\n'.join([header]+body)+"\n", encoding='utf-8')

def rebuild_all():
    with ALL.open('w', encoding='utf-8') as out:
        for f in sorted((ROOT/'json').rglob('*.json')):
            out.write(f'# {f.relative_to(ROOT)}\n')
            out.write(f.read_text(encoding='utf-8'))
            out.write('\n\n')

def sync_commands(files):
    lines = CMDS.read_text(encoding='utf-8').splitlines()
    # Remove previous SingleCell block if present
    new = []
    in_sc = False
    for ln in lines:
        if ln.startswith('######## SingleCell ########'):
            in_sc = True
            continue
        if in_sc and ln.startswith('######## '):
            in_sc = False
            new.append(ln)
            continue
        if not in_sc:
            new.append(ln)
    # Append fresh block
    new.append('')
    new.append('######## SingleCell ########')
    for f in sorted(files):
        new.append(f'SearchDIA("{f.resolve()}")')
    new.append('')
    CMDS.write_text('\n'.join(new)+"\n", encoding='utf-8')

def main():
    files = sorted(SC_DIR.glob('*.json'))
    rows = []
    for fp in files:
        data = json.loads(fp.read_text(encoding='utf-8'))
        rows.append(build_row(fp, data))
    write_runs(rows)
    rebuild_all()
    sync_commands(files)

if __name__ == '__main__':
    main()

