#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # September-17-2025
JSON_ROOT = ROOT / 'json'
RUNS_TSV = ROOT / 'runs.tsv'
OUT_TSV = ROOT / 'runs_with_msdata.tsv'

def load_json(p: Path):
    with p.open('r', encoding='utf-8') as fh:
        return json.load(fh)

def is_windows_abs(p: str) -> bool:
    return len(p) > 2 and p[1] == ':' and (p[2] == '\\' or p[2] == '/')

def is_posix_abs(p: str) -> bool:
    return p.startswith('/')

def find_json_col(cols):
    for i, c in enumerate(cols):
        if c.endswith('.json'):
            return i
    return None

def check_json_file(fp: Path):
    rel = fp.relative_to(ROOT)
    issues = []
    try:
        data = load_json(fp)
    except Exception as e:
        issues.append(f'FAIL parse {rel}: {e}')
        return issues
    paths = data.get('paths', {})
    ms_data = paths.get('ms_data')
    results = paths.get('results')
    if not isinstance(ms_data, str):
        issues.append(f'FAIL {rel}: paths.ms_data missing or not a string')
    else:
        if ms_data.strip() == '.' or ms_data.strip() == '':
            issues.append(f'FAIL {rel}: ms_data is "." or empty')
        elif not (is_windows_abs(ms_data) or is_posix_abs(ms_data)):
            issues.append(f'FAIL {rel}: ms_data is not absolute: {ms_data}')
        # Trailing slash requirement only for Windows style to avoid confusion
        elif is_windows_abs(ms_data) and not ms_data.endswith('\\'):
            issues.append(f'WARN {rel}: ms_data should end with \\: {ms_data}')

    if not isinstance(results, str):
        issues.append(f'FAIL {rel}: paths.results missing or not a string')
    else:
        # Expect results under the same experiment folder as the file's parent
        # e.g., json/<Experiment>/<File>.json -> results startswith C:\...\<Experiment>\<Experiment>_<RunID>
        experiment = fp.parent.name
        exp_marker = f'\\{experiment}\\{experiment}_'
        if exp_marker not in results:
            issues.append(f'WARN {rel}: results does not follow ParentFolder_RunID under {experiment}: {results}')

    # Light library sanity
    lib = paths.get('library', '')
    if not lib.endswith('.poin'):
        issues.append(f'FAIL {rel}: library not .poin: {lib}')

    return issues

def annotate_runs():
    if not RUNS_TSV.exists():
        print(f'No runs.tsv at {RUNS_TSV}', file=sys.stderr)
        return 2

    lines = RUNS_TSV.read_text(encoding='utf-8').splitlines()
    out_lines = []
    errors = 0
    for ln in lines:
        if not ln.strip():
            out_lines.append(ln)
            continue
        cols = ln.split('\t')
        jcol = find_json_col(cols)
        if jcol is None:
            # passthrough
            out_lines.append(ln + '\t')
            continue
        jrel = cols[jcol]
        jpath = (ROOT / jrel).resolve()
        ms_val = ''
        try:
            data = load_json(jpath)
            ms_val = data.get('paths', {}).get('ms_data', '')
        except Exception as e:
            ms_val = f'ERROR:{e}'
            errors += 1
        out_lines.append(ln + '\t' + ms_val)

    OUT_TSV.write_text('\n'.join(out_lines) + '\n', encoding='utf-8')
    return 0 if errors == 0 else 1

def main():
    # Sweep all JSONs for issues
    all_jsons = sorted(JSON_ROOT.rglob('*.json'))
    any_fail = False
    reports = []
    for fp in all_jsons:
        issues = check_json_file(fp)
        reports.extend(issues)
        if any(i.startswith('FAIL') for i in issues):
            any_fail = True

    if reports:
        print('\n'.join(reports))

    rc = annotate_runs()
    # Non-zero if we saw any FAIL or annotate had errors
    if any_fail or rc != 0:
        sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())
