#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JSON_ROOT = ROOT / 'json' / 'SingleCell'
RUNS_TSV = ROOT / 'runs.tsv'
ALL_JSON = ROOT / 'all_json_concatenated.txt'
CMDS = ROOT / 'commands_by_experiment.txt'
OUT_TSV = ROOT / 'runs_with_msdata.tsv'

LIB_FILENAMES = {
    'Standard': 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025.poin',
    'ARAB': 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_ARAB.poin',
    'EntrapR1': 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_entrapR1.poin',
    'EntrapR2': 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_entrapR2.poin',
    'EntrapR3': 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_entrapR3.poin',
}

def make_params(ms_data_dir: str, lib_path: str, results_dir: str, mbr: bool):
    return {
        "paths": {
            "ms_data": ms_data_dir,
            "library": lib_path,
            "results": results_dir,
        },
        "output": {"write_csv": True, "write_decoys": False, "delete_temp": True, "plots_per_page": 12},
        "logging": {"debug_console_level": 0},
        "global": {
            "scoring": {"q_value_threshold": 0.01},
            "match_between_runs": mbr,
            "ms1_quant": False,
        },
        "first_search": {"fragment_settings": {"min_score": 15}},
        "acquisition": {"nce": 26, "quad_transmission": {"fit_from_data": True}},
        "maxLFQ": {"run_to_run_normalization": False},
    }

def ensure_files(fps):
    for p in fps:
        p.parent.mkdir(parents=True, exist_ok=True)

def scan_amount_times(base: Path):
    amounts = []
    for child in sorted(base.iterdir()):
        if child.is_dir():
            amounts.append(child)
    result = {}
    for amt in amounts:
        times = [d for d in sorted(amt.iterdir()) if d.is_dir() and 'ms' in d.name]
        result[amt.name] = [t.name for t in times]
    return result

def append_runs(rows):
    with RUNS_TSV.open('a', encoding='utf-8') as fh:
        for r in rows:
            fh.write('\t'.join(r) + '\n')

def rebuild_all_json():
    with ALL_JSON.open('w', encoding='utf-8') as out:
        for f in sorted((ROOT / 'json').rglob('*.json')):
            out.write(f'# {f.relative_to(ROOT)}\n')
            out.write(f.read_text(encoding='utf-8'))
            out.write('\n\n')

def append_commands(json_files):
    # Commands will point to the JSONs within this repo to keep it simple
    lines = []
    lines.append('\n######## SingleCell ########\n')
    for f in json_files:
        abs_path = str(f.resolve())
        # Escape backslashes not needed on POSIX; keep as-is
        lines.append(f'SearchDIA("{abs_path}")\n')
    with CMDS.open('a', encoding='utf-8') as fh:
        fh.writelines(lines)

def main():
    ap = argparse.ArgumentParser(description='Generate SingleCell JSONs and update artifacts')
    ap.add_argument('--ms-base', required=True, help='Absolute path to SingleCell data root (contains 10ng, 250pg, 500pg, 1000pg, MBR example)')
    ap.add_argument('--libs-base', required=True, help='Absolute path to SPEC_LIBS/September-17-2025 folder containing SCP libraries')
    ap.add_argument('--results-base', required=False, help='Base folder for results (default: <ms-base>)')
    args = ap.parse_args()

    ms_base = Path(args.ms_base)
    libs_base = Path(args.libs_base)
    results_base = Path(args.results_base) if args.results_base else ms_base

    # Discover amounts and times from filesystem
    amt_times = scan_amount_times(ms_base)
    if not amt_times:
        raise SystemExit(f'No amount folders found under {ms_base}')

    json_files = []
    run_rows = []

    for amount, times in amt_times.items():
        for t in times:
            is_24 = t.lower().startswith('24ms') or t.lower() == '24ms'
            variants = ['Standard'] + (['ARAB', 'EntrapR1', 'EntrapR2', 'EntrapR3'] if is_24 else [])
            mbr_flags = [True]
            if amount.lower().startswith('mbr'):
                mbr_flags = [True, False]
            for mbr in mbr_flags:
                for var in variants:
                    ms_dir = str((ms_base / amount / t).resolve())
                    lib_path = str((libs_base / LIB_FILENAMES[var]).resolve())
                    run_id = f'{amount}_{t}_{("MBRoff" if not mbr else var)}' if amount.lower().startswith('mbr') and var == 'Standard' else f'{amount}_{t}_{var}'
                    results_dir = str((results_base / 'SingleCell' / f'SingleCell_{run_id}').resolve())
                    params = make_params(ms_dir, lib_path, results_dir, mbr)
                    # filename
                    fn = f'SingleCell__Human_SCP__{amount}_{t}__{("Standard__MBRoff" if amount.lower().startswith("mbr") and var=="Standard" and not mbr else var)}.json'
                    jf = JSON_ROOT / fn
                    ensure_files([jf])
                    jf.write_text(json.dumps(params, indent=2), encoding='utf-8')
                    json_files.append(jf)
                    # runs.tsv row
                    row = [
                        'SingleCell', 'SingleCell', 'Human', 'SCP',
                        ('Standard-MBRoff' if amount.lower().startswith('mbr') and var=='Standard' and not mbr else var),
                        'true', '2', 'false', '', 'true' if mbr else 'false',
                        str(jf.relative_to(ROOT)).replace('\\', '/'),
                        results_dir,
                    ]
                    run_rows.append(row)

    # Append to runs.tsv
    append_runs(run_rows)
    # Rebuild concatenated JSONs
    rebuild_all_json()
    # Re-run checker to refresh runs_with_msdata.tsv
    try:
        from check_and_annotate import main as check_main
        check_main()
    except Exception:
        pass
    # Append commands
    append_commands(json_files)

if __name__ == '__main__':
    main()

