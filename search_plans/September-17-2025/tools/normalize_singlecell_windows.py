#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SC_DIR = ROOT / 'json' / 'SingleCell'

WIN_BASE = "C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER"
LIB_WIN_BASE = WIN_BASE + "\\SPEC_LIBS\\September-17-2025\\"

def fix_one(fp: Path):
    data = json.loads(fp.read_text(encoding='utf-8'))
    paths = data.get('paths', {})
    ms = paths.get('ms_data', '')
    res = paths.get('results', '')
    lib = paths.get('library', '')

    # Derive amount/time/variant from filename
    # SingleCell__Human_SCP__<Amount>_<Time>__<Variant>.json
    stem = fp.stem
    parts = stem.split('__')
    amount_time = parts[2]
    variant = '__'.join(parts[3:]) if len(parts) > 3 else 'Standard'
    if '_' in amount_time:
        amount, time = amount_time.split('_', 1)
    else:
        amount, time = amount_time, ''

    # Build Windows ms_data
    if amount == 'MBR example':
        ms_data = WIN_BASE + "\\SingleCell\\MBR example\\"
    else:
        ms_data = WIN_BASE + "\\SingleCell\\" + amount + "\\" + time + "\\"

    # Build Windows library
    # Map variant -> library filename
    lib_name = None
    if variant == 'Standard' or variant == 'Standard__MBRoff':
        lib_name = 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025.poin'
    elif variant == 'ARAB':
        lib_name = 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_ARAB.poin'
    elif variant in ('EntrapR1','EntrapR2','EntrapR3'):
        v = variant[-1]
        lib_name = f'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025_entrapR{v}.poin'
    else:
        # Fallback to Standard
        lib_name = 'altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025.poin'
    lib_win = LIB_WIN_BASE + lib_name

    # Build Windows results
    run_id = f"{amount}_{time}_{variant}" if amount != 'MBR example' else f"MBR example_{variant}"
    results = WIN_BASE + "\\SingleCell\\SingleCell_" + run_id

    data['paths']['ms_data'] = ms_data
    data['paths']['library'] = lib_win
    data['paths']['results'] = results
    fp.write_text(json.dumps(data, indent=2), encoding='utf-8')

def main():
    if not SC_DIR.exists():
        return
    for fp in sorted(SC_DIR.glob('*.json')):
        fix_one(fp)

if __name__ == '__main__':
    main()
