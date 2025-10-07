#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SC_DIR = ROOT / 'json' / 'SingleCell'

WIN_BASE = r"C:\\Users\\n.t.wamsley\\Documents\\PIONEER_PAPER"
LIB_WIN_BASE = WIN_BASE + r"\\SPEC_LIBS\\September-17-2025\\altimeter_human_len7o40_ch2o3_mc1_SCP_Sep17-2025.poin"

AMOUNTS = ['10ng','250pg','500pg','1000pg']
TIMES = ['3ms','6ms','12ms']  # non-24ms

def make_params(amount: str, time: str):
    ms_data = WIN_BASE + rf"\\SingleCell\\{amount}\\{time}\\"
    results = WIN_BASE + rf"\\SingleCell\\SingleCell_{amount}_{time}_Standard"
    return {
        "paths": {
            "ms_data": ms_data,
            "library": LIB_WIN_BASE,
            "results": results,
        },
        "output": {"write_csv": True, "write_decoys": False, "delete_temp": True, "plots_per_page": 12},
        "logging": {"debug_console_level": 0},
        "global": {"scoring": {"q_value_threshold": 0.01}, "match_between_runs": True, "ms1_quant": False},
        "first_search": {"fragment_settings": {"min_score": 15}},
        "acquisition": {"nce": 26, "quad_transmission": {"fit_from_data": True}},
        "maxLFQ": {"run_to_run_normalization": False},
    }

def main():
    SC_DIR.mkdir(parents=True, exist_ok=True)
    for amount in AMOUNTS:
        for time in TIMES:
            fn = f'SingleCell__Human_SCP__{amount}_{time}__Standard.json'
            fp = SC_DIR / fn
            fp.write_text(json.dumps(make_params(amount, time), indent=2), encoding='utf-8')

if __name__ == '__main__':
    main()

