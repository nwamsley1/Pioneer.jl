#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # September-17-2025
JSON_ROOT = ROOT / 'json'

EXPERIMENTS = {
    'MtacThreeProteomeAlternating',
    'MtacThreeProteomeStandard',
    'MtacYeastAlternating3M',
    'MtacYeastAlternating5M',
    'MtacYeastStandard3M',
    'MtacYeastStandard5M',
    'OlsenAstralThreeProteome200ng',
    'OlsenExplorisThreeProteome500ng',
    'SciexThreeProteomeNswath4',
    'EWZ_MBR',
    'YEAST_KO',
}

def run_id_from_filename(name: str) -> str:
    # name like: <Experiment>__<LibClass>__<Variant>.json
    stem = name.rsplit('.', 1)[0]
    parts = stem.split('__')
    variant = parts[-1] if len(parts) >= 2 else stem
    # Normalize: join variant pieces with underscores, and map HuberDelta_XXXX -> HuberXXXX
    rid = variant.replace('__', '_')
    rid = rid.replace('HuberDelta_', 'Huber')
    return rid

def normalize_file(fp: Path):
    exp = fp.parent.name
    try:
        with fp.open('r', encoding='utf-8') as fh:
            data = json.load(fh)
    except Exception as e:
        print(f"SKIP parse error: {fp}: {e}")
        return
    # Normalize ms_data
    data.setdefault('paths', {})
    ms = f"C\\\\Users\\\\n.t.wamsley\\\\Documents\\\\PIONEER_PAPER\\\\{exp}\\\\"
    data['paths']['ms_data'] = ms
    # Normalize results
    rid = run_id_from_filename(fp.name)
    res = f"C\\\\Users\\\\n.t.wamsley\\\\Documents\\\\PIONEER_PAPER\\\\{exp}\\\\{exp}_{rid}"
    data['paths']['results'] = res
    with fp.open('w', encoding='utf-8') as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)

def main():
    targets = []
    for exp in EXPERIMENTS:
        d = JSON_ROOT / exp
        if d.exists():
            targets.extend(sorted(d.glob('*.json')))
    # Also handle Proteomes/YeastHuman for EWZ_MBR
    yh = JSON_ROOT / 'Proteomes' / 'YeastHuman'
    if yh.exists():
        for fp in sorted(yh.glob('*.json')):
            with fp.open('r', encoding='utf-8') as fh:
                data = json.load(fh)
            data.setdefault('paths', {})
            data['paths']['ms_data'] = 'C\\\\Users\\\\n.t.wamsley\\\\Documents\\\\PIONEER_PAPER\\\\EWZ_MBR\\\\'
            # Keep results already set for EWZ_MBR Proteomes files; they were correct
            with fp.open('w', encoding='utf-8') as fh:
                json.dump(data, fh, ensure_ascii=False, indent=2)
    # Normalize targets
    for fp in targets:
        normalize_file(fp)

if __name__ == '__main__':
    main()
