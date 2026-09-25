"""Regenerate with the pinned MannLabs/directlfq checkout on PYTHONPATH.

Reference: c1b5b650b61557c00fafc1f436eeeff6552b25de (Apache-2.0).
Requires numpy, pandas, numba, and multiprocess; not needed to run Julia tests.
"""
import json
import subprocess
import directlfq
from pathlib import Path
import numpy as np
import pandas as pd
from directlfq.protein_intensity_estimation import calculate_peptide_and_protein_intensities

revision = subprocess.check_output([
    "git", "-C", str(Path(directlfq.__file__).resolve().parents[1]), "rev-parse", "HEAD"
], text=True).strip()
assert revision == "c1b5b650b61557c00fafc1f436eeeff6552b25de", revision
rng = np.random.default_rng(312)
cases = []

def add(name, values):
    values = np.asarray(values, dtype=float)
    ids = np.arange(1, len(values) + 1)
    valid_rows = ~np.isnan(values).all(axis=1)
    df = pd.DataFrame(values[valid_rows], index=ids[valid_rows])
    profile, aligned = calculate_peptide_and_protein_intensities(1, df, 10, 1)
    def clean(x):
        return None if not np.isfinite(x) else float(x)
    cases.append(dict(name=name, input=[[clean(x) for x in row] for row in values],
                      expected=[clean(x) for x in profile] if profile is not None else [None] * values.shape[1],
                      selected=[int(x) for x in aligned.index],
                      counts=np.isfinite(aligned.to_numpy()).sum(axis=0).tolist()))

add('uniform', [[1,2,3], [2,3,4]])
add('single_run', [[10], [12], [15]])
add('single_precursor', [[10, 11, np.nan, 12]])
add('single_observation', [[10, np.nan, np.nan]])
add('disconnected', [[10,11,np.nan,np.nan], [12,13,np.nan,np.nan], [np.nan,np.nan,18,20]])
add('ties', np.full((12, 6), 20.0))
add('all_missing_row', [[np.nan,np.nan], [10,12]])
for nprec in (2, 7, 10, 11, 30, 100, 101, 125):
    for nruns in (3, 21):
        values = 20 + rng.normal(size=(nprec, nruns))
        values[rng.random(values.shape) < .35] = np.nan
        add(f'random_{nprec}_{nruns}', values)
# Equal completeness and equal summed log intensity exercise stable precursor ordering.
add('cap_ties', np.tile([19.,20.,21.], (105, 1)))
def format_case(case):
    fields = []
    for key, value in case.items():
        if key == "input":
            rows = ",\n".join("      " + json.dumps(row) for row in value)
            encoded = "[\n" + rows + "\n    ]"
        else:
            encoded = json.dumps(value)
        fields.append("    " + json.dumps(key) + ": " + encoded)
    return "  {\n" + ",\n".join(fields) + "\n  }"

Path(__file__).with_name('reference.json').write_text(
    "[\n" + ",\n".join(map(format_case, cases)) + "\n]\n"
)
