import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
samples = [json.loads(line) for line in (RESULTS / 'samples.jsonl').read_text().splitlines()]
environment = json.loads((RESULTS / 'environment.json').read_text())
inputs = json.loads((RESULTS / 'inputs.json').read_text())

groups = defaultdict(list)
for row in samples:
    if row['kind'] == 'training_warmup':
        continue
    groups[(row['kind'], row['scenario'], row.get('mode', ''), row['rows'], row['threads'])].append(row)

summary = []
for (kind, scenario, mode, rows, threads), values in sorted(groups.items()):
    time_key = 'total_seconds' if kind == 'training' else 'seconds'
    seconds = [v[time_key] for v in values]
    record = dict(kind=kind, scenario=scenario, mode=mode, rows=rows, threads=threads,
                  trials=len(values), median_seconds=statistics.median(seconds),
                  min_seconds=min(seconds), max_seconds=max(seconds),
                  median_cpu_seconds=statistics.median(v['cpu_seconds'] for v in values),
                  median_allocated_bytes=statistics.median(v['allocated_bytes'] for v in values),
                  max_prediction_differences=max(v['prediction_differences'] for v in values),
                  max_prediction_abs_difference=max(v['prediction_max_abs'] for v in values),
                  features=values[0]['features'])
    if kind == 'training':
        for field in ('dataset_seconds', 'training_seconds', 'detach_seconds'):
            record['median_' + field] = statistics.median(v[field] for v in values)
        record['row_mode_prediction_differences'] = max(v.get('row_mode_prediction_differences', 0) for v in values)
        record['row_mode_prediction_max_abs'] = max(v.get('row_mode_prediction_max_abs', 0) for v in values)
    summary.append(record)

fields = list(dict.fromkeys(k for row in summary for k in row))
with (ROOT / 'summary.csv').open('w') as output:
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(summary)
(ROOT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

for scenario in sorted({s['scenario'] for s in summary}):
    print('\nSCENARIO', scenario)
    for s in summary:
        if s['scenario'] != scenario:
            continue
        if s['kind'] == 'training' or (s['kind'] == 'prediction_complete' and s['rows'] == 500000):
            print(s['kind'], s['mode'], s['rows'], s['threads'],
                  'median', round(s['median_seconds'], 4),
                  'max prediction differences', s['max_prediction_differences'])
print('\nJob', environment['job'], 'host', environment['host'], 'Julia threads', environment['Julia_default_threads'])
print('Source rows', inputs['source_rows'], 'train pool', inputs['unique_source_train_rows'], 'test pool', inputs['unique_source_test_rows'])
print('Completed', (RESULTS / 'COMPLETED').exists())
