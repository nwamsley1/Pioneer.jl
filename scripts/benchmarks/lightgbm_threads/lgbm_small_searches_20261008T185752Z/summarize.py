import json
import pathlib
import statistics

root = pathlib.Path(__file__).resolve().parent
samples = [json.loads(line) for line in (root / "results/samples.jsonl").read_text().splitlines()]
measured = [row for row in samples if row["kind"] == "pass1" and row["repeat"] > 0]
metric_keys = ["n_total", "n_fold0", "n_fold1", "n_pool_fold0", "n_pool_fold1",
               "selected_iteration", "pool_targets_q01", "pool_decoys_q01", "features"]
references = {row["case"]: row for row in samples
              if row["kind"] == "pass1" and row["repeat"] == 0 and row["threads"] == 24}
for row in measured:
    reference = references[row["case"]]
    assert row["sidecars_identical"] == (row["sidecar_values"] == reference["sidecar_values"])
    assert row["metrics_identical"] == all(row[key] == reference[key] for key in metric_keys)
    assert all(call["threads"] == row["threads"] for call in row["fit_calls"] + row["predict_calls"])
if (root / "results/COMPLETED").exists():
    assert len(measured) == 40 and len(references) == 4
summary = []
for case in sorted({row["case"] for row in measured}):
    by_threads = {threads: sorted(
        [row for row in measured if row["case"] == case and row["threads"] == threads],
        key=lambda row: row["repeat"],
    ) for threads in (24, 48)}
    pairs = {row["repeat"]: row for row in by_threads[24]}
    changes = [100 * (row["seconds"] / pairs[row["repeat"]]["seconds"] - 1)
               for row in by_threads[48] if row["repeat"] in pairs]
    entry = {"case": case, "paired_change_percent": changes,
             "median_paired_change_percent": statistics.median(changes) if changes else None,
             "variants": []}
    for threads, rows in by_threads.items():
        if not rows:
            continue
        entry["variants"].append({
            "threads": threads, "trials": len(rows), "files": rows[0]["files"],
            "rows": rows[0]["n_total"], "features": len(rows[0]["features"]),
            "median_seconds": statistics.median(row["seconds"] for row in rows),
            "min_seconds": min(row["seconds"] for row in rows),
            "max_seconds": max(row["seconds"] for row in rows),
            "median_fit_seconds": statistics.median(row["fit_seconds"] for row in rows),
            "median_prediction_seconds": statistics.median(row["prediction_seconds"] for row in rows),
            "median_cpu_seconds": statistics.median(row["cpu_seconds"] for row in rows),
            "median_gc_seconds": statistics.median(row["gc_seconds"] for row in rows),
            "median_allocated_bytes": statistics.median(row["allocated_bytes"] for row in rows),
            "selected_iterations": sorted({row["selected_iteration"] for row in rows}),
            "target_counts_q01": sorted({row["pool_targets_q01"] for row in rows}),
            "all_metrics_identical": all(row["metrics_identical"] for row in rows),
            "all_sidecars_identical": all(row["sidecars_identical"] for row in rows),
            "fit_row_counts": sorted({call["rows"] for row in rows for call in row["fit_calls"]}),
        })
    summary.append(entry)
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
