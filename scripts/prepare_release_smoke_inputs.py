#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare absolute-path config files for pre-release smoke validation."
    )
    parser.add_argument("--repo-root", required=True, help="Path to the repository root")
    parser.add_argument("--output-dir", required=True, help="Directory where smoke inputs should be written")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    predict_config = json.loads((repo_root / "data" / "precompile" / "build_ecoli_altimeter.json").read_text(encoding="utf-8"))
    predict_config["library_path"] = str(output_dir / "predicted_library")
    predict_config["out_dir"] = str(output_dir)
    predict_config["lib_name"] = "predicted_library"
    predict_config["new_lib_name"] = "predicted_library"
    predict_config["out_name"] = "predicted_library.tsv"
    predict_config["max_koina_requests"] = 4
    predict_config["max_koina_batch"] = 250
    predict_config["fasta_paths"] = [str((repo_root / "data" / "precompile" / "ecoli.fasta").resolve())]
    calibration_raw_file = str((repo_root / "data" / "ecoli_test" / "raw" / "ecoli_filtered_01.arrow").resolve())
    predict_config["calibration_raw_file"] = calibration_raw_file
    predict_config["library_params"]["calibration_raw_file"] = calibration_raw_file
    write_json(output_dir / "build_predict.json", predict_config)

    search_config = json.loads((repo_root / "test" / "integration" / "search_ecoli.json").read_text(encoding="utf-8"))
    search_config["paths"]["library"] = str((repo_root / "test" / "integration" / "ecoli_lib.poin").resolve())
    search_config["paths"]["ms_data"] = str((repo_root / "data" / "ecoli_test" / "raw").resolve())
    search_config["paths"]["results"] = str((output_dir / "search_results").resolve())
    search_config["output"]["write_csv"] = False
    search_config["output"]["write_decoys"] = False
    write_json(output_dir / "search.json", search_config)


if __name__ == "__main__":
    main()
