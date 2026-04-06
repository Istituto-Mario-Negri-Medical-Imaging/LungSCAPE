#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import sys
import shutil
import subprocess
from pathlib import Path


# =========================
# CONFIGURATION (EDIT HERE)

DATASET_ID = 201
CONFIGURATION = "3d_fullres"

FOLDS = ["0", "1", "2", "3", "4"]  # Adjust if your model has fewer folds

OUTPUT_SUFFIX = "_airways201.nii.gz"

# =========================

def log(msg: str) -> None:
    print(f"[INFO] {msg}")


def err(msg: str) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)


def which_or_fail(cmd: str) -> str:
    """Ensure a required executable is available in PATH."""
    path = shutil.which(cmd)
    if not path:
        raise RuntimeError(f"Command not found: '{cmd}'")
    return path


def validate_paths(base_dir: Path) -> dict:
    """Build and validate expected project paths."""
    paths = {
        "base": base_dir,
        "downloads": base_dir / "downloads",
        "nnUNet_raw": base_dir / "nnUNet_raw",
        "nnUNet_preprocessed": base_dir / "nnUNet_preprocessed",
        "nnUNet_results": base_dir / "nnUNet_results",
        "INPUT": base_dir / "INPUT",
        "OUTPUT": base_dir / "OUTPUT",
    }

    required_dirs = ["nnUNet_raw", "nnUNet_preprocessed", "nnUNet_results", "OUTPUT"]
    for k in required_dirs:
        if not paths[k].exists():
            raise RuntimeError(f"Required directory not found: {paths[k]}")

    return paths


def check_input_folder(input_dir: Path) -> None:
    """
    Basic sanity check on input filenames.
    nnU-Net typically expects names like case001_0000.nii.gz for single-channel input.
    """
    files = [p for p in input_dir.iterdir() if p.is_file() and p.suffix in (".gz", ".nii")]
    if not files:
        raise RuntimeError(f"No input files found in {input_dir}")

    log(f"Found {len(files)} file(s) in input folder.")
    for p in sorted(files):
        log(f" - {p.name}")


def run_prediction(paths: dict, input_dir: Path) -> None:
    """Run nnUNetv2_predict with configured dataset, config and folds."""
    predict_cmd = which_or_fail("nnUNetv2_predict")

    os.environ["nnUNet_raw"] = str(paths["nnUNet_raw"])
    os.environ["nnUNet_preprocessed"] = str(paths["nnUNet_preprocessed"])
    os.environ["nnUNet_results"] = str(paths["nnUNet_results"])

    cmd = [
        predict_cmd,
        "-i", str(input_dir),
        "-o", str(paths["OUTPUT"]),
        "-d", str(DATASET_ID),
        "-c", CONFIGURATION,
        "-f", *FOLDS,
    ]

    log("Running airways segmentation...")
    log("Command: " + " ".join(cmd))

    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Segmentation failed (exit code {result.returncode}).")

    log("Airways segmentation completed successfully.")


def deploy_outputs(output_dir: Path, data_dir: Path, nnunet_input_dir: Path) -> None:
    """Rename and copy prediction outputs into the DATA directory.

    Reads case_mapping.json from nnUNet_INPUT/ to map caseNNN -> PatientID,
    then saves each output as DATA/PatientID/PatientID_airways201.nii.gz.
    """
    mapping_file = nnunet_input_dir / "case_mapping.json"
    if not mapping_file.exists():
        raise RuntimeError(
            f"case_mapping.json not found in {nnunet_input_dir}.\n"
            "Make sure prepare_dataset.py has been run and --project-dir is correct."
        )
    with open(mapping_file) as f:
        mapping = json.load(f)

    output_files = sorted(p for p in output_dir.iterdir()
                          if p.is_file() and p.name.startswith("case") and p.name.endswith(".nii.gz"))
    if not output_files:
        raise RuntimeError(f"No output files found in {output_dir}")

    deployed, skipped = 0, 0
    for src in output_files:
        case_id = src.name[: -len(".nii.gz")]
        patient_id = mapping.get(case_id)

        if patient_id is None:
            log(f"[WARN] No mapping entry for '{case_id}', skipping")
            skipped += 1
            continue

        patient_dir = data_dir / patient_id
        if not patient_dir.is_dir():
            log(f"[WARN] Patient folder not found: {patient_dir}, skipping")
            skipped += 1
            continue

        dest = patient_dir / f"{patient_id}{OUTPUT_SUFFIX}"
        shutil.copy2(src, dest)
        log(f"Deployed: {src.name} -> {dest.relative_to(data_dir.parent)}")
        deployed += 1

    log(f"Deployed {deployed} file(s) to DATA/.")
    if skipped:
        log(f"[WARN] Skipped {skipped} file(s) (no mapping or missing patient folder).")


def safe_rmtree(path: Path) -> None:
    """Remove a directory tree if it exists."""
    if path.exists():
        shutil.rmtree(path)
        log(f"Removed: {path}")


def cleanup(paths: dict) -> None:
    """
    Remove all nnU-Net folders and downloaded files after prediction.
    This includes nnUNet_results, so the model must be reinstalled before the next run.
    INPUT and OUTPUT are kept.
    """
    safe_rmtree(paths["nnUNet_raw"])
    safe_rmtree(paths["nnUNet_preprocessed"])
    safe_rmtree(paths["nnUNet_results"])
    safe_rmtree(paths["downloads"])

    log("Full cleanup completed.")
    log("The pretrained model was removed (nnUNet_results). Run setup_model201.py again before the next segmentation.")


def main():
    parser = argparse.ArgumentParser(
        description="Run nnU-Net airways segmentation (Dataset201) using INPUT/OUTPUT folders inside the project base directory."
    )
    parser.add_argument(
        "base_dir",
        help="Project base directory created by setup_model201.py",
    )
    parser.add_argument(
        "--project-dir",
        default=None,
        metavar="PROJECT_DIR",
        dest="project_dir",
        help=(
            "Directory containing DATA/ and nnUNet_INPUT/ (created by prepare_dataset.py). "
            "Defaults to the current working directory. "
            "When the expected structure is found, nnUNet_INPUT/ is used directly as the "
            "nnU-Net input folder, and outputs are automatically deployed to "
            "DATA/PatientID/PatientID_airways201.nii.gz using case_mapping.json."
        ),
    )
    parser.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Do not remove nnU-Net folders after successful segmentation (useful for debugging/reuse)",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    project_dir = Path(args.project_dir).expanduser().resolve() if args.project_dir else Path.cwd()

    nnunet_input_dir = project_dir / "nnUNet_INPUT"
    data_dir = project_dir / "DATA"

    explicit_project = args.project_dir is not None
    auto_mode = nnunet_input_dir.is_dir() and data_dir.is_dir()

    if explicit_project and not auto_mode:
        err(
            f"Expected nnUNet_INPUT/ and DATA/ not found in {project_dir}.\n"
            "Make sure this is the output directory from prepare_dataset.py."
        )
        sys.exit(1)

    try:
        paths = validate_paths(base_dir)

        if auto_mode:
            log(f"Project directory: {project_dir}")
            log(f"Using nnUNet_INPUT/ as input: {nnunet_input_dir}")
            input_dir = nnunet_input_dir
        else:
            log("Auto-mode inactive: using INPUT/ folder as-is.")
            input_dir = paths["INPUT"]
            if not input_dir.exists():
                raise RuntimeError(f"Required directory not found: {input_dir}")

        check_input_folder(input_dir)
        run_prediction(paths, input_dir)

        print("\n=== DONE ✅ ===")

        if auto_mode:
            deploy_outputs(paths["OUTPUT"], data_dir, nnunet_input_dir)
            print(f"Airways outputs automatically deployed to {data_dir}")
        else:
            print(f"Airways segmentation output saved in: {paths['OUTPUT']}")
            print("Each output must be renamed to {PatientID}_airways201.nii.gz and placed")
            print("in the corresponding patient folder within the DATA directory.")

        if args.no_cleanup:
            print("Cleanup skipped (--no-cleanup).")
        else:
            cleanup(paths)
            print("Cleanup completed: nnUNet_raw, nnUNet_preprocessed, nnUNet_results, and downloads were removed.")
            print("Run setup_model201.py again before the next segmentation.")

    except Exception as e:
        err(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
