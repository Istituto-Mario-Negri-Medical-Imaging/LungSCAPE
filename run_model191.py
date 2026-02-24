#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path


# =========================
# CONFIGURATION (EDIT HERE)

DATASET_ID = 191
CONFIGURATION = "3d_fullres"

FOLDS = ["0", "1", "2", "3", "4"]  # Adjust if your model has fewer folds

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

    required_dirs = ["nnUNet_raw", "nnUNet_preprocessed", "nnUNet_results", "INPUT", "OUTPUT"]
    for k in required_dirs:
        if not paths[k].exists():
            raise RuntimeError(f"Required directory not found: {paths[k]}")

    return paths


def check_input_folder(input_dir: Path) -> None:
    """
    Basic sanity check on input filenames.
    nnU-Net typically expects names like case001_0000.nii.gz for single-channel input.
    """
    files = [p for p in input_dir.iterdir() if p.is_file()]
    if not files:
        raise RuntimeError(f"No input files found in {input_dir}")

    log(f"Found {len(files)} file(s) in INPUT folder.")
    for p in sorted(files):
        log(f" - {p.name}")


def run_prediction(paths: dict) -> None:
    """Run nnUNetv2_predict with configured dataset, config and folds."""
    predict_cmd = which_or_fail("nnUNetv2_predict")

    env = os.environ.copy()
    env["nnUNet_raw"] = str(paths["nnUNet_raw"])
    env["nnUNet_preprocessed"] = str(paths["nnUNet_preprocessed"])
    env["nnUNet_results"] = str(paths["nnUNet_results"])

    cmd = [
        predict_cmd,
        "-i", str(paths["INPUT"]),
        "-o", str(paths["OUTPUT"]),
        "-d", str(DATASET_ID),
        "-c", CONFIGURATION,
        "-f", *FOLDS,
    ]

    log("Running segmentation...")
    log("Command: " + " ".join(cmd))

    result = subprocess.run(cmd, env=env, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Segmentation failed (exit code {result.returncode}).")

    log("Segmentation completed successfully.")


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
    log("The pretrained model was removed (nnUNet_results). Run setup_model.py again before the next segmentation.")


def main():
    parser = argparse.ArgumentParser(
        description="Run nnU-Net segmentation using INPUT/OUTPUT folders inside the project base directory."
    )
    parser.add_argument(
        "base_dir",
        help="Project base directory created by setup_model.py",
    )
    parser.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Do not remove nnU-Net folders after successful segmentation (useful for debugging/reuse)",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()

    try:
        paths = validate_paths(base_dir)
        check_input_folder(paths["INPUT"])
        run_prediction(paths)

        print("\n=== DONE ✅ ===")
        print(f"Segmentation output saved in: {paths['OUTPUT']. It must be renamed [idpat]_HighAttenuation.nii.gz and placed in the patient folder for next processing.}")

        if args.no_cleanup:
            print("Cleanup skipped (--no-cleanup).")
        else:
            cleanup(paths)
            print("Cleanup completed: nnUNet_raw, nnUNet_preprocessed, nnUNet_results, and downloads were removed.")
            print("Run setup_model.py again before the next segmentation.")

    except Exception as e:
        err(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
