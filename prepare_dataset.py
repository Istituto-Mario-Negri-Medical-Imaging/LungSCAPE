#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LungScape - Step 0: Dataset Preparation

Converts DICOM datasets to NIfTI format, resamples to standard resolution,
and organizes the working directory structure expected by the LungScape pipeline.

Usage:
    python prepare_dataset.py <dicom_dir> <output_dir> [options]

Output structure:
    <output_dir>/
    ├── DATA/
    │   ├── Patient001/
    │   │   └── Patient001.nii.gz
    │   ├── Patient002/
    │   │   └── Patient002.nii.gz
    │   └── ...
    └── nnUNet_INPUT/          (if --no-nnunet-input is not set)
        ├── case001_0000.nii.gz
        ├── case002_0000.nii.gz
        ├── ...
        └── case_mapping.json

Requirements:
    - dcm2niix (system command, must be in PATH)
    - antspyx (pip install antspyx)
    - nibabel (pip install nibabel)
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import ants
import nibabel as nib


# -- Constants ----------------------------------------------------------------

TARGET_MATRIX = 512         # Target in-plane matrix size (512 x 512)
TARGET_Z_SPACING = 0.5      # Target slice spacing in mm


# -- Conversion helpers -------------------------------------------------------

def ants_to_nibabel(ants_img):
    """Convert an ANTsImage to a nibabel Nifti1Image via a temp file."""
    with TemporaryDirectory() as tmpdir:
        tmpfile = os.path.join(tmpdir, "tmp.nii")
        ants_img.to_filename(tmpfile)
        nib_img = nib.load(tmpfile)
        # Force data into memory before the temp file is deleted
        _ = nib_img.get_fdata()
    return nib_img


def convert_dicom_to_nifti(dicom_dir, output_dir):
    """Run dcm2niix on a DICOM directory.

    Returns the path to the generated .nii file, or None on failure.
    dcm2niix is called with ``-f %f`` so the output filename matches the
    DICOM folder name.
    """
    os.makedirs(output_dir, exist_ok=True)

    result = subprocess.run(
        ["dcm2niix", "-z", "n", "-f", "%f", "-o", output_dir, str(dicom_dir)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"  [ERROR] dcm2niix failed:\n{result.stderr}", file=sys.stderr)
        return None

    # Find the generated .nii file(s)
    nii_files = sorted(Path(output_dir).glob("*.nii"))
    if not nii_files:
        print(f"  [ERROR] No .nii files produced by dcm2niix", file=sys.stderr)
        return None

    if len(nii_files) > 1:
        print(f"  [WARN] Multiple NIfTI files produced, using first: {nii_files[0].name}")

    return str(nii_files[0])


def resample_volume(nii_path):
    """Resample a NIfTI volume to 512x512 matrix and 0.5 mm z-spacing.

    The resampling follows two stages (matching the original pipeline logic):
      1. If the in-plane matrix is not 512x512, resample to 512x512 voxels
         (use_voxels=True, nearest-neighbor interpolation).
      2. Resample z-spacing to TARGET_Z_SPACING mm (use_voxels=False, nearest-
         neighbor interpolation).

    Returns a nibabel Nifti1Image with the resampled data.
    """
    img_nib = nib.load(nii_path)
    img_ants = ants.image_read(nii_path)
    shape = img_nib.shape

    # Stage 1: resample in-plane matrix to 512x512 if needed
    if shape[0] != TARGET_MATRIX or shape[1] != TARGET_MATRIX:
        print(f"  Resampling matrix {shape[0]}x{shape[1]} -> {TARGET_MATRIX}x{TARGET_MATRIX}")
        new_shape = (TARGET_MATRIX, TARGET_MATRIX, shape[2])
        img_ants = ants.resample_image(img_ants, new_shape, use_voxels=True, interp_type=0)

        # Reload header to get updated spacing
        img_nib = ants_to_nibabel(img_ants)

    # Stage 2: resample z-spacing to TARGET_Z_SPACING mm
    sx, sy, sz = img_nib.header.get_zooms()[:3]
    if abs(sz - TARGET_Z_SPACING) > 0.01:
        print(f"  Resampling z-spacing {sz:.3f} mm -> {TARGET_Z_SPACING} mm")
        new_spacing = (float(sx), float(sy), TARGET_Z_SPACING)
        img_ants = ants.resample_image(img_ants, new_spacing, use_voxels=False, interp_type=0)

    return ants_to_nibabel(img_ants)


# -- Main pipeline ------------------------------------------------------------

def process_patient(dicom_folder, data_dir, patient_id, skip_resampling=False):
    """Process a single patient: DICOM -> NIfTI -> resample -> save.

    Creates ``data_dir/patient_id/patient_id.nii.gz``.
    Returns the path to the saved NIfTI, or None on failure.
    """
    print(f"\n  Converting DICOM -> NIfTI ...")

    with TemporaryDirectory() as tmpdir:
        nii_path = convert_dicom_to_nifti(dicom_folder, tmpdir)
        if nii_path is None:
            return None

        # Resample
        if skip_resampling:
            print(f"  Skipping resampling (--skip-resampling)")
            result_img = nib.load(nii_path)
        else:
            result_img = resample_volume(nii_path)

    # Save to DATA structure
    patient_dir = os.path.join(data_dir, patient_id)
    os.makedirs(patient_dir, exist_ok=True)
    output_path = os.path.join(patient_dir, f"{patient_id}.nii.gz")
    nib.save(result_img, output_path)
    print(f"  Saved: {output_path}")

    return output_path


def prepare_nnunet_input(data_dir, nnunet_input_dir):
    """Create nnU-Net INPUT folder with caseNNN_0000.nii.gz naming.

    Also generates a case_mapping.json file mapping case IDs to patient IDs.
    """
    os.makedirs(nnunet_input_dir, exist_ok=True)
    mapping = {}

    patient_dirs = sorted(
        d for d in Path(data_dir).iterdir()
        if d.is_dir()
    )

    for idx, patient_dir in enumerate(patient_dirs, start=1):
        patient_id = patient_dir.name
        ct_path = patient_dir / f"{patient_id}.nii.gz"

        if not ct_path.exists():
            print(f"  [WARN] CT not found for {patient_id}, skipping nnU-Net input")
            continue

        case_id = f"case{idx:03d}"
        dest = os.path.join(nnunet_input_dir, f"{case_id}_0000.nii.gz")
        shutil.copy2(str(ct_path), dest)
        mapping[case_id] = patient_id

    # Save mapping
    mapping_file = os.path.join(nnunet_input_dir, "case_mapping.json")
    with open(mapping_file, "w") as f:
        json.dump(mapping, f, indent=2)

    print(f"\nnnU-Net INPUT prepared: {len(mapping)} cases")
    print(f"  Mapping saved to: {mapping_file}")


def main():
    parser = argparse.ArgumentParser(
        description="LungScape Step 0: Convert DICOM dataset to NIfTI and prepare working directories."
    )
    parser.add_argument(
        "dicom_dir",
        help="Directory containing patient DICOM subfolders",
    )
    parser.add_argument(
        "output_dir",
        help="Output base directory (DATA/ and nnUNet_INPUT/ will be created here)",
    )
    parser.add_argument(
        "--skip-resampling",
        action="store_true",
        help="Skip resampling (keep original matrix and spacing)",
    )
    parser.add_argument(
        "--no-nnunet-input",
        action="store_true",
        help="Do not prepare nnU-Net INPUT folder",
    )
    args = parser.parse_args()

    dicom_dir = Path(args.dicom_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    data_dir = output_dir / "DATA"
    nnunet_input_dir = output_dir / "nnUNet_INPUT"

    # Validate input
    if not dicom_dir.is_dir():
        print(f"[ERROR] DICOM directory not found: {dicom_dir}", file=sys.stderr)
        sys.exit(1)

    # Check dcm2niix is available
    if shutil.which("dcm2niix") is None:
        print("[ERROR] dcm2niix not found in PATH. Install it first.", file=sys.stderr)
        print("  Ubuntu/Debian: sudo apt install dcm2niix", file=sys.stderr)
        print("  conda: conda install -c conda-forge dcm2niix", file=sys.stderr)
        sys.exit(1)

    # Find patient folders
    patient_folders = sorted(
        d for d in dicom_dir.iterdir()
        if d.is_dir()
    )

    if not patient_folders:
        print(f"[ERROR] No patient subfolders found in {dicom_dir}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("  LungScape - Step 0: Dataset Preparation")
    print("=" * 60)
    print(f"\n  DICOM source:  {dicom_dir}")
    print(f"  Output:        {output_dir}")
    print(f"  Patients:      {len(patient_folders)}")
    print(f"  Resampling:    {'OFF' if args.skip_resampling else f'{TARGET_MATRIX}x{TARGET_MATRIX}, z={TARGET_Z_SPACING}mm'}")
    print(f"  nnU-Net INPUT: {'OFF' if args.no_nnunet_input else 'ON'}")

    # Process each patient
    os.makedirs(data_dir, exist_ok=True)
    successful = 0
    failed = []

    for idx, patient_folder in enumerate(patient_folders, start=1):
        patient_id = patient_folder.name
        print(f"\n[{idx}/{len(patient_folders)}] {patient_id}")
        print("-" * 40)

        try:
            result = process_patient(
                patient_folder, str(data_dir), patient_id,
                skip_resampling=args.skip_resampling,
            )
            if result is not None:
                successful += 1
            else:
                failed.append(patient_id)
        except Exception as e:
            print(f"  [ERROR] {e}", file=sys.stderr)
            failed.append(patient_id)
            continue

    # Prepare nnU-Net input
    if not args.no_nnunet_input and successful > 0:
        print(f"\n{'=' * 60}")
        print("  Preparing nnU-Net INPUT folder")
        print("=" * 60)
        prepare_nnunet_input(str(data_dir), str(nnunet_input_dir))

    # Summary
    print(f"\n{'=' * 60}")
    print("  Summary")
    print("=" * 60)
    print(f"  Successful: {successful}/{len(patient_folders)}")
    if failed:
        print(f"  Failed:     {', '.join(failed)}")
    print(f"\n  DATA directory: {data_dir}")
    if not args.no_nnunet_input and successful > 0:
        print(f"  nnU-Net INPUT:  {nnunet_input_dir}")

    print(f"\n  Next steps:")
    print(f"    1. Run TotalSegmentator:")
    print(f"       python run_TotalSegmentator.py --data-dir {data_dir}")
    print(f"    2. Run nnU-Net (optional):")
    print(f"       python run_model191.py <nnunet_base_dir>")
    print(f"    3. Run MATLAB pipeline:")
    print(f"       Set DATA_DIRECTORY = '{data_dir}' in LungScape.m")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
