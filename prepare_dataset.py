#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LungScape - Step 0: Dataset Preparation

Converts a dataset (DICOM or NIfTI) to a standard NIfTI format, resamples to
the target resolution, and organizes the working directory structure expected
by the LungScape pipeline.

Usage:
    # From DICOM (default):
    python prepare_dataset.py <dicom_dir> <output_dir> [options]

    # From NIfTI:
    python prepare_dataset.py <nifti_dir> <output_dir> --nifti [options]

Input formats:
    DICOM mode (default):
        <dicom_dir>/
        ├── Patient001/      (folder with DICOM files)
        ├── Patient002/
        └── ...

    NIfTI mode (--nifti):
        Flat layout:
            <nifti_dir>/
            ├── Patient001.nii.gz
            ├── Patient002.nii.gz
            └── ...
        Subdir layout:
            <nifti_dir>/
            ├── Patient001/
            │   └── Patient001.nii.gz   (or any single .nii/.nii.gz file)
            ├── Patient002/
            └── ...

Output structure:
    <output_dir>/
    ├── DATA/
    │   ├── Patient001/
    │   │   └── Patient001.nii.gz
    │   ├── Patient002/
    │   │   └── Patient002.nii.gz
    │   └── ...
    └── nnUNet_INPUT/          (unless --no-nnunet-input)
        ├── case001_0000.nii.gz
        ├── case002_0000.nii.gz
        ├── ...
        └── case_mapping.json

Requirements:
    - dcm2niix (system command, must be in PATH) -- DICOM mode only
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
import numpy as np


# -- Constants ----------------------------------------------------------------

TARGET_MATRIX = 512         # Target in-plane matrix size (512 x 512)
TARGET_Z_SPACING = 0.5      # Target slice spacing in mm


# -- Conversion helpers -------------------------------------------------------

def ants_to_nibabel(ants_img):
    """Convert an ANTsImage to a nibabel Nifti1Image without temporary files.

    Builds the NIfTI affine directly from ANTs spacing/direction/origin,
    converting from ANTs/ITK LPS convention to NIfTI RAS convention.
    """
    data = ants_img.numpy()
    spacing = np.array(ants_img.spacing)
    direction = np.array(ants_img.direction).reshape(3, 3)
    origin = np.array(ants_img.origin)

    # Affine in ANTs/ITK LPS convention: each column j of [:3,:3] is the
    # direction of image axis j scaled by its voxel spacing.
    affine = np.eye(4)
    affine[:3, :3] = direction * spacing
    affine[:3, 3] = origin

    # LPS (ANTs/ITK) -> RAS (NIfTI): negate the x and y rows
    affine[0] *= -1
    affine[1] *= -1

    return nib.Nifti1Image(data, affine)


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

def process_patient_dicom(dicom_folder, data_dir, patient_id, skip_resampling=False):
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


def process_patient_nifti(nii_path, data_dir, patient_id, skip_resampling=False):
    """Process a single patient from an existing NIfTI: resample -> save.

    Creates ``data_dir/patient_id/patient_id.nii.gz``.
    Returns the path to the saved NIfTI, or None on failure.
    """
    nii_path = str(nii_path)
    print(f"\n  Loading NIfTI: {os.path.basename(nii_path)} ...")

    if skip_resampling:
        print(f"  Skipping resampling (--skip-resampling)")
        result_img = nib.load(nii_path)
    else:
        result_img = resample_volume(nii_path)

    patient_dir = os.path.join(data_dir, patient_id)
    os.makedirs(patient_dir, exist_ok=True)
    output_path = os.path.join(patient_dir, f"{patient_id}.nii.gz")
    nib.save(result_img, output_path)
    print(f"  Saved: {output_path}")

    return output_path


def find_nifti_patients(nifti_dir):
    """Find NIfTI patient files in a directory.

    Supports two layouts:
    - Flat:   <nifti_dir>/Patient001.nii.gz  (patient_id = filename stem)
    - Subdir: <nifti_dir>/Patient001/        (patient_id = subdir name;
              looks for PatientID.nii.gz first, then any .nii.gz/.nii file)

    Returns a sorted list of (nii_path, patient_id) tuples.
    """
    nifti_dir = Path(nifti_dir)

    def is_nifti(p):
        return p.is_file() and (p.name.endswith(".nii.gz") or p.name.endswith(".nii"))

    # Try flat layout first
    flat_files = sorted(f for f in nifti_dir.iterdir() if is_nifti(f))
    if flat_files:
        patients = []
        for f in flat_files:
            # Strip .nii.gz or .nii extension to get patient_id
            patient_id = f.name
            for ext in (".nii.gz", ".nii"):
                if patient_id.endswith(ext):
                    patient_id = patient_id[: -len(ext)]
                    break
            patients.append((f, patient_id))
        return patients

    # Subdir layout
    patients = []
    for subdir in sorted(nifti_dir.iterdir()):
        if not subdir.is_dir():
            continue
        patient_id = subdir.name

        # Prefer PatientID.nii.gz / PatientID.nii
        found = None
        for ext in (".nii.gz", ".nii"):
            candidate = subdir / f"{patient_id}{ext}"
            if candidate.exists():
                found = candidate
                break

        # Fall back to any single NIfTI file in the subdir
        if found is None:
            candidates = sorted(f for f in subdir.iterdir() if is_nifti(f))
            if not candidates:
                print(f"  [WARN] No NIfTI file found in {subdir}, skipping")
                continue
            if len(candidates) > 1:
                print(f"  [WARN] Multiple NIfTI files in {subdir.name}/, "
                      f"using: {candidates[0].name}")
            found = candidates[0]

        patients.append((found, patient_id))

    return patients


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
        description=(
            "LungScape Step 0: Prepare dataset (DICOM or NIfTI) for the pipeline. "
            "Resamples volumes to 512x512 / 0.5 mm z-spacing and organizes the "
            "DATA/ and nnUNet_INPUT/ directory structure."
        )
    )
    parser.add_argument(
        "input_dir",
        help=(
            "Input directory. DICOM mode (default): directory with one subfolder per "
            "patient, each containing DICOM files. NIfTI mode (--nifti): directory "
            "containing .nii.gz/.nii files (flat) or patient subfolders with NIfTI files."
        ),
    )
    parser.add_argument(
        "output_dir",
        help="Output base directory (DATA/ and nnUNet_INPUT/ will be created here)",
    )
    parser.add_argument(
        "--nifti",
        action="store_true",
        help="Input directory contains NIfTI files instead of DICOM folders",
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

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    data_dir = output_dir / "DATA"
    nnunet_input_dir = output_dir / "nnUNet_INPUT"

    if not input_dir.is_dir():
        print(f"[ERROR] Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("  LungScape - Step 0: Dataset Preparation")
    print("=" * 60)

    # ------------------------------------------------------------------ NIfTI
    if args.nifti:
        patients = find_nifti_patients(input_dir)
        if not patients:
            print(f"[ERROR] No NIfTI files found in {input_dir}", file=sys.stderr)
            sys.exit(1)

        print(f"\n  Input (NIfTI): {input_dir}")
        print(f"  Output:        {output_dir}")
        print(f"  Patients:      {len(patients)}")
        print(f"  Resampling:    {'OFF' if args.skip_resampling else f'{TARGET_MATRIX}x{TARGET_MATRIX}, z={TARGET_Z_SPACING}mm'}")
        print(f"  nnU-Net INPUT: {'OFF' if args.no_nnunet_input else 'ON'}")

        os.makedirs(data_dir, exist_ok=True)
        successful = 0
        failed = []

        for idx, (nii_path, patient_id) in enumerate(patients, start=1):
            print(f"\n[{idx}/{len(patients)}] {patient_id}")
            print("-" * 40)
            try:
                result = process_patient_nifti(
                    nii_path, str(data_dir), patient_id,
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

    # ------------------------------------------------------------------ DICOM
    else:
        # Check dcm2niix is available
        if shutil.which("dcm2niix") is None:
            print("[ERROR] dcm2niix not found in PATH. Install it first.", file=sys.stderr)
            print("  Ubuntu/Debian: sudo apt install dcm2niix", file=sys.stderr)
            print("  conda: conda install -c conda-forge dcm2niix", file=sys.stderr)
            sys.exit(1)

        patient_folders = sorted(d for d in input_dir.iterdir() if d.is_dir())
        if not patient_folders:
            print(f"[ERROR] No patient subfolders found in {input_dir}", file=sys.stderr)
            sys.exit(1)

        print(f"\n  Input (DICOM): {input_dir}")
        print(f"  Output:        {output_dir}")
        print(f"  Patients:      {len(patient_folders)}")
        print(f"  Resampling:    {'OFF' if args.skip_resampling else f'{TARGET_MATRIX}x{TARGET_MATRIX}, z={TARGET_Z_SPACING}mm'}")
        print(f"  nnU-Net INPUT: {'OFF' if args.no_nnunet_input else 'ON'}")

        os.makedirs(data_dir, exist_ok=True)
        successful = 0
        failed = []

        for idx, patient_folder in enumerate(patient_folders, start=1):
            patient_id = patient_folder.name
            print(f"\n[{idx}/{len(patient_folders)}] {patient_id}")
            print("-" * 40)
            try:
                result = process_patient_dicom(
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

    # ------------------------------------------------------------------ Common
    n_patients = len(patients) if args.nifti else len(patient_folders)

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
    print(f"  Successful: {successful}/{n_patients}")
    if failed:
        print(f"  Failed:     {', '.join(failed)}")
    print(f"\n  DATA directory: {data_dir}")
    if not args.no_nnunet_input and successful > 0:
        print(f"  nnU-Net INPUT:  {nnunet_input_dir}")

    print(f"\n  Next steps:")
    print(f"    1. Run TotalSegmentator (lobes, lungs, vessels):")
    print(f"       python run_TotalSegmentator.py --data-dir {data_dir}")
    print(f"    2. Run nnU-Net Model201 (airways segmentation):")
    print(f"       python setup_model201.py <model201_base_dir>   # first time only")
    print(f"       python run_model201.py <model201_base_dir> --project-dir {output_dir}")
    print(f"       (env vars are set internally; outputs deployed to DATA/ via case_mapping.json)")
    print(f"    3. Run nnU-Net Model191 (high attenuation):")
    print(f"       python setup_model191.py <model191_base_dir>   # first time only")
    print(f"       python run_model191.py <model191_base_dir> --project-dir {output_dir}")
    print(f"       (env vars are set internally; outputs deployed to DATA/ via case_mapping.json)")
    print(f"    4. Run MATLAB pipeline:")
    print(f"       Set DATA_DIRECTORY = '{data_dir}' in LungScape.m")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
