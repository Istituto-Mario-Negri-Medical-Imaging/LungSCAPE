import os
import argparse
import subprocess
import nibabel as nib
import numpy as np


LOBE_FILES = [
    ("lung_upper_lobe_left.nii.gz", 1),
    ("lung_lower_lobe_left.nii.gz", 2),
    ("lung_upper_lobe_right.nii.gz", 3),
    ("lung_middle_lobe_right.nii.gz", 4),
    ("lung_lower_lobe_right.nii.gz", 5),
]

ROI_SUBSET = [name.replace(".nii.gz", "") for name, _ in LOBE_FILES]


def aggregate_lobes(patient_dir):
    """Merge individual lobe masks into a single labelmap (labels 1-5)
    and a binary whole-lung mask."""
    first_path = os.path.join(patient_dir, LOBE_FILES[0][0])
    nii = nib.load(first_path)
    affine = nii.affine
    result = np.zeros(nii.shape, dtype=np.uint8)

    for filename, label in LOBE_FILES:
        mask = nib.load(os.path.join(patient_dir, filename)).get_fdata().astype(np.uint8)
        result[mask != 0] = label

    lobes = nib.Nifti1Image(result, affine)
    lungs = nib.Nifti1Image((result > 0).astype(np.uint8), affine)
    return lobes, lungs


def process_patient(ct_path, patient_dir, patient_id, use_ts_airways=False):
    """Run TotalSegmentator (lobes + lung_vessels) and produce final outputs."""
    os.makedirs(patient_dir, exist_ok=True)

    # --- Run TotalSegmentator ---
    subprocess.run([
        "TotalSegmentator",
        "-i", ct_path,
        "-o", patient_dir,
        "--roi_subset", *ROI_SUBSET,
    ], check=True)

    subprocess.run([
        "TotalSegmentator",
        "-i", ct_path,
        "-o", patient_dir,
        "--ta", "lung_vessels",
    ], check=True)

    # --- Build aggregated volumes ---
    lobes_nii, lungs_nii = aggregate_lobes(patient_dir)

    # --- Save final outputs with patient-prefixed names ---
    nib.save(lobes_nii, os.path.join(patient_dir, f"{patient_id}_lobesTS.nii.gz"))
    nib.save(lungs_nii, os.path.join(patient_dir, f"{patient_id}_lungsTS.nii.gz"))

    os.rename(
        os.path.join(patient_dir, "lung_vessels.nii.gz"),
        os.path.join(patient_dir, f"{patient_id}_vesselsTS.nii.gz"),
    )

    trachea_bronchia_path = os.path.join(patient_dir, "lung_trachea_bronchia.nii.gz")
    if use_ts_airways:
        os.rename(
            trachea_bronchia_path,
            os.path.join(patient_dir, f"{patient_id}_airwaysTS.nii.gz"),
        )
    else:
        # Airways will be provided by nnU-Net Model201; discard TotalSegmentator output
        if os.path.exists(trachea_bronchia_path):
            os.remove(trachea_bronchia_path)

    # --- Remove intermediate per-lobe files ---
    for filename, _ in LOBE_FILES:
        path = os.path.join(patient_dir, filename)
        if os.path.exists(path):
            os.remove(path)

    airways_source = "_airwaysTS" if use_ts_airways else "(airways from Model201, not TS)"
    print(f"  Saved: {patient_id}_lobesTS / _lungsTS / _vesselsTS / {airways_source}")


def find_patients_from_data_dir(data_dir):
    """Find CT NIfTI files inside a DATA/ directory structure.

    Expects: DATA/PatientID/PatientID.nii.gz
    Returns list of (ct_path, patient_dir, patient_id) tuples.
    """
    patients = []
    for entry in sorted(os.scandir(data_dir), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        patient_id = entry.name
        patient_dir = entry.path
        # Look for PatientID.nii.gz or PatientID.nii
        for ext in (".nii.gz", ".nii"):
            ct_path = os.path.join(patient_dir, f"{patient_id}{ext}")
            if os.path.isfile(ct_path):
                patients.append((ct_path, patient_dir, patient_id))
                break
        else:
            print(f"  [WARN] No CT found for {patient_id}, skipping")
    return patients


def find_patients_from_flat_dir(input_dir, output_dir):
    """Find CT NIfTI files in a flat directory (original behavior).

    Returns list of (ct_path, patient_dir, patient_id) tuples.
    """
    nifti_files = sorted(
        f for f in os.listdir(input_dir)
        if f.endswith((".nii", ".nii.gz"))
    )
    patients = []
    for ct_file in nifti_files:
        patient_id = ct_file.split(".nii")[0]
        ct_path = os.path.join(input_dir, ct_file)
        patient_dir = os.path.join(output_dir, patient_id)
        patients.append((ct_path, patient_dir, patient_id))
    return patients


def main():
    parser = argparse.ArgumentParser(
        description="Run TotalSegmentator lung segmentation and aggregate lobes."
    )
    # Two modes: --data-dir (new) or input_dir + output_dir (legacy)
    parser.add_argument("input_dir", nargs="?", default=None,
                        help="Directory containing NIfTI CT images (legacy mode)")
    parser.add_argument("output_dir", nargs="?", default=None,
                        help="Directory where patient folders will be created (legacy mode)")
    parser.add_argument("--data-dir",
                        help="DATA directory with patient subfolders (from prepare_dataset.py). "
                             "Each subfolder must contain PatientID.nii.gz.")
    parser.add_argument("--use-ts-airways",
                        action="store_true",
                        help="Use TotalSegmentator for airways segmentation instead of nnU-Net Model201 "
                             "(saves {PatientID}_airwaysTS.nii.gz). By default airways are discarded "
                             "here and should be generated with run_model201.py.")
    args = parser.parse_args()

    # Determine mode
    if args.data_dir:
        if not os.path.isdir(args.data_dir):
            print(f"[ERROR] DATA directory not found: {args.data_dir}")
            return
        patients = find_patients_from_data_dir(args.data_dir)
    elif args.input_dir and args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        patients = find_patients_from_flat_dir(args.input_dir, args.output_dir)
    else:
        parser.error("Provide either --data-dir or both input_dir and output_dir")
        return

    if not patients:
        print("No patients found to process.")
        return

    for ct_path, patient_dir, patient_id in patients:
        print(f"Processing {patient_id}...")
        process_patient(ct_path, patient_dir, patient_id, use_ts_airways=args.use_ts_airways)
        print(f"Done: {patient_id}\n")


if __name__ == "__main__":
    main()
