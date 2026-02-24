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


def process_patient(ct_path, patient_dir, patient_id):
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
    os.rename(
        os.path.join(patient_dir, "lung_trachea_bronchia.nii.gz"),
        os.path.join(patient_dir, f"{patient_id}_airwaysTS.nii.gz"),
    )

    # --- Remove intermediate per-lobe files ---
    for filename, _ in LOBE_FILES:
        path = os.path.join(patient_dir, filename)
        if os.path.exists(path):
            os.remove(path)

    print(f"  Saved: {patient_id}_lobesTS / _lungsTS / _vesselsTS / _airwaysTS")


def main():
    parser = argparse.ArgumentParser(
        description="Run TotalSegmentator lung segmentation and aggregate lobes."
    )
    parser.add_argument("input_dir", help="Directory containing NIfTI CT images")
    parser.add_argument("output_dir", help="Directory where patient folders will be created")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    nifti_files = sorted(
        f for f in os.listdir(args.input_dir)
        if f.endswith((".nii", ".nii.gz"))
    )

    if not nifti_files:
        print(f"No NIfTI files found in {args.input_dir}")
        return

    for ct_file in nifti_files:
        patient_id = ct_file.split(".nii")[0]
        ct_path = os.path.join(args.input_dir, ct_file)
        patient_dir = os.path.join(args.output_dir, patient_id)

        print(f"Processing {patient_id}...")
        process_patient(ct_path, patient_dir, patient_id)
        print(f"Done: {patient_id}\n")


if __name__ == "__main__":
    main()
