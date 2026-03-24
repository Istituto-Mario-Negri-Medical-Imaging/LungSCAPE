#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_data_ready.py - LungScape pre-MATLAB data validation

Scans the DATA/ directory and verifies that every patient folder contains
the files required by LungScape.m (loadPatientData / loadPatientVolumes).

Usage:
    python check_data_ready.py <data_dir>

Arguments:
    data_dir    Path to the DATA/ directory produced by the Python pipeline.

Exit codes:
    0   All patients are READY (no required file is missing).
    1   One or more patients are BLOCKED (required file missing), or
        data_dir does not exist / contains no patient folders.

Output example:
    Patient001  [READY]    CT  lungs  lobes      arteries  veins  airways(201)  highAtten
    Patient002  [BLOCKED]  CT  lungs  -(no lobes) arteries  veins  ✗ airways missing
    Patient003  [READY]    CT  lungs  -(no lobes) arteries  veins  airways(TS)

    Summary: 2 ready, 1 blocked.
    Blocked patients:
      Patient002 — airways segmentation missing
                   (expected airways201 from run_model201.py
                    or airwaysTS from run_TotalSegmentator.py --use-ts-airways)
"""

import sys
from pathlib import Path


# ── File definitions ──────────────────────────────────────────────────────────

REQUIRED = [
    ("ct",       "{id}.nii.gz",            "CT volume"),
    ("lungs",    "{id}_lungsTS.nii.gz",     "lung mask (TotalSegmentator)"),
    ("arteries", "{id}_arteriesTS.nii.gz",  "arteries (TotalSegmentator >= 2.13.0)"),
    ("veins",    "{id}_veinsTS.nii.gz",     "veins (TotalSegmentator >= 2.13.0)"),
]

# Airways: Model201 preferred, airwaysTS accepted as fallback
AIRWAYS_201 = "{id}_airways201.nii.gz"
AIRWAYS_TS  = "{id}_airwaysTS.nii.gz"

OPTIONAL = [
    ("lobes",       "{id}_lobesTS.nii.gz",             "lobe segmentation"),
    ("highAtten",   "{id}_highAttenuation191.nii.gz",   "high attenuation (Model191)"),
    ("consolidation", "{id}_consolidation.nii.gz",      "external consolidation mask"),
]


# ── Per-patient check ─────────────────────────────────────────────────────────

def check_patient(patient_dir: Path) -> dict:
    """
    Check one patient folder and return a result dict with keys:
        name        patient ID (folder name)
        ready       bool — all required files present
        missing     list of (key, description) for missing required files
        airways_src 'model201' | 'ts' | None
        optional    dict key -> bool
    """
    pid = patient_dir.name
    result = {
        "name":       pid,
        "ready":      True,
        "missing":    [],
        "airways_src": None,
        "optional":   {},
    }

    # Required files
    for key, pattern, description in REQUIRED:
        path = patient_dir / pattern.format(id=pid)
        if not path.exists():
            result["ready"] = False
            result["missing"].append((key, description))

    # Airways (at least one of the two)
    has_201 = (patient_dir / AIRWAYS_201.format(id=pid)).exists()
    has_ts  = (patient_dir / AIRWAYS_TS.format(id=pid)).exists()
    if has_201:
        result["airways_src"] = "model201"
    elif has_ts:
        result["airways_src"] = "ts"
    else:
        result["ready"] = False
        result["missing"].append((
            "airways",
            f"airways segmentation missing\n"
            f"           (expected airways201 from run_model201.py\n"
            f"            or airwaysTS from run_TotalSegmentator.py --use-ts-airways)"
        ))

    # Optional files
    for key, pattern, _ in OPTIONAL:
        result["optional"][key] = (patient_dir / pattern.format(id=pid)).exists()

    return result


# ── Formatting helpers ────────────────────────────────────────────────────────

def _tag(present: bool, label: str) -> str:
    return label if present else f"-(no {label})"

def _airways_tag(src) -> str:
    if src == "model201":
        return "airways(201)"
    if src == "ts":
        return "airways(TS) "
    return "✗ airways   "

def format_patient_line(r: dict) -> str:
    status = "[READY]  " if r["ready"] else "[BLOCKED]"
    opt = r["optional"]
    cols = [
        "CT",
        "lungs",
        _tag(opt.get("lobes",      False), "lobes     "),
        "arteries",
        "veins",
        _airways_tag(r["airways_src"]),
        _tag(opt.get("highAtten",  False), "highAtten "),
        _tag(opt.get("consolidation", False), "consolidation"),
    ]
    return f"  {r['name']:<20} {status}  {'  '.join(cols)}"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0)

    data_dir = Path(sys.argv[1])

    if not data_dir.exists():
        print(f"[ERROR] DATA directory not found: {data_dir}")
        sys.exit(1)

    patient_dirs = sorted(p for p in data_dir.iterdir() if p.is_dir())

    if not patient_dirs:
        print(f"[ERROR] No patient folders found in {data_dir}")
        sys.exit(1)

    print()
    print("=" * 90)
    print("  LungScape — Pre-MATLAB Data Validation")
    print(f"  DATA directory: {data_dir}")
    print("=" * 90)
    print()
    print(f"  {'Patient':<20} {'Status':<11}  "
          f"{'CT':<4}  {'lungs':<7}  {'lobes':<10}  "
          f"{'arteries':<10}  {'veins':<7}  {'airways':<14}  "
          f"{'highAtten':<11}  consolidation")
    print("  " + "-" * 86)

    results = []
    for pd in patient_dirs:
        r = check_patient(pd)
        results.append(r)
        print(format_patient_line(r))

    n_ready   = sum(1 for r in results if r["ready"])
    n_blocked = sum(1 for r in results if not r["ready"])

    print()
    print("=" * 90)
    print(f"  Summary: {n_ready} ready, {n_blocked} blocked  "
          f"(total: {len(results)} patients)")

    if n_blocked > 0:
        print()
        print("  Blocked patients:")
        for r in results:
            if not r["ready"]:
                for _, desc in r["missing"]:
                    print(f"    {r['name']} — {desc}")
        print()
        print("  Run the missing pipeline steps and re-run check_data_ready.py to confirm.")
        print("=" * 90)
        print()
        sys.exit(1)
    else:
        print()
        print("  All patients are ready. You can now run the MATLAB pipeline:")
        print(f"    Set DATA_DIRECTORY = '{data_dir}' in LungScape.m and run LungScape.")
        print("=" * 90)
        print()
        sys.exit(0)


if __name__ == "__main__":
    main()
