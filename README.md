# LungSCAPE - Lung Segmentation: a Comprehensive and Automatic PipelinE
Alberto Arrigoni - Istituto di Ricerche Farmacologiche Mario Negri IRCCS

<p align="center">
  <img src="images/lungscape_logo.png" alt="Logo" width="600">
</p>

## Overview

LungScape is a two-stage pipeline for comprehensive segmentation of lung structures from chest HRCT scans. Deep learning models (TotalSegmentator, nnU-Net) generate initial segmentations in Python, which are then refined and combined by a MATLAB processing engine using morphological operations, vesselness filtering, and graph-based analysis.

The pipeline segments:
- **Airways** and trachea (with wall detection)
- **Pulmonary vessels** classified by diameter (large, medium, small) and by type (arteries, veins), with wall detection
- **Lung lobes** and interlobar fissures
- **Pathological tissues**: reticulation/consolidations, ground glass opacities (GGO), air trapping
- **Regional anatomical analysis**: Superior-Inferior, Posterior-Anterior, Left-Right divisions, and lobe-based analysis

## Pipeline Architecture

```
DICOM Dataset
       |
       v
 +--------------------------+
 | STEP 0: DATASET PREP     |
 |  prepare_dataset.py      |  --> NIfTI conversion, resampling,
 |                          |      DATA/ structure, nnU-Net INPUT/
 +--------------------------+
       |
       v
 +--------------------------+
 | PYTHON PREPROCESSING     |
 |  1. TotalSegmentator     |  --> lobes, lungs, artery/vein masks
 |  2. nnU-Net Model201     |  --> airways mask (default)
 |  [or TS --use-ts-airways]|
 |  3. nnU-Net Model191     |  --> high attenuation abnormalities mask
 |     (optional)           |
 +--------------------------+
       |
       v
 +--------------------------+
 | MATLAB PROCESSING        |
 |  LungScape.m (15 steps)  |  --> refined segmentations, label maps,
 |                          |      vessel classification, regional maps
 +--------------------------+
       |
       v
 NRRD + MAT output files
```

## Project Structure

```
LungScape/
├── LungScape.m                     # Main MATLAB pipeline (17 steps)
├── setup_matlab.m                  # MATLAB environment setup & dependency checker
├── run_pipeline.sh                 # Python pipeline orchestrator (Steps 0-3)
├── check_data_ready.py             # Pre-MATLAB data validation (verify DATA/ completeness)
├── README.md
├── requirements.txt                # Python dependencies
├── prepare_dataset.py              # Step 0: DICOM→NIfTI conversion + resampling
├── run_TotalSegmentator.py         # Step 1: generate lobes, lungs, vessels masks
├── setup_model201.py               # Step 2a: download and install nnU-Net airways model
├── run_model201.py                 # Step 2b: run nnU-Net airways segmentation (default)
├── setup_model191.py               # Step 3a: download and install nnU-Net high attenuation model
├── run_model191.py                 # Step 3b: run nnU-Net high attenuation segmentation (optional)
├── config/
│   └── getProcessingParameters.m         # Centralized parameter management
├── io/
│   ├── loadPatientData.m                 # Scan DATA directory, load patient info
│   ├── loadPatientVolumes.m              # Load NIfTI volumes into MATLAB
│   └── saveResults.m                     # Save NRRD and MAT results
├── preprocessing/
│   ├── preprocessLungs.m                 # Erode margins, separate L/R lungs
│   ├── refineLobeSegmentation.m          # Fill gaps in lobe segmentation
│   └── classifyDistanceFromBorder.m      # Distance-based parenchyma classification
├── segmentation/
│   ├── segmentVessels.m                  # Jerman vesselness filter + reconstruction
│   ├── refineVessels.m                   # Multi-threshold vessel refinement
│   ├── classifyVesselsByDiameter.m       # Diameter-based vessel classification
│   ├── refineAirways.m                   # Airway cleanup, honeycombing detection
│   ├── segmentTrachea.m                  # Trachea identification and separation
│   ├── segmentFissures.m                 # Interlobar fissure detection
│   ├── segmentInjuries.m                 # GGO, consolidation, air trapping
│   ├── refineAirSeg.m                    # Low-attenuation refinement
│   └── finalizeSegmentations.m           # Wall segmentation and final masking
└── utils/
    ├── computeFractionalAnisotropy.m     # FA from eigenvalues
    ├── filterByAnisotropy.m              # FA-based shape filtering
    ├── filterByEigenVectorOrientation.m  # Orientation-based filtering
    ├── reconstructBySeeds.m              # Morphological reconstruction
    ├── computeVesselDiameter.m           # Local diameter map computation
    ├── connectVesselBranches.m           # Graph-based branch connectivity
    ├── classifyVesselsByType.m           # Artery/vein classification (topological + distance)
    ├── createRegionalMaps.m              # SI, PA, LR regional division
    ├── createLabelMap.m                  # Final label map generation
    ├── binaryImageGraph3Weighted.m      # Weighted 3-D pixel graph (Image Graphs extension)
    └── binaryImageGraphWeighted.m       # Weighted 2-D pixel graph (Image Graphs extension)
```

## Requirements

### Python

Install with:
```bash
pip install -r requirements.txt
```

| Package | Version | Purpose |
|---------|---------|---------|
| TotalSegmentator | >= 2.13.0 | Lung, lobe, artery/vein/airway segmentation (separate outputs require v2.13.0+) |
| nibabel | >= 4.0 | NIfTI file I/O |
| numpy | >= 1.21 | Array operations |
| nnunetv2 | >= 2.3.1 | High attenuation abnormality segmentation (also installed transitively by TotalSegmentator) |
| antspyx | >= 0.4.0 | Volume resampling (Step 0) |

TotalSegmentator will automatically install PyTorch, nnunetv2, and its own dependencies.

**System dependencies** (must be installed separately):
- [dcm2niix](https://github.com/rordenlab/dcm2niix) - DICOM to NIfTI conversion (Step 0)
  - Ubuntu/Debian: `sudo apt install dcm2niix`
  - conda: `conda install -c conda-forge dcm2niix`

### MATLAB

- **MATLAB R2019b** or later
- **Toolboxes**: Image Processing Toolbox, Parallel Computing Toolbox
- **External libraries** (must be on the MATLAB path):
  - [NIfTI and ANALYZE tools](https://mathworks.com/matlabcentral/fileexchange/8797) (`load_untouch_nii`)
  - [anisodiff3D](https://it.mathworks.com/matlabcentral/fileexchange/14995-anisotropic-diffusion-perona-malik) - 3D anisotropic diffusion filter
  - [vesselness3D](https://it.mathworks.com/matlabcentral/fileexchange/63171-jerman-enhancement-filter) - Jerman vessel enhancement filter (Jerman et al., IEEE TMI 2016). **After installation, the MEX file must be compiled** by running `mex eig3volume.c` from within the library folder in MATLAB (requires a configured C compiler; run `mex -setup` if needed).
  - [nrrdWriter](https://it.mathworks.com/matlabcentral/fileexchange/48621-nrrdwriter-filename-matrix-pixelspacing-origin-encoding) - NRRD file export utility
  - [Image Graphs](https://it.mathworks.com/matlabcentral/fileexchange/53614-image-graphs) - pixel neighbor graph analysis (must be on the MATLAB path). The repository includes `binaryImageGraph3Weighted.m` and `binaryImageGraphWeighted.m` in `utils/`, which are weighted extensions of Image Graphs functions and depend on it.

To verify that all MATLAB dependencies are correctly installed, run the provided setup script:

```matlab
cd /path/to/LungSCAPE_codebase
setup_matlab
```

The script checks each toolbox and library, reports what is missing, and prints the corresponding download URL and installation instructions. Re-run it after resolving each issue to confirm.

## Quick Start

### Option A — Automated (recommended)

Run the entire Python pipeline with a single command using `run_pipeline.sh`:

```bash
# Activate your Python environment first
conda activate lungscape   # or: source venv/bin/activate

./run_pipeline.sh \
  --input  /path/to/dicom \
  --output /path/to/output \
  --model201-dir /path/to/model201 \
  --model191-dir /path/to/model191
```

The script runs Steps 0–3 in sequence, handles model setup automatically on first run (skips download if model is already installed), and saves a full log to `<output_dir>/pipeline.log`.

After completion, open `LungScape.m`, set `DATA_DIRECTORY = '/path/to/output/DATA'`, and run the MATLAB pipeline.

**Key options:**

| Option | Description |
|--------|-------------|
| `--skip-model191` | Skip Model191 (patients without significant high attenuation pathology) |
| `--nifti` | Input is NIfTI instead of DICOM |
| `--skip-resampling` | Skip resampling in Step 0 |
| `--use-ts-airways` | Use TotalSegmentator airways (skips Model201 steps) |
| `--no-cleanup` | Keep nnU-Net folders after inference |

Run `./run_pipeline.sh --help` for the full reference.

---

### Option B — Manual (step by step)

### Step 0 - Dataset Preparation

Use `prepare_dataset.py` to resample volumes to the target resolution and organize the working directory structure. Supports both DICOM and NIfTI input.

**From DICOM (default):**

```bash
python prepare_dataset.py <dicom_dir> <output_dir>
```

- `<dicom_dir>`: directory with one subfolder per patient, each containing DICOM files

For each patient the script:
1. Converts DICOM to NIfTI using `dcm2niix`
2. Resamples the in-plane matrix to 512x512 if needed (ANTsPy, nearest-neighbor)
3. Resamples z-spacing to 0.5 mm (ANTsPy, nearest-neighbor)
4. Saves the result to `DATA/PatientID/PatientID.nii.gz`
5. Prepares `nnUNet_INPUT/` with nnU-Net naming convention (`caseNNN_0000.nii.gz`) and a `case_mapping.json`

**From NIfTI (--nifti):**

```bash
python prepare_dataset.py <nifti_dir> <output_dir> --nifti
```

The NIfTI input directory can have two layouts:

```
# Flat layout (one file per patient):
<nifti_dir>/
├── Patient001.nii.gz
├── Patient002.nii.gz
└── ...

# Subdir layout (one folder per patient):
<nifti_dir>/
├── Patient001/
│   └── Patient001.nii.gz   (or any single .nii/.nii.gz file)
├── Patient002/
└── ...
```

Resampling is applied in both modes unless `--skip-resampling` is set.

**Options:**
- `--nifti`: input directory contains NIfTI files (instead of DICOM folders)
- `--skip-resampling`: keep original matrix size and spacing
- `--no-nnunet-input`: skip nnU-Net INPUT folder preparation

**Output structure:**

```
<output_dir>/
├── DATA/
│   ├── Patient001/
│   │   └── Patient001.nii.gz
│   └── Patient002/
│       └── Patient002.nii.gz
└── nnUNet_INPUT/
    ├── case001_0000.nii.gz
    ├── case002_0000.nii.gz
    └── case_mapping.json
```

### Step 1 - Python Preprocessing: TotalSegmentator

Run TotalSegmentator on the prepared DATA directory (from Step 0) or on a directory of NIfTI files.

```bash
# Using DATA directory from Step 0 (recommended):
python run_TotalSegmentator.py --data-dir <output_dir>/DATA

# Legacy mode (flat directory of NIfTI files):
python run_TotalSegmentator.py <input_dir> <output_dir>
```

When using `--data-dir`, the script finds CT volumes inside patient subfolders (`PatientID/PatientID.nii.gz`) and writes segmentation outputs to the same folders.

For each CT volume, the script:
1. Runs TotalSegmentator with lobe-specific ROI extraction (5 lobes)
2. Runs TotalSegmentator with the `lung_vessels` task for separate artery/vein segmentation (requires TS >= 2.13.0)
3. Aggregates the individual lobe masks into a single labelmap (labels 1-5)
4. Creates a binary whole-lung mask
5. Renames all outputs with patient-prefixed names

By default, TotalSegmentator airways output is discarded (airways are provided by nnU-Net Model201 in Step 2). To use TotalSegmentator for airways instead, pass `--use-ts-airways`:

```bash
python run_TotalSegmentator.py --data-dir <output_dir>/DATA --use-ts-airways
```

**Output files per patient:**

| File | Description |
|------|-------------|
| `{PatientID}_lobesTS.nii.gz` | Lobe segmentation (1=upper left, 2=lower left, 3=upper right, 4=middle right, 5=lower right) |
| `{PatientID}_lungsTS.nii.gz` | Binary whole-lung mask |
| `{PatientID}_arteriesTS.nii.gz` | Pulmonary artery segmentation (TS >= 2.13.0) |
| `{PatientID}_veinsTS.nii.gz` | Pulmonary vein segmentation (TS >= 2.13.0) |
| `{PatientID}_airwaysTS.nii.gz` | Airways (only if `--use-ts-airways` is specified) |

### Step 2 - Python Preprocessing: nnU-Net Airways Model (Model201)

This step segments airways using a pretrained nnU-Net v2 model (Dataset201). This is the **default** airways source. Skip this step only if you ran Step 1 with `--use-ts-airways`.

**2a. Setup the model (first time only):**

```bash
python setup_model201.py <model201_base_dir>
```

This will:
- Download the pretrained airways model from Zenodo
- Create the nnU-Net directory structure (`nnUNet_raw/`, `nnUNet_preprocessed/`, `nnUNet_results/`, `INPUT/`, `OUTPUT/`)
- Install the model
- Generate an activation script for environment variables

**2b. Run the segmentation:**

```bash
python run_model201.py <model201_base_dir> --project-dir <output_dir>
```

The script automatically:
1. Populates `INPUT/` with the files from `nnUNet_INPUT/` (created by Step 0)
2. Sets the required nnU-Net environment variables internally
3. Runs the segmentation
4. Deploys each output to `DATA/PatientID/PatientID_airways201.nii.gz` using `case_mapping.json`

Use `--no-cleanup` to keep nnU-Net folders for reuse:
```bash
python run_model201.py <model201_base_dir> --project-dir <output_dir> --no-cleanup
```

> **Note on `activate_nnunet_paths.sh`:** setup generates this script for manual reference, but running it is not required — `run_model201.py` sets the nnU-Net environment variables (`nnUNet_raw`, `nnUNet_preprocessed`, `nnUNet_results`) internally before invoking the inference. Make sure the correct Python environment (with nnunetv2 and TotalSegmentator installed) is active before running.

### Step 3 - Python Preprocessing: nnU-Net High Attenuation Model (Model191) *(optional)*

This step segments high attenuation abnormalities (consolidations, dense opacities) using a pretrained nnU-Net v2 model (Dataset191). Skip this step if your patients do not present significant pathology.

**3a. Setup the model (first time only):**

```bash
python setup_model191.py <model191_base_dir>
```

**3b. Run the segmentation:**

```bash
python run_model191.py <model191_base_dir> --project-dir <output_dir>
```

The script automatically populates `INPUT/`, sets the nnU-Net environment variables internally, runs the segmentation, and deploys each output to `DATA/PatientID/PatientID_highAttenuation191.nii.gz`. The same note on `activate_nnunet_paths.sh` from Step 2b applies here.

### Step 4 - Verify Data Directory

Before launching the MATLAB pipeline, verify that all required files are in place:

```bash
python check_data_ready.py <output_dir>/DATA
```

This script checks each patient folder for required and optional files, reports their status, and lists any blocking issues. When using `run_pipeline.sh`, this check is run automatically at the end.

Example output:
```
  Patient001           [READY]    CT  lungs  lobes       arteries  veins  airways(201)  highAtten
  Patient002           [BLOCKED]  CT  lungs  -(no lobes) arteries  veins  ✗ airways missing

  Summary: 1 ready, 1 blocked.
  Blocked patients:
    Patient002 — airways segmentation missing
                 (expected airways201 from run_model201.py
                  or airwaysTS from run_TotalSegmentator.py --use-ts-airways)
```

After preprocessing, the `DATA` directory should have this structure:

```
DATA/
├── Patient001/
│   ├── Patient001.nii.gz                   (CT volume - required)
│   ├── Patient001_lobesTS.nii.gz           (from TotalSegmentator - optional)
│   ├── Patient001_lungsTS.nii.gz           (from TotalSegmentator - required)
│   ├── Patient001_arteriesTS.nii.gz        (from TotalSegmentator >= 2.13.0 - required)
│   ├── Patient001_veinsTS.nii.gz           (from TotalSegmentator >= 2.13.0 - required)
│   ├── Patient001_airways201.nii.gz        (from nnU-Net Model201 - required by default)
│   ├── Patient001_airwaysTS.nii.gz         (from TotalSegmentator - required if --use-ts-airways)
│   ├── Patient001_consolidation.nii.gz     (external consolidation mask - optional)
│   └── Patient001_highAttenuation191.nii.gz   (from nnU-Net Model191 - optional)
├── Patient002/
│   └── ...
```

The MATLAB pipeline automatically looks for `airways201` first; if not found, it falls back to `airwaysTS`.

### Step 5 - Run MATLAB Pipeline

```matlab
% Set the DATA path in LungScape.m, then run:
LungScape
```

The pipeline processes each patient through 15 steps (plus two sub-steps):

| Step | Function | Description |
|------|----------|-------------|
| 1 | `loadPatientVolumes` + `anisodiff3D` | Load CT and segmentation volumes; anisotropic diffusion pre-filter |
| 2 | `preprocessLungs` | Reorient volumes, erode margins, separate L/R lungs, compute convex hull |
| 3 | `refineLobeSegmentation` | Fill gaps in TotalSegmentator lobe output |
| 4 | `segmentFissures` | Detect interlobar fissures |
| 5 | `classifyDistanceFromBorder` | Classify parenchyma by distance from pleural surface |
| 6 | `refineAirways` | Remove spurious airways, detect honeycombing |
| 7 | `segmentInjuries` | Segment GGO, consolidations, air trapping |
| 8 | `segmentVessels` | Jerman vesselness filter + graph-based reconstruction |
| 9 | `refineVessels` | Multi-threshold refinement, loop removal |
| 10 | `segmentTrachea` | Trachea identification and wall separation |
| 10b | `postTracheaVesselRefinement` | Second-pass bronchovascular recovery (post-trachea) |
| 11 | `refineAirSeg` | Refine low-attenuation segmentation |
| 12 | `finalizeSegmentations` | Wall detection (airways, vessels, trachea) |
| 13 | `classifyVesselsByDiameter` | Classify vessels by diameter: large, medium, small |
| 13b | `classifyVesselsByType` | Classify vessels by type: arteries, veins, undetermined |
| 14 | `createRegionalMaps` | Create SI, PA, LR regional maps |
| 15 | `createLabelMap` + `saveResults` | Combine and save final outputs |

### Step 6 - View Results

Results are saved in `<PatientFolder>/Results/`:

| File | Description |
|------|-------------|
| `TotalLabelMap.nrrd` | Complete segmentation (labels 0-10, see table below) |
| `VesselsLabelMap.nrrd` | Vessel diameter classification (1=large, 2=medium, 3=small) |
| `VesselsTypeMap.nrrd` | Vessel type classification (1=artery, 2=vein, 3=undetermined) |
| `Vol_<Patient>.nrrd` | Filtered CT volume |
| `<Patient>.mat` | Full MATLAB workspace (all maps, masks, parameters, metadata) |

## Label Map Encoding

| Label | Structure | Description |
|-------|-----------|-------------|
| 0 | Background | Outside lungs |
| 1 | Healthy parenchyma | Normal lung tissue |
| 2 | GGO | Ground glass opacity |
| 3 | Reticulation/Consolidation | Dense opacification |
| 4 | Air trapping | Cysts, honeycombing, air-filled pathology |
| 5 | Vessels | Pulmonary vessels |
| 6 | Airways | Bronchi |
| 7 | Trachea | Main airway |
| 8 | Airway walls | Bronchial walls |
| 9 | Trachea walls | Tracheal walls |
| 10 | Vessel walls | Vascular walls |

## Configuration

All processing hyperparameters are centralized in `config/getProcessingParameters.m` and organized by stage:

| Category | Key parameters |
|----------|---------------|
| **Filtering** | Anisotropic diffusion: kappa, gamma, radius, iterations |
| **HU thresholds** | Air (-980), GGO (-700 to -200), consolidation (< -300), healthy (< -600) |
| **Airways** | Min volume, honeycomb diameter (3.5-5 mm), FA threshold (0.90), wall distance |
| **Vessels** | Jerman scales (0.5-1.5 mm), tau (0.5), FA thresholds (0.92-0.95), large vessel cutoff (> 8 mm); type classification max distance (10 vox) |
| **Distance classification** | Close/midfar/far thresholds for multi-zone refinement |
| **Morphology** | Structuring element sizes (voxel-aware, auto-adjusted for scanner resolution) |

## Citation

If you use this code, please cite this github page!

## License
Apache-2.0 license

## Contact

For questions or issues:
- Alberto Arrigoni

## Changelog

### Version 2.4 (2026-03)
- `run_model201.py` / `run_model191.py`: added `--project-dir` argument (default: current working directory). When `nnUNet_INPUT/` and `DATA/` are found in the project directory, input files are automatically copied from `nnUNet_INPUT/` to `INPUT/` before inference, and outputs are automatically deployed to `DATA/PatientID/` using `case_mapping.json` — no manual renaming or file-moving required.

### Version 2.3 (2026-03)
- Updated TotalSegmentator requirement to >= 2.5.0: `lung_vessels` task now produces separate `lung_arteries.nii.gz` and `lung_veins.nii.gz` outputs
- `run_TotalSegmentator.py`: replaced single vessel rename with two renames (`_arteriesTS`, `_veinsTS`); `lung_airways.nii.gz` replaces `lung_trachea_bronchia.nii.gz`
- `loadPatientData.m` / `loadPatientVolumes.m`: load arteries and veins separately, mask to lung parenchyma, then combine into `volumes.vessels` for backward compatibility
- `segmentVessels.m`: bronchovascular recovery (`recoverVesselsNearAirways`) now constrained to arterial territory; venous voxels in loop/end classification bypass FA filtering via new `filterVenousVessels` helper
- `refineVessels.m` (`addLargeVessels`): replaced HU-threshold approach with TS-based morphological reconstruction from artery+vein masks, removing 7 obsolete parameters
- `postTracheaVesselRefinement.m`: arterial mask constraint applied to bronchovascular recovery
- Added Step 14b: `classifyVesselsByType` — two-phase artery/vein classification (topological propagation via `imreconstruct` from TS anchors, then distance fallback for conflicts and isolated components)
- New output file: `VesselsTypeMap.nrrd` (0=non-vessel, 1=artery, 2=vein, 3=undetermined); `TotalLabelMap.nrrd` label 5 unchanged for backward compatibility
- `saveResults.m`: saves `VesselsTypeMap.nrrd` and adds artery/vein/undetermined volume summary

### Version 2.2 (2026-03)
- Added nnU-Net Model201 for airways segmentation (`setup_model201.py`, `run_model201.py`)
- Model201 is now the **default** airways source; TotalSegmentator airways available via `--use-ts-airways`
- `loadPatientData.m` now prioritizes `airways201` over `airwaysTS` with automatic fallback
- Updated pipeline step numbering and documentation accordingly

### Version 2.1 (2026-02)
- Added Step 0: `prepare_dataset.py` for DICOM-to-NIfTI conversion, resampling (512x512, 0.5mm z-spacing), and automatic working directory setup
- `run_TotalSegmentator.py` now supports `--data-dir` mode to work directly on the DATA/ structure
- Automatic nnU-Net INPUT folder preparation with case mapping

### Version 2.0 (2026-01)
- Major refactoring of original monolithic script into modular architecture
- Added Python preprocessing scripts (`run_TotalSegmentator.py`, `setup_model191.py`, `run_model191.py`)
- Centralized parameter management in `getProcessingParameters.m`
- Complete vessel segmentation pipeline with Jerman filter and graph-based connectivity
- Vessel diameter classification (large/medium/small)
- Fissure segmentation
- Airway, vessel, and trachea wall detection
- Distance-based parenchyma classification for multi-threshold refinement
- Regional anatomical maps (SI, PA, LR)

### Version 1.0
- Original implementation by Alberto Arrigoni
