# LungScape - Lung Segmentation: a Comprehensive and Automatic PipelinE
Alberto Arrigoni - Istituto di Ricerche Farmacologiche Mario Negri IRCCS

![Logo](images/lungscape_logo.png)

## Overview

LungScape is a two-stage pipeline for comprehensive segmentation of lung structures from chest HRCT scans. Deep learning models (TotalSegmentator, nnU-Net) generate initial segmentations in Python, which are then refined and combined by a MATLAB processing engine using morphological operations, vesselness filtering, and graph-based analysis.

The pipeline segments:
- **Airways** and trachea (with wall detection)
- **Pulmonary vessels** classified by diameter (large, medium, small), with wall detection
- **Lung lobes** and interlobar fissures
- **Pathological tissues**: reticulation/consolidations, ground glass opacities (GGO), air trapping
- **Regional anatomical analysis**: Superior-Inferior, Posterior-Anterior, Left-Right divisions, and lobe-based analysis

## Pipeline Architecture

```
CT NIfTI volumes
       |
       v
 +--------------------------+
 | PYTHON PREPROCESSING     |
 |  1. TotalSegmentator     |  --> lobes, lungs, vessels, airways masks
 |  2. nnU-Net Model191     |  --> high attenuation abnormalities mask
 +--------------------------+
       |
       v
 +--------------------------+
 | MATLAB PROCESSING        |
 |  LungScape.m (16 steps)  |  --> refined segmentations, label maps,
 |                          |      vessel classification, regional maps
 +--------------------------+
       |
       v
 NRRD + MAT output files
```

## Project Structure

```
LungScape/
├── LungScape.m                     # Main MATLAB pipeline (16 steps)
├── README.md
├── requirements.txt                # Python dependencies
├── run_TotalSegmentator.py         # Step 1: generate initial segmentations
├── setup_model191.py               # Step 2a: download and install nnU-Net model
├── run_model191.py                 # Step 2b: run nnU-Net high attenuation segmentation
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
    ├── createRegionalMaps.m              # SI, PA, LR regional division
    └── createLabelMap.m                  # Final label map generation
```

## Requirements

### Python

Install with:
```bash
pip install -r requirements.txt
```

| Package | Version | Purpose |
|---------|---------|---------|
| TotalSegmentator | >= 2.0 | Lung, lobe, vessel, airway segmentation |
| nibabel | >= 4.0 | NIfTI file I/O |
| numpy | >= 1.21 | Array operations |
| nnunetv2 | >= 2.0 | High attenuation abnormality segmentation |

TotalSegmentator will automatically install PyTorch and its own dependencies.

### MATLAB

- **MATLAB R2019b** or later
- **Toolboxes**: Image Processing Toolbox, Parallel Computing Toolbox
- **External libraries** (must be on the MATLAB path):
  - [NIfTI and ANALYZE tools](https://mathworks.com/matlabcentral/fileexchange/8797) (`load_untouch_nii`)
  - `anisodiff3D` - 3D anisotropic diffusion filter
  - `vesselness3D` - Jerman vessel enhancement filter (Jerman et al., IEEE TMI 2016)
  - `nrrdWriter` - NRRD file export utility
  - [Image Graphs](https://it.mathworks.com/matlabcentral/fileexchange/53614-image-graphs) - pixel neighbor graph analysis

## Quick Start

### Step 1 - Python Preprocessing: TotalSegmentator

Run TotalSegmentator on your CT volumes to generate the initial segmentation masks required by the MATLAB pipeline.

```bash
python run_TotalSegmentator.py <input_dir> <output_dir>
```

- `<input_dir>`: directory containing the CT NIfTI files (`.nii` or `.nii.gz`)
- `<output_dir>`: directory where patient folders will be created (this will become the `DATA` folder for MATLAB)

For each CT volume, the script:
1. Runs TotalSegmentator with lobe-specific ROI extraction (5 lobes)
2. Runs TotalSegmentator with the `lung_vessels` task for vessel and airway segmentation
3. Aggregates the individual lobe masks into a single labelmap (labels 1-5)
4. Creates a binary whole-lung mask
5. Renames all outputs with patient-prefixed names

**Output files per patient:**

| File | Description |
|------|-------------|
| `{PatientID}_lobesTS.nii.gz` | Lobe segmentation (1=upper left, 2=lower left, 3=upper right, 4=middle right, 5=lower right) |
| `{PatientID}_lungsTS.nii.gz` | Binary whole-lung mask |
| `{PatientID}_vesselsTS.nii.gz` | Pulmonary vessel segmentation |
| `{PatientID}_airwaysTS.nii.gz` | Airway and trachea segmentation |

### Step 2 - Python Preprocessing: nnU-Net High Attenuation Model

This step segments high attenuation abnormalities (consolidations, dense opacities) using a pretrained nnU-Net v2 model (Dataset191). Skip this step if your patients do not present significant pathology.

**2a. Setup the model (first time only):**

```bash
python setup_model191.py <base_dir>
```

This will:
- Download the pretrained model from Zenodo
- Create the nnU-Net directory structure (`nnUNet_raw/`, `nnUNet_preprocessed/`, `nnUNet_results/`, `INPUT/`, `OUTPUT/`)
- Install the model
- Generate an activation script for environment variables

**2b. Run the segmentation:**

Place CT scans in `<base_dir>/INPUT/` with nnU-Net naming convention (`case001_0000.nii.gz`, `case002_0000.nii.gz`, ...), then:

```bash
# Load nnU-Net environment variables (in a new terminal)
source <base_dir>/activate_nnunet_paths.sh

python run_model191.py <base_dir>
```

Use `--no-cleanup` to keep nnU-Net folders for reuse:
```bash
python run_model191.py <base_dir> --no-cleanup
```

**Important**: the output must be renamed to `{PatientID}_HighAttenuation.nii.gz` and placed in the corresponding patient folder within the `DATA` directory.

### Step 3 - Organize Data Directory

After preprocessing, the `DATA` directory should have this structure:

```
DATA/
├── Patient001/
│   ├── Patient001.nii.gz                   (CT volume - required)
│   ├── Patient001_lobesTS.nii.gz           (from TotalSegmentator - optional)
│   ├── Patient001_lungsTS.nii.gz           (from TotalSegmentator - required)
│   ├── Patient001_vesselsTS.nii.gz         (from TotalSegmentator - required)
│   ├── Patient001_airwaysTS.nii.gz         (from TotalSegmentator - required)
│   ├── Patient001_consolidation.nii.gz     (external consolidation mask - optional)
│   └── Patient001_HighAttenuation.nii.gz   (from nnU-Net Model191 - optional)
├── Patient002/
│   └── ...
```

### Step 4 - Run MATLAB Pipeline

```matlab
% Set the DATA path in LungScape.m, then run:
LungScape
```

The pipeline processes each patient through 16 steps:

| Step | Function | Description |
|------|----------|-------------|
| 1 | `loadPatientVolumes` | Load CT and segmentation volumes |
| 2 | `anisodiff3D` | Anisotropic diffusion filtering of CT |
| 3 | `preprocessLungs` | Erode margins, separate L/R lungs, compute convex hull |
| 4 | `refineLobeSegmentation` | Fill gaps in TotalSegmentator lobe output |
| 5 | `segmentFissures` | Detect interlobar fissures |
| 6 | `classifyDistanceFromBorder` | Classify parenchyma by distance from pleural surface |
| 7 | `refineAirways` | Remove spurious airways, detect honeycombing |
| 8 | `segmentInjuries` | Segment GGO, consolidations, air trapping |
| 9 | `segmentVessels` | Jerman vesselness filter + graph-based reconstruction |
| 10 | `refineVessels` | Multi-threshold refinement, loop removal |
| 11 | `segmentTrachea` | Trachea identification and wall separation |
| 12 | `refineAirSeg` | Refine low-attenuation segmentation |
| 13 | `finalizeSegmentations` | Wall detection (airways, vessels, trachea) |
| 14 | `classifyVesselsByDiameter` | Classify vessels: large, medium, small |
| 15 | `createRegionalMaps` | Create SI, PA, LR regional maps |
| 16 | `createLabelMap` + `saveResults` | Combine and save final outputs |

### Step 5 - View Results

Results are saved in `<PatientFolder>/InterimResults/`:

| File | Description |
|------|-------------|
| `TotalLabelMap.nrrd` | Complete segmentation (labels 0-10, see table below) |
| `VesselsLabelMap.nrrd` | Vessel diameter classification (1=large, 2=medium, 3=small) |
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
| **Vessels** | Jerman scales (0.5-1.5 mm), tau (0.5), FA thresholds (0.92-0.95), large vessel cutoff (> 8 mm) |
| **Distance classification** | Close/midfar/far thresholds for multi-zone refinement |
| **Morphology** | Structuring element sizes (voxel-aware, auto-adjusted for scanner resolution) |

## Citation

If you use this code, please cite this github page!

## License

[Specify license]

## Contact

For questions or issues:
- Alberto Arrigoni

## Changelog

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
