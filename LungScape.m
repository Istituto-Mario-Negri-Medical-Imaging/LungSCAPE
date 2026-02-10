%% LungScape - Lung CT Segmentation Pipeline
%
% Purpose:
%   Comprehensive segmentation of lung structures from chest CT scans including:
%   - Airways and trachea
%   - Pulmonary vessels (also classified by size)
%   - Lung lobes and fissures
%   - Pathological tissues (reticulations/consolidations, GGO, air trapping)
%   - Regional anatomical analysis
%
% Directory Structure:
%   LungSegMN/
%   ├── config/          - Parameter configuration
%   ├── io/              - Data loading and saving
%   ├── preprocessing/   - Lung mask preprocessing
%   ├── segmentation/    - Core segmentation algorithms
%   └── utils/           - Utility functions (FA, filtering, etc.)
%
% Usage:
%   1. Place patient data in DATA/ directory
%   2. Run this script: LungScape.m
%   3. Results saved to <PatientFolder>/RESULTS/
%
% Input Data Structure (per patient):
%   DATA/Patient001/
%       ├── Patient001/
%       │   ├── Patient001.nii.gz              (CT volume)
%       │   ├── Patient001_lobesTS.nii.gz      (Lobes coarse segmentation)
%       │   ├── Patient001_lungsTS.nii.gz      (Lungs coarse segmentation)
%       │   ├── Patient001_vesselsTS.nii.gz    (Vessels coarse segmentation)
%       │   ├── Patient001_airwaysTS.nii.gz    (Airways coarse segmentation)
%       │   ├── Patient001_consolidation.nii.gz
%       │   ├── Patient001_HighAttenuation.nii.gz
%       └── Patient002/
%           └── ...
%
% Output:
%   <PatientFolder>/RESULTS/
%   ├── LungLabelMap.nrrd  - Complete segmentation (labels 0-10)
%   ├── Vol_<Patient>.nrrd  - Processed CT volume
%   └── <Patient>.mat       - MATLAB workspace with all results
%
% Dependencies:
%   - NIfTI toolbox (load_untouch_nii)
%   - anisodiff3D (anisotropic diffusion filtering)
%   - vesselness3D (Jerman vessel enhancement)
%   - nrrdWriter (NRRD file export)
%
% Author: Alberto Arrigoni, Istituto di Ricerche Farmacologiche Mario Negri IRCCS
%
% See also: loadPatientData, processPatient, getProcessingParameters

% %% Clear workspace and initialize
clc
clear

%% Add paths to all modules
addpath(genpath(pwd));

%% Configuration
DATA_DIRECTORY = '/mnt/sdb1/lungs_longitudinale/LungScapeTEST/DATA';

%% Load patient dataset
fprintf('========================================================================\n');
fprintf('  LungScape - Lung Segmentation: Comprehensive & Automatic PipelinE\n');
fprintf('=========================================================================\n\n');

%% 1 - Load patient data
try
    patientList = loadPatientData(convertCharsToStrings(DATA_DIRECTORY));
catch ME
    error('Failed to load patient data: %s', ME.message);
end

numPatients = length(patientList);
fprintf('Processing %d patients\n\n', numPatients);

%% 2 - Process each patient
for patientIdx = 1:numPatients
    fprintf('======================================\n');
    fprintf('Patient %d/%d: %s\n', patientIdx, numPatients, patientList(patientIdx).name);
    fprintf('======================================\n\n');

    try
        % Process this patient
        processPatient(patientList(patientIdx));

    catch ME
        warning('Error processing patient %s: %s\n%s', ...
            patientList(patientIdx).name, ME.message, getReport(ME));
        fprintf('Skipping to next patient...\n\n');
        continue;
    end

    fprintf('Patient %s completed successfully\n\n', patientList(patientIdx).name);
end

fprintf('======================================\n');
fprintf('All patients processed!\n');
fprintf('======================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Processing Function

function processPatient(patientData)
%PROCESSPATIENT Main processing pipeline for a single patient
%
%   This function orchestrates the entire segmentation pipeline:
%   1. Load volumes and get parameters
%   2. Preprocess CT and lung masks
%   3. Segment structures (airways, vessels, pathology)
%   4. Create regional maps and label map
%   5. Save all results
%
%   See also: loadPatientVolumes, getProcessingParameters, saveResults

    %% Step 1: Load data
    fprintf('----------------------------------------\n');
    fprintf('STEP 1: Loading data...\n');
    fprintf('----------------------------------------\n');
    volumes = loadPatientVolumes(patientData);
    params = getProcessingParameters(volumes.nii);

    %% Step 2: Preprocess CT volume
    fprintf('----------------------------------------\n');
    fprintf('STEP 2: Preprocessing CT volume...\n');
    fprintf('----------------------------------------\n');
    fprintf('  Applying anisotropic diffusion filter...\n');

    % Create masked volume for filtering
    VolSq = volumes.ct;
    VolSq(volumes.lungs == 0) = min(VolSq(:)) - 1;

    % Apply anisotropic diffusion
    anisodiff = anisodiff3D(VolSq, params.filtering.kappa, ...
        params.filtering.gamma, params.filtering.radius, ...
        params.filtering.iterations, params.voxel.dimensions);

    fprintf('  Filtering complete\n\n');

    %% Step 3: Preprocess lungs
    fprintf('----------------------------------------\n');
    fprintf('STEP 3: Preprocessing lungs...\n');
    fprintf('----------------------------------------\n');

    % Reorient volumes to standard orientation
    volumes.ct = imrotate3(volumes.ct, 90, [0 0 1], 'crop');
    volumes.airways = imrotate3(volumes.airways, 90, [0 0 1], 'crop');
    volumes.lungs = imrotate3(volumes.lungs, 90, [0 0 1], 'crop');
    volumes.vessels = imrotate3(volumes.vessels, 90, [0 0 1], 'crop');
    volumes.consolidation = imrotate3(volumes.consolidation, 90, [0 0 1], 'crop');
    volumes.injury = imrotate3(volumes.injury, 90, [0 0 1], 'crop');
    
    volumes.filtered = imrotate3(anisodiff, 90, [0 0 1], 'crop');

    % Preprocess lung masks
    [lungsProcessed, cardinalPoints] = preprocessLungs(volumes, params);

    %% Step 4: Refine lobe segmentation (if available)
    if volumes.hasLobes
        fprintf('----------------------------------------\n');
        fprintf('STEP 4: Refining lobe segmentation...\n');
        fprintf('----------------------------------------\n');
        volumes.lobes = imrotate3(volumes.lobes, 90, [0 0 1], 'crop');

        lobesRefined = refineLobeSegmentation(volumes.lobes, ...
            lungsProcessed.binary, params);
        results.lobes = lobesRefined;
    end

    %% Step 5: Classify distance from pleural surface
    fprintf('----------------------------------------\n');
    fprintf('STEP 5: Classifying parenchymal regions...\n');
    fprintf('----------------------------------------\n');
    [distanceLabelmap, distanceMasks] = classifyDistanceFromBorder(...
        lungsProcessed.convexhull, lungsProcessed.binary, params);

    %% Step 6: Segment and refine airways
    fprintf('----------------------------------------\n');
    fprintf('STEP 6: Segmenting and refining airways...\n');
    fprintf('----------------------------------------\n');
    [airwaysSeg, ~] = refineAirways(volumes, lungsProcessed, distanceMasks, params);
    results.airways = airwaysSeg;

    %% Step 7: Segment and refine injuries (low/high attenuation injuries)
    fprintf('----------------------------------------\n');
    fprintf('STEP 7: Segmenting pathological tissue...\n');
    fprintf('----------------------------------------\n');
    injurySeg = segmentInjuries(volumes,  lungsProcessed, airwaysSeg, params);
    results.injuries = injurySeg; 
    % REVISE BRIGHTNESS MAPPING APPROACH AND SEPARATION BETWEEN GGO AND CONSOLIDATION

    %% Step 8: Segment vessels with Jerman filter
    fprintf('----------------------------------------\n');
    fprintf('STEP 8: Segmenting vessels with Jerman filter...\n');
    fprintf('----------------------------------------\n');
    vesselsSeg = segmentVessels(volumes, lungsProcessed, airwaysSeg, ...
        injurySeg, distanceMasks, params);

    %% Step 9: Refine vessels
    fprintf('----------------------------------------\n');
    fprintf('STEP 9: Post-processing vessels...\n');
    fprintf('----------------------------------------\n');
    [vesselsRefined, largeVessels] = refineVessels(vesselsSeg, volumes, ...
        lungsProcessed, airwaysSeg, injurySeg, distanceMasks, params);
    results.vessels = vesselsRefined;
    results.largeVessels = largeVessels;

    %% Step 10: Classify vessels by diameter
    fprintf('----------------------------------------\n');
    fprintf('STEP 10: Classifying vessels by diameter...\n');
    fprintf('----------------------------------------\n');
    vesselsClassified = classifyVesselsByDiameter(vesselsRefined, params);
    results.vesselsclass = vesselsClassified;

     %% Step 11: Segment fissures
    fprintf('----------------------------------------\n');
    fprintf('STEP 11: Segmenting fissures...\n');
    fprintf('----------------------------------------\n');
    if volumes.hasLobes
        [fissures, fissureInfo] = segmentFissures(volumes, lungsProcessed, ...
            lobesRefined, params);
    else
        fissures = false(size(vesselsRefined));
        fissureInfo = struct();
    end    

    %% Step 12: Segment walls
    fprintf('----------------------------------------\n');
    fprintf('STEP 12: Segmenting tubular structure walls...\n');
    fprintf('----------------------------------------\n');
    airwayWalls = segmentWalls(airwaysSeg.airways, 'airway', params);
    vesselWalls = segmentWalls(vesselsRefined, 'vessel', params);
    largeVesselWalls = segmentWalls(largeVessels, 'large_vessel', params);

    % Trachea walls
    if isfield(airwaysSeg, 'trachea') && any(airwaysSeg.trachea(:))
        tracheaWalls = segmentWalls(airwaysSeg.trachea, 'airway', params);
    else
        tracheaWalls = false(size(vesselsRefined));
    end

    %% Step 13: Create regional maps
    fprintf('----------------------------------------\n');
    fprintf('STEP 13: Creating regional analysis maps...\n');
    fprintf('----------------------------------------\n');
    lungs_R = lungsProcessed.binary & (volumes.lungs == 1);
    lungs_L = lungsProcessed.binary & (volumes.lungs == 2);
    regionalMaps = createRegionalMaps(lungs_R, lungs_L, lungsProcessed.binary, cardinalPoints);

    %% Step 14: Create final label map
    fprintf('----------------------------------------\n');
    fprintf('STEP 14: Creating final label map...\n');
    fprintf('----------------------------------------\n');

    % Consolidation segmentation (dense opacities from injury segmentation)
    consolidationSeg = injurySeg.dense;
    consolidationSeg = consolidationSeg & ~vesselsRefined;
    consolidationSeg = consolidationSeg & ~airwaysSeg.airways;
    consolidationSeg = bwareaopen(consolidationSeg, 20, 6);

    % Air segmentation (use pre-computed from injurySeg)
    airSeg = injurySeg.lowAttenuation;
    airSeg = airSeg & ~airwaysSeg.airways;

    % Compromised tissue (GGO - ground glass opacity)
    compromisedSeg = injurySeg.ggo;
    compromisedSeg = compromisedSeg & ~consolidationSeg;
    compromisedSeg = compromisedSeg & ~vesselsRefined;
    compromisedSeg = compromisedSeg & ~airwaysSeg.airways;
    compromisedSeg = bwareaopen(compromisedSeg, 30, 6);

    % Create segmentation structure
    segmentations.lungs = lungsProcessed.binary;
    if volumes.hasLobes
        segmentations.lobes = lobesRefined;
    end
    segmentations.airways = airwaysSeg.airways;
    if isfield(airwaysSeg, 'trachea')
        segmentations.trachea = airwaysSeg.trachea;
    else
        segmentations.trachea = false(size(vesselsRefined));
    end
    segmentations.vessels = vesselsRefined;
    segmentations.largeVessels = largeVessels;
    segmentations.vesselsClassified = vesselsClassified;
    segmentations.consolidation = consolidationSeg;
    segmentations.compromised = compromisedSeg;
    segmentations.air = airSeg;
    segmentations.airwayWalls = airwayWalls;
    segmentations.vesselWalls = vesselWalls;
    segmentations.largeVesselWalls = largeVesselWalls;
    segmentations.tracheaWalls = tracheaWalls;
    segmentations.fissures = fissures;

    totalLabelMap = createLabelMap(segmentations, lungsProcessed.binary);

    %% Save results
    fprintf('----------------------------------------\n');
    fprintf('Saving results...\n');
    fprintf('----------------------------------------\n');

    % Build results structure for saveResults
    results.labelMap = totalLabelMap;
    results.vesselClassified = vesselsClassified;
    results.segmentsMapSI = regionalMaps.SI;
    results.segmentsMapPA = regionalMaps.PA;
    results.segmentsMapLR = regionalMaps.LR;
    results.distalClass = distanceLabelmap;
    results.ctProcessed = volumes.ct;
    results.lungs = volumes.lungs;
    results.lungsbin = lungsProcessed.binary;
    if volumes.hasLobes
        results.lobesVolume = lobesRefined;
    else
        results.lobesVolume = [];
    end

    saveResults(results, patientData, params);

    fprintf('\n Patient %s completed successfully!\n', patientData.name);

end