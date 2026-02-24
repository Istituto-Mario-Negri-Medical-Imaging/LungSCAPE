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
%   LungScape/
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

    %% Start total timer
    tTotal = tic;
    stepTimes = struct();

    %% Step 0: Initialize parallel pool (shared across all steps)
    params_preload = struct('parallel', struct('workers', 5, 'useParallel', true));
    initParallelPool(params_preload);

    %% Step 1: Load data
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 1: Loading data...\n');
    fprintf('----------------------------------------\n');
    volumes = loadPatientVolumes(patientData);
    params = getProcessingParameters(volumes.nii);
    stepTimes.step01 = toc(tStep);
    fprintf('  [Step 1 completed in %.1f s]\n\n', stepTimes.step01);

    %% Step 2: Preprocess CT volume
    tStep = tic;
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
    stepTimes.step02 = toc(tStep);
    fprintf('  [Step 2 completed in %.1f s]\n\n', stepTimes.step02);

    %% Step 3: Preprocess lungs
    tStep = tic;
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
    stepTimes.step03 = toc(tStep);
    fprintf('  [Step 3 completed in %.1f s]\n\n', stepTimes.step03);

    %% Step 4: Refine lobe segmentation (if available)
    tStep = tic;
    if volumes.hasLobes
        fprintf('----------------------------------------\n');
        fprintf('STEP 4: Refining lobe segmentation...\n');
        fprintf('----------------------------------------\n');
        volumes.lobes = imrotate3(volumes.lobes, 90, [0 0 1], 'crop');

        lobesRefined = refineLobeSegmentation(volumes.lobes, ...
            lungsProcessed.binary, params);
        results.lobes = lobesRefined;
    end
    stepTimes.step04 = toc(tStep);
    fprintf('  [Step 4 completed in %.1f s]\n\n', stepTimes.step04);

    %% Step 5: Segment fissures (needed for vessel refinement)
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 5: Segmenting fissures...\n');
    fprintf('----------------------------------------\n');
    if volumes.hasLobes
        [fissures, ~] = segmentFissures(volumes, lungsProcessed, ...
            lobesRefined, params);
    else
        fissures = false(size(lungsProcessed.binary));
    end
    stepTimes.step05 = toc(tStep);
    fprintf('  [Step 5 completed in %.1f s]\n\n', stepTimes.step05);

    %% Step 6: Classify distance from pleural surface
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 6: Classifying parenchymal regions...\n');
    fprintf('----------------------------------------\n');
    [distanceLabelmap, distanceMasks] = classifyDistanceFromBorder(...
        lungsProcessed.convexhull, lungsProcessed.binary, params);
    stepTimes.step06 = toc(tStep);
    fprintf('  [Step 6 completed in %.1f s]\n\n', stepTimes.step06);

    %% Step 7: Segment and refine airways
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 7: Segmenting and refining airways...\n');
    fprintf('----------------------------------------\n');
    [airwaysSeg, ~] = refineAirways(volumes, lungsProcessed, distanceMasks, params);
    results.airways = airwaysSeg;
    stepTimes.step07 = toc(tStep);
    fprintf('  [Step 7 completed in %.1f s]\n\n', stepTimes.step07);

    %% Step 8: Segment and refine injuries (low/high attenuation injuries)
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 8: Segmenting pathological tissue...\n');
    fprintf('----------------------------------------\n');
    injurySeg = segmentInjuries(volumes,  lungsProcessed, airwaysSeg, params);
    results.injuries = injurySeg;
    % REVISE BRIGHTNESS MAPPING APPROACH AND SEPARATION BETWEEN GGO AND CONSOLIDATION
    stepTimes.step08 = toc(tStep);
    fprintf('  [Step 8 completed in %.1f s]\n\n', stepTimes.step08);

    %% Step 9: Segment vessels with Jerman filter
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 9: Segmenting vessels with Jerman filter...\n');
    fprintf('----------------------------------------\n');
    vesselsSeg = segmentVessels(volumes, lungsProcessed, airwaysSeg, ...
        injurySeg, distanceMasks, params);
    stepTimes.step09 = toc(tStep);
    fprintf('  [Step 9 completed in %.1f s]\n\n', stepTimes.step09);

    %% Step 10: Refine vessels
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 10: Post-processing vessels...\n');
    fprintf('----------------------------------------\n');
    [vesselsRefined, largeVessels, addConsolidation] = refineVessels(vesselsSeg, volumes, ...
        lungsProcessed, airwaysSeg, injurySeg, distanceMasks, fissures, params);
    results.vessels = vesselsRefined;
    injurySeg.dense = injurySeg.dense | addConsolidation;
    stepTimes.step10 = toc(tStep);
    fprintf('  [Step 10 completed in %.1f s]\n\n', stepTimes.step10);

    %% Step 11: Segment trachea and finalize bronchi
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 11: Segmenting trachea and finalizing bronchi...\n');
    fprintf('----------------------------------------\n');
    [airwaysSeg, injurySeg.lowAttenuation] = segmentTrachea(volumes, ...
        lungsProcessed, airwaysSeg, injurySeg.lowAttenuation, params);
    results.airways = airwaysSeg;
    stepTimes.step11 = toc(tStep);
    fprintf('  [Step 11 completed in %.1f s]\n\n', stepTimes.step11);

    %% Step 12: Refine low-attenuation segmentation (AirSeg)
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 12: Refining low-attenuation segmentation...\n');
    fprintf('----------------------------------------\n');
    injurySeg.lowAttenuation = refineAirSeg(volumes, lungsProcessed, ...
        injurySeg.lowAttenuation, injurySeg, airwaysSeg.airways, ...
        vesselsRefined, params);
    stepTimes.step12 = toc(tStep);
    fprintf('  [Step 12 completed in %.1f s]\n\n', stepTimes.step12);

    %% Step 13: Finalize segmentations (walls + masking + fibrotic bands)
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 13: Finalizing segmentations...\n');
    fprintf('----------------------------------------\n');

    % Prepare inputs for finalization
    denseSeg = injurySeg.dense;
    denseSeg = denseSeg & ~vesselsRefined;
    denseSeg = denseSeg & ~airwaysSeg.airways;
    denseSeg = bwareaopen(denseSeg, 20, 6);

    airSeg = injurySeg.lowAttenuation;
    airSeg = airSeg & ~airwaysSeg.airways;

    ggoSeg = injurySeg.ggo;
    ggoSeg = ggoSeg & ~denseSeg;
    ggoSeg = ggoSeg & ~vesselsRefined;
    ggoSeg = ggoSeg & ~airwaysSeg.airways;
    ggoSeg = bwareaopen(ggoSeg, 30, 6);

    trachea = [];
    if isfield(airwaysSeg, 'trachea')
        trachea = airwaysSeg.trachea;
    end
    if isempty(trachea)
        trachea = false(size(lungsProcessed.binary));
    end

    segs = finalizeSegmentations(vesselsRefined, largeVessels, ...
        airwaysSeg.airways, trachea, denseSeg, ggoSeg, airSeg, ...
        lungsProcessed.binary, distanceMasks, fissures, volumes, params);
    stepTimes.step13 = toc(tStep);
    fprintf('  [Step 13 completed in %.1f s]\n\n', stepTimes.step13);

    %% Step 14: Classify vessels by diameter
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 14: Classifying vessels by diameter...\n');
    fprintf('----------------------------------------\n');
    vesselsClassified = classifyVesselsByDiameter(segs.vessels, params);
    results.vesselsclass = vesselsClassified;
    stepTimes.step14 = toc(tStep);
    fprintf('  [Step 14 completed in %.1f s]\n\n', stepTimes.step14);

    %% Step 15: Create regional maps
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 15: Creating regional analysis maps...\n');
    fprintf('----------------------------------------\n');
    regionalMaps = createRegionalMaps(lungsProcessed.right, lungsProcessed.left, lungsProcessed.binary, cardinalPoints);
    stepTimes.step15 = toc(tStep);
    fprintf('  [Step 15 completed in %.1f s]\n\n', stepTimes.step15);

    %% Step 16: Create final label map
    tStep = tic;
    fprintf('----------------------------------------\n');
    fprintf('STEP 16: Creating final label map...\n');
    fprintf('----------------------------------------\n');

    % Create segmentation structure from finalized results
    segmentations.lungs = lungsProcessed.binary;
    if volumes.hasLobes
        segmentations.lobes = lobesRefined;
    end
    segmentations.airways = segs.airways;
    segmentations.trachea = segs.trachea;
    segmentations.vessels = segs.vessels;
    segmentations.dense = segs.dense;
    segmentations.ggo = segs.ggo;
    segmentations.air = segs.airSeg;
    segmentations.airwayWalls = segs.airwayWalls;
    segmentations.vesselWalls = segs.vesselWalls;
    segmentations.tracheaWalls = segs.tracheaWalls;

    totalLabelMap = createLabelMap(segmentations, lungsProcessed.binary);
    stepTimes.step16 = toc(tStep);
    fprintf('  [Step 16 completed in %.1f s]\n\n', stepTimes.step16);

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

    %% Cleanup parallel pool
    shutdownParallelPool();

    %% Print timing summary
    totalTime = toc(tTotal);
    fprintf('\n========== TIMING SUMMARY ==========\n');
    stepNames = {'Load data', 'Preprocess CT', 'Preprocess lungs', ...
        'Refine lobes', 'Segment fissures', 'Distance classification', ...
        'Refine airways', 'Segment pathology', 'Segment vessels (Jerman)', ...
        'Refine vessels', 'Segment trachea', 'Refine air seg', ...
        'Classify vessel diameter', 'Finalize segmentations', ...
        'Regional maps', 'Create label map'};
    fields = fieldnames(stepTimes);
    for i = 1:length(fields)
        fprintf('  Step %2d: %6.1f s  (%s)\n', i, stepTimes.(fields{i}), stepNames{i});
    end
    fprintf('  -----------------------------------\n');
    fprintf('  TOTAL:   %6.1f s  (%.1f min)\n', totalTime, totalTime/60);
    fprintf('====================================\n');

    fprintf('\n Patient %s completed successfully!\n', patientData.name);

end

%% Parallel Pool Management

function initParallelPool(params)
%INITPARALLELPOOL Create parallel pool if needed (once per patient)
    if params.parallel.useParallel && params.parallel.workers > 0
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool(params.parallel.workers);
            fprintf('  Parallel pool created (%d workers)\n', params.parallel.workers);
        end
    end
end

function shutdownParallelPool()
%SHUTDOWNPARALLELPOOL Delete parallel pool at end of patient processing
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
        fprintf('  Parallel pool closed\n');
    end
end