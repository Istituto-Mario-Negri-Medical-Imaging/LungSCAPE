function saveResults(results, patientData, params)
%SAVERESULTS Save segmentation results for a patient
%
%   saveResults(results, patientData, params)
%
%   Saves all segmentation results including the total label map, regional
%   maps, vessel classification, and processed CT volume in NRRD and MAT formats.
%
%   Inputs:
%       results     - Structure with segmentation results:
%           .labelMap          - Total label map (uint8)
%           .vesselClassified  - Vessel diameter classification (uint8)
%           .vesselType        - Vessel type map (uint8): 1=artery,2=vein,3=undetermined
%           .segmentsMapSI     - Superior-Inferior regional map (uint8)
%           .segmentsMapPA     - Posterior-Anterior regional map (uint8)
%           .segmentsMapLR     - Left-Right regional map (uint8)
%           .lobesVolume       - Lobe segmentation (uint8)
%           .distalClass       - Distance classification (uint8)
%           .ctProcessed       - Filtered CT volume (int16)
%       patientData - Patient data structure from loadPatientData
%       params      - Processing parameters structure
%
%   Output Files:
%       <PatientFolder>/Results/
%           - TotalLabelMap.nrrd     : Complete segmentation label map
%           - VesselsLabelMap.nrrd   : Vessel diameter classification (1=large,2=medium,3=small)
%           - VesselsTypeMap.nrrd    : Vessel type (1=artery, 2=vein, 3=undetermined)
%           - Vol_<Patient>.nrrd     : Original/processed CT volume
%           - <Patient>.mat          : MATLAB workspace with all variables
%
%   Spatial alignment:
%       NRRD files are saved with the spatial origin from the input NIfTI header
%       (params.voxel.origin). The 90-degree rotation applied during processing
%       is inverted before saving so that NRRD outputs align correctly with the
%       original CT NIfTI when overlaid in 3D Slicer or other viewers.
%       The MAT workspace retains variables in the internal (rotated) coordinate
%       frame used during processing.
%
%   Label Map Values:
%       0 - Background (outside lungs)
%       1 - Healthy parenchyma
%       2 - Ground glass opacity (GGO)
%       3 - Consolidation
%       4 - Air trapping/cysts
%       5 - Vessels (all; see VesselsTypeMap for artery/vein breakdown)
%       6 - Airways (bronchi)
%       7 - Trachea
%       8 - Airway walls
%       9 - Trachea walls
%      10 - Vessel walls
%
%   Example:
%       results = processPatient(volumes, params);
%       saveResults(results, patientData, params);
%
%   See also: loadPatientData, processPatient, nrrdWriter

fprintf('Saving results for patient: %s\n', patientData.name);

%% Create output directory
outputDir = fullfile(patientData.folder, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('  Created output directory: %s\n', outputDir);
end

%% Undo internal 90-degree rotation before saving NRRD outputs
%
% During processing (LungScape Step 2) all volumes are rotated 90 degrees
% around the z-axis (imrotate3, [0 0 1]) to bring them into a standard
% processing orientation.  Before writing NRRD files we invert this rotation
% so that the output data is back in the original NIfTI coordinate frame.
% This allows the NRRD files to be correctly overlaid on the input CT NIfTI
% in 3D Slicer and other viewers using the spatial origin stored in the header.
%
% 'nearest' interpolation is used to preserve integer label values exactly.

labelMapOut       = imrotate3(results.labelMap,       -90, [0 0 1], 'nearest', 'crop');
vesselClassOut    = imrotate3(results.vesselClassified,-90, [0 0 1], 'nearest', 'crop');
ctOut             = imrotate3(results.ctProcessed,     -90, [0 0 1], 'nearest', 'crop');
if isfield(results, 'vesselType')
    vesselTypeOut = imrotate3(results.vesselType,      -90, [0 0 1], 'nearest', 'crop');
end

origin = params.voxel.origin;

%% Save Total Label Map (NRRD)
fprintf('  Saving total label map...\n');
labelMapFile = fullfile(outputDir, 'TotalLabelMap.nrrd');
nrrdWriter(labelMapFile, uint8(labelMapOut), ...
    params.voxel.dimensions, origin, 'raw');

%% Save Vessel Classification Label Map (NRRD)
fprintf('  Saving vessel classification label map...\n');
vesselMapFile = fullfile(outputDir, 'VesselsLabelMap.nrrd');
nrrdWriter(vesselMapFile, uint8(vesselClassOut), ...
    params.voxel.dimensions, origin, 'raw');

%% Save Vessel Type Map (NRRD)
if isfield(results, 'vesselType')
    fprintf('  Saving vessel type map (artery/vein)...\n');
    vesselTypeFile = fullfile(outputDir, 'VesselsTypeMap.nrrd');
    nrrdWriter(vesselTypeFile, uint8(vesselTypeOut), ...
        params.voxel.dimensions, origin, 'raw');
end

%% Save Processed CT Volume (NRRD)
fprintf('  Saving CT volume...\n');
ctFile = fullfile(outputDir, sprintf('Vol_%s.nrrd', patientData.name));
nrrdWriter(ctFile, int16(ctOut), ...
    params.voxel.dimensions, origin, 'raw');

%% Save MATLAB Workspace
fprintf('  Saving MATLAB workspace...\n');
matFile = fullfile(outputDir, sprintf('%s.mat', patientData.name));

% Prepare variables to save
saveVars = struct();
saveVars.TotalLabelMap = results.labelMap;
saveVars.Vessels_classified = results.vesselClassified;
if isfield(results, 'vesselType')
    saveVars.Vessels_type = results.vesselType;
end
saveVars.SegmentsMapSI = results.segmentsMapSI;
saveVars.SegmentsMapPA = results.segmentsMapPA;
saveVars.SegmentsMapLR = results.segmentsMapLR;
saveVars.distal_classification = results.distalClass;
saveVars.Vol = results.ctProcessed;
saveVars.lungs = results.lungs;
saveVars.lungsbin = results.lungsbin;
saveVars.VoxelVolume = params.voxel.volume;
saveVars.px = params.voxel.px;
saveVars.py = params.voxel.py;
saveVars.vz = params.voxel.vz;
saveVars.PatientName = patientData.name;
saveVars.ProcessingDate = datetime('now');
saveVars.Parameters = params;

if ~isempty(results.lobesVolume)
    saveVars.LobesVolume = results.lobesVolume;
end

% Save to MAT file
save(matFile, '-struct', 'saveVars', '-v7.3');

fprintf('  Results saved successfully to: %s\n\n', outputDir);

%% Display summary statistics
displaySummaryStatistics(results, params);

end

%% Helper Function

function displaySummaryStatistics(results, params)
%DISPLAYSUMMARYSTATISTICS Display summary of segmentation results

fprintf('=== Segmentation Summary ===\n');

labelMap    = results.labelMap;
voxelVolume = params.voxel.volume;  % mm³

% Volumes in mL for each tissue class
healthyVol       = sum(labelMap(:) == 1) * voxelVolume / 1000;
ggoVol           = sum(labelMap(:) == 2) * voxelVolume / 1000;
consolidationVol = sum(labelMap(:) == 3) * voxelVolume / 1000;
airVol           = sum(labelMap(:) == 4) * voxelVolume / 1000;
vesselVol        = sum(labelMap(:) == 5) * voxelVolume / 1000;
airwayVol        = sum(labelMap(:) == 6) * voxelVolume / 1000;
tracheaVol       = sum(labelMap(:) == 7) * voxelVolume / 1000;

totalLungVol = healthyVol + ggoVol + consolidationVol + airVol + ...
               vesselVol + airwayVol;

fprintf('  Total Lung Volume:    %8.1f mL\n', totalLungVol);
fprintf('  Healthy Parenchyma:   %8.1f mL (%5.1f%%)\n', ...
    healthyVol, 100*healthyVol/totalLungVol);
fprintf('  Ground Glass Opacity: %8.1f mL (%5.1f%%)\n', ...
    ggoVol, 100*ggoVol/totalLungVol);
fprintf('  Consolidation:        %8.1f mL (%5.1f%%)\n', ...
    consolidationVol, 100*consolidationVol/totalLungVol);
fprintf('  Air Trapping:         %8.1f mL (%5.1f%%)\n', ...
    airVol, 100*airVol/totalLungVol);
fprintf('  Vessels (total):      %8.1f mL (%5.1f%%)\n', ...
    vesselVol, 100*vesselVol/totalLungVol);

% Artery / vein breakdown from VesselsTypeMap (if available)
if isfield(results, 'vesselType') && ~isempty(results.vesselType)
    arteryVol = sum(results.vesselType(:) == 1) * voxelVolume / 1000;
    veinVol   = sum(results.vesselType(:) == 2) * voxelVolume / 1000;
    undetVol  = sum(results.vesselType(:) == 3) * voxelVolume / 1000;
    fprintf('    Arteries:           %8.1f mL (%5.1f%% of vessels)\n', ...
        arteryVol, 100*arteryVol/max(vesselVol,1e-6));
    fprintf('    Veins:              %8.1f mL (%5.1f%% of vessels)\n', ...
        veinVol,   100*veinVol/max(vesselVol,1e-6));
    if undetVol > 0
        fprintf('    Undetermined:       %8.1f mL (%5.1f%% of vessels)\n', ...
            undetVol, 100*undetVol/max(vesselVol,1e-6));
    end
end

fprintf('  Airways:              %8.1f mL (%5.1f%%)\n', ...
    airwayVol, 100*airwayVol/totalLungVol);
fprintf('  Trachea:              %8.1f mL\n', tracheaVol);
fprintf('============================\n\n');

end
