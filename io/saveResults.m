function saveResults(results, patientData, params)
%SAVERESULTS Save segmentation results for a patient
%
%   saveResults(results, patientData, params)
%
%   Saves all segmentation results including the total label map, regional
%   maps, vessel classification, and processed CT volume in NRRD and MAT formats.
%
%   Inputs:
%       results                - Structure with segmentation results:
%           .labelMap            - Total label map (uint8)
%           .vesselClassified    - Vessel diameter classification (uint8)
%           .segmentsMapSI       - Superior-Inferior regional map (uint8)
%           .segmentsMapPA       - Posterior-Anterior regional map (uint8)
%           .segmentsMapLR       - Left-Right regional map (uint8)
%           .lobesVolume         - Lobe segmentation (uint8)
%           .distalClass         - Distance classification (uint8)
%           .ctProcessed         - Filtered CT volume (int16)
%       patientData            - Patient data structure from loadPatientData
%       params                 - Processing parameters structure
%
%   Output Files:
%       <PatientFolder>/InterimResults/
%           - TotalLabelMap.nrrd     : Complete segmentation label map
%           - Vol_<Patient>.nrrd     : Original/processed CT volume
%           - <Patient>.mat          : MATLAB workspace with all variables
%
%   Label Map Values:
%       0 - Background (outside lungs)
%       1 - Healthy parenchyma
%       2 - Ground glass opacity (GGO)
%       3 - Reticulation/Consolidation
%       4 - Air trapping/cysts
%       5 - Vessels
%       6 - Airways (bronchi)
%       7 - Trachea
%       8 - Vessels walls
%       9 - Airway walls
%      10 - Trachea walls
%
%   Example:
%       results = processPatient(volumes, params);
%       saveResults(results, patientData, params);
%
%   See also: loadPatientData, processPatient, nrrdWriter

fprintf('Saving results for patient: %s\n', patientData.name);

%% Create output directory
outputDir = fullfile(patientData.folder, 'InterimResults');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('  Created output directory: %s\n', outputDir);
end

%% Save Total Label Map (NRRD)
fprintf('  Saving total label map...\n');
labelMapFile = fullfile(outputDir, 'TotalLabelMap.nrrd');
nrrdWriter(labelMapFile, uint8(results.labelMap), ...
    params.voxel.dimensions, [0, 0, 0], 'raw');

%% Save Processed CT Volume (NRRD)
fprintf('  Saving CT volume...\n');
ctFile = fullfile(outputDir, sprintf('Vol_%s.nrrd', patientData.name));
nrrdWriter(ctFile, int16(results.ctProcessed), ...
    params.voxel.dimensions, [0, 0, 0], 'raw');

%% Save MATLAB Workspace
fprintf('  Saving MATLAB workspace...\n');
matFile = fullfile(outputDir, sprintf('%s.mat', patientData.name));

% Prepare variables to save
saveVars = struct();
saveVars.TotalLabelMap = results.labelMap;
saveVars.Vessels_classified = results.vesselClassified;
saveVars.SegmentsMapSI = results.segmentsMapSI;
saveVars.SegmentsMapPA = results.segmentsMapPA;
saveVars.SegmentsMapLR = results.segmentsMapLR;
saveVars.distal_classification = results.distalClass;
saveVars.Vol = results.ctProcessed;
saveVars.lungs = results.lungs;
saveVars.lungsbin = results.lungsbin;
saveVars.VoxelVolume = params.voxel.volume;
saveVars.px = params.voxel.px;
saveVars.vx = params.voxel.vz;
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

% Count voxels for each label
labelMap = results.labelMap;
voxelVolume = params.voxel.volume;  % mm³

% Calculate volumes for each tissue type (in mL)
healthyVol = sum(labelMap(:) == 1) * voxelVolume / 1000;
ggoVol = sum(labelMap(:) == 2) * voxelVolume / 1000;
consolidationVol = sum(labelMap(:) == 3) * voxelVolume / 1000;
airVol = sum(labelMap(:) == 4) * voxelVolume / 1000;
vesselVol = sum(labelMap(:) == 5) * voxelVolume / 1000;
airwayVol = sum(labelMap(:) == 6) * voxelVolume / 1000;
tracheaVol = sum(labelMap(:) == 7) * voxelVolume / 1000;
largeVesselVol = sum(labelMap(:) == 8) * voxelVolume / 1000;

totalLungVol = healthyVol + ggoVol + consolidationVol + airVol + ...
               vesselVol + airwayVol + largeVesselVol;

fprintf('  Total Lung Volume:    %8.1f mL\n', totalLungVol);
fprintf('  Healthy Parenchyma:   %8.1f mL (%5.1f%%)\n', ...
    healthyVol, 100*healthyVol/totalLungVol);
fprintf('  Ground Glass Opacity: %8.1f mL (%5.1f%%)\n', ...
    ggoVol, 100*ggoVol/totalLungVol);
fprintf('  Consolidation:        %8.1f mL (%5.1f%%)\n', ...
    consolidationVol, 100*consolidationVol/totalLungVol);
fprintf('  Air Trapping:         %8.1f mL (%5.1f%%)\n', ...
    airVol, 100*airVol/totalLungVol);
fprintf('  Vessels (Small/Med):  %8.1f mL (%5.1f%%)\n', ...
    vesselVol, 100*vesselVol/totalLungVol);
fprintf('  Vessels (Large):      %8.1f mL (%5.1f%%)\n', ...
    largeVesselVol, 100*largeVesselVol/totalLungVol);
fprintf('  Airways:              %8.1f mL (%5.1f%%)\n', ...
    airwayVol, 100*airwayVol/totalLungVol);
fprintf('  Trachea:              %8.1f mL\n', tracheaVol);
fprintf('============================\n\n');

end
