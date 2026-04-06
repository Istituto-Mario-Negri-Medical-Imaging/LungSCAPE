function saveResults(results, patientData, params)
%SAVERESULTS Save segmentation results for a patient
%
%   saveResults(results, patientData, params)
%
%   Saves all segmentation results as NIfTI files (.nii.gz) sharing the same
%   spatial header as the input CT. Output volumes are un-rotated before
%   saving (inverting the 90-degree rotation applied in Step 2) so that they
%   align perfectly with the input CT in any NIfTI-compatible viewer.
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
%       patientData - Patient data structure from loadPatientData
%       params      - Processing parameters structure
%
%   Output Files:
%       <PatientFolder>/Results/
%           - TotalLabelMap.nii.gz      : Complete segmentation label map
%           - VesselCalibreMap.nii.gz   : Vessel diameter classification (1=large,2=medium,3=small)
%           - VesselsTypeMap.nii.gz     : Vessel type (1=artery, 2=vein, 3=undetermined)
%           - <Patient>.mat             : MATLAB workspace with all variables
%
%   Label Map Values:
%       0 - Background (outside lungs)
%       1 - Healthy parenchyma
%       2 - Ground glass opacity (GGO)
%       3 - Consolidation
%       4 - Air trapping/cysts
%       5 - Vessels (all sizes incl. large hilar; see VesselsTypeMap for artery/vein breakdown)
%       6 - Airways (bronchi)
%       7 - Trachea
%       8 - Airway walls
%       9 - Trachea walls
%
%   Example:
%       results = processPatient(volumes, params);
%       saveResults(results, patientData, params);
%
%   See also: loadPatientData, processPatient, load_untouch_nii, save_untouch_nii

fprintf('Saving results for patient: %s\n', patientData.name);

%% Create output directory
outputDir = fullfile(patientData.folder, 'Results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('  Created output directory: %s\n', outputDir);
end

%% Load CT NIfTI as spatial header template
% All outputs will share the same header as the input CT so that they are
% perfectly aligned with the original scan in any NIfTI viewer.
ctNii = load_untouch_nii(patientData.files.ct);

%% Undo internal 90-degree rotation before saving
%
% During processing (LungScape Step 2) all volumes are rotated 90 degrees
% around the z-axis (imrotate3, [0 0 1]). Before saving we invert this
% rotation with -90 degrees so the output data is back in the original
% NIfTI coordinate frame.
% 'nearest' interpolation preserves integer label values exactly.

labelMapOut    = imrotate3(results.labelMap,        -90, [0 0 1], 'nearest', 'crop');
vesselClassOut = imrotate3(results.vesselClassified, -90, [0 0 1], 'nearest', 'crop');
if isfield(results, 'vesselType')
    vesselTypeOut = imrotate3(results.vesselType,   -90, [0 0 1], 'nearest', 'crop');
end

%% Save Total Label Map
fprintf('  Saving total label map...\n');
labelMapFile = fullfile(outputDir, 'TotalLabelMap.nii.gz');
saveNifti(ctNii, uint8(labelMapOut), labelMapFile);

%% Save Vessel Calibre Map
fprintf('  Saving vessel calibre map...\n');
vesselMapFile = fullfile(outputDir, 'VesselCalibreMap.nii.gz');
saveNifti(ctNii, uint8(vesselClassOut), vesselMapFile);

%% Save Vessel Type Map
if isfield(results, 'vesselType')
    fprintf('  Saving vessel type map (artery/vein)...\n');
    vesselTypeFile = fullfile(outputDir, 'VesselsTypeMap.nii.gz');
    saveNifti(ctNii, uint8(vesselTypeOut), vesselTypeFile);
end

%% Save MATLAB Workspace
fprintf('  Saving MATLAB workspace...\n');
matFile = fullfile(outputDir, sprintf('%s.mat', patientData.name));

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

save(matFile, '-struct', 'saveVars', '-v7.3');

fprintf('  Results saved successfully to: %s\n\n', outputDir);

%% Display summary statistics
displaySummaryStatistics(results, params);

end


%% Helper: Save a volume as NIfTI reusing the CT header
function saveNifti(templateNii, imgData, filepath)
%SAVENIFTI Save a volume as NIfTI using a template header.
%
%   The template header is copied verbatim (preserving sform/qform,
%   origin, voxel spacing, etc.) and only the data-type fields and image
%   data are replaced. This guarantees spatial alignment with the input CT.

nii = templateNii;
nii.img = imgData;

% Update data-type fields for uint8 label maps
nii.hdr.dime.datatype = 2;    % DT_UINT8
nii.hdr.dime.bitpix   = 8;
nii.hdr.dime.scl_slope = 1;
nii.hdr.dime.scl_inter = 0;
nii.hdr.dime.cal_max  = double(max(imgData(:)));
nii.hdr.dime.cal_min  = 0;

save_untouch_nii(nii, filepath);
end


%% Helper: Display summary statistics
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
