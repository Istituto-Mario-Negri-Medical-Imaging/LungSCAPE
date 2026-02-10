function volumes = loadPatientVolumes(patientData)
%LOADPATIENTVOLUMES Load all NIfTI volumes for a patient
%
%   volumes = loadPatientVolumes(patientData)
%
%   Loads all required and optional NIfTI volumes for a patient, including
%   CT volume and various segmentation masks.
%
%   Input:
%       patientData - Structure from loadPatientData with fields:
%           .name    - Patient name
%           .files   - Structure with file paths
%
%   Output:
%       volumes - Structure with fields:
%           .nii           - NIfTI file to get info
%           .ct            - CT volume (int16, HU values)
%           .lungs         - Lung segmentation (uint8: 0=bg, 1=right, 2=left)
%           .lobes         - Lobe segmentation (uint8: 1-5 for lobes)
%           .airways       - Airway segmentation (uint8 binary)
%           .vessels       - Vessel segmentation (uint8 binary)
%           .consolidation - Consolidation/dense opacities mask (uint8 binary)
%           .injury        - high attenuation abnormalities (uint8 binary)
%           .hasLobes      - Boolean indicating if lobes are available
%
%   Example:
%       patients = loadPatientData('DATA/');
%       volumes = loadPatientVolumes(patients(1));
%       processedVolume = anisodiff3D(volumes.ct, ...);
%
%   See also: loadPatientData, load_untouch_nii

fprintf('Loading volumes for patient: %s\n', patientData.name);

%% Load CT Volume
fprintf('  Loading CT volume...\n');
ctNii = load_untouch_nii(patientData.files.ct);
volumes.nii = ctNii;
volumes.ct = int16(ctNii.img);

% Apply rescale intercept if present
if isfield(ctNii.hdr.dime, 'scl_inter')
    volumes.ct = volumes.ct + ctNii.hdr.dime.scl_inter;
end

% Clip minimum HU value
volumes.ct(volumes.ct < -1024) = -1024;

%% Load Lungs Segmentation (Required)
fprintf('  Loading lungs segmentation...\n');
lungsNii = load_untouch_nii(patientData.files.lungs);
volumes.lungs = uint8(lungsNii.img);

%% Load Airways (Required)
fprintf('  Loading airway segmentation...\n');
airwayNii = load_untouch_nii(patientData.files.airway);
volumes.airways = uint8(airwayNii.img);

%% Load Vessels (Required)
fprintf('  Loading vessels segmentation...\n');
vesselsNii = load_untouch_nii(patientData.files.vessels);
volumes.vessels = uint8(vesselsNii.img);

%% Load Lobes (Optional)
if isfield(patientData.files, 'lobes')
    fprintf('  Loading lobes segmentation...\n');
    lobesNii = load_untouch_nii(patientData.files.lobes);
    % Clean up small misclassified regions
    volumes.lobes = cleanLobeSegmentation(uint8(lobesNii.img));
    volumes.hasLobes = true;
else
    fprintf('  No lobes segmentation found\n');
    volumes.lobes = [];
    volumes.hasLobes = false;
end

%% Load full opacities
if isfield(patientData.files, 'consolidation')
    fprintf('  Loading consolidation mask...\n');
    consolidationNii = load_untouch_nii(patientData.files.consolidation);
    volumes.consolidation = uint8(consolidationNii.img);
else
    fprintf('  No consolidation mask found\n');
    volumes.consolidation = zeros(size(volumes.ct), 'uint8');
end

%% Load injury (High Attenuation abnormalities) 
if isfield(patientData.files, 'injury')
    fprintf('  Loading injury mask...\n');
    injuryNii = load_untouch_nii(patientData.files.injury);
    volumes.injury = uint8(injuryNii.img);
else
    fprintf('  No injury mask found\n');
    volumes.injury = zeros(size(volumes.ct), 'uint8');
end

fprintf('  Volumes loaded successfully\n\n');

end

%% Helper Function

function lobesCleaned = cleanLobeSegmentation(lobes)
%CLEANLOBESEGMENTATION Remove small misclassified particles in lobe segmentation

minLobeVolume = 5000;  % Minimum volume in voxels
connectivity = 26;

lobesCleaned = zeros(size(lobes), 'uint8');

for lobeLabel = 1:5
    lobeMask = (lobes == lobeLabel);
    lobeMask = bwareaopen(lobeMask, minLobeVolume, connectivity);
    lobesCleaned(lobeMask) = lobeLabel;
end

end
