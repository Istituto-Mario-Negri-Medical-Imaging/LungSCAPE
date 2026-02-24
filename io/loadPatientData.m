function patientList = loadPatientData(dataDirectory)
%LOADPATIENTDATA Load patient dataset information
%
%   patientList = loadPatientData(dataDirectory)
%
%   Scans the data directory and returns a structure array containing
%   information about each patient and their associated files.
%
%   Input:
%       dataDirectory - Path to directory containing patient folders
%
%   Output:
%       patientList - Structure array with fields:
%           .name        - Patient identifier/folder name
%           .folder      - Full path to patient folder
%           .files       - Structure of file paths (CT, masks, etc.)
%
%   Expected Directory Structure:
%       dataDirectory/
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
%   Example:
%       patients = loadPatientData('DATA');
%       for i = 1:length(patients)
%           fprintf('Processing: %s\n', patients(i).name);
%           volumes = loadPatientVolumes(patients(i));
%       end
%
%   See also: loadPatientVolumes, dir

% Validate input
if ~exist(dataDirectory, 'dir')
    error('loadPatientData:DirectoryNotFound', ...
          'Data directory not found: %s', dataDirectory);
end

% Get list of patient folders (exclude . and ..)
dirContents = dir(dataDirectory);
dirContents = dirContents(~ismember({dirContents.name}, {'.', '..'}));

% Filter to only include directories
patientFolders = dirContents([dirContents.isdir]);
numPatients = length(patientFolders);

if numPatients == 0
    warning('loadPatientData:NoPatients', 'No patient folders found in %s', dataDirectory);
    patientList = struct([]);
    return;
end

fprintf('Found %d patient folders\n', numPatients);

% Initialize patient list structure
patientList(numPatients) = struct(...
    'name', '', ...
    'folder', '', ...
    'files', struct());

% Process each patient
for i = 1:numPatients
    patientName = patientFolders(i).name;
    patientFolder = fullfile(patientFolders(i).folder, patientName);

    patientList(i).name = patientName;
    patientList(i).folder = patientFolder;

    % Find all files in patient folder
    files = dir(patientFolder);
    files = files(~ismember({files.name}, {'.', '..'}));

    % Locate specific files
    try
        patientList(i).files = findPatientFiles(files, patientName, patientFolder);
    catch ME
        warning('loadPatientData:FileError', ...
                'Error finding files for patient %s: %s', patientName, ME.message);
        continue;
    end

    % Note: DICOM header reading removed. If voxel size or DICOM
    % information is required, compute or provide it elsewhere.

    fprintf('  [%d/%d] Loaded: %s\n', i, numPatients, patientName);
end

end

%% Helper Functions

function fileStruct = findPatientFiles(files, patientName, patientFolder)
%FINDPATIENTFILES Locate all required files for a patient

% Initialize file structure
fileStruct = struct();

% Define file patterns to search for
patterns = struct(...
    'ct',            sprintf('%s.nii.gz', patientName), ...
    'lobes',         'lobesTS', ...
    'lungs',         'lungsTS', ...
    'vessels',       'vesselsTS', ...
    'airway',        'airwaysTS', ...
    'consolidation', 'consolidation', ...
    'injury',        'HighAttenuation');

% Find CT volume (required)
ctFiles = files(contains({files.name}, patterns.ct));
if isempty(ctFiles)
    error('CT volume not found: %s', patterns.ct);
end
fileStruct.ct = fullfile(patientFolder, ctFiles(1).name);

% Find lungs segmentation (required)
lungsFiles = files(contains({files.name}, patterns.lungs));
if isempty(lungsFiles)
    error('Lungs segmentation not found');
end
fileStruct.lungs = fullfile(patientFolder, lungsFiles(1).name);

% Find vessels segmentation (required)
vesselsFiles = files(contains({files.name}, patterns.vessels));
if isempty(vesselsFiles)
    error('Vessels segmentation not found');
end
fileStruct.vessels = fullfile(patientFolder, vesselsFiles(1).name);

% Find airway segmentation (required)
airwayFiles = files(contains({files.name}, patterns.airway));
if isempty(airwayFiles)
    error('Airway segmentation not found');
end
fileStruct.airway = fullfile(patientFolder, airwayFiles(1).name);

%%% Optional

% Find lobes segmentation
lobesFiles = files(contains({files.name}, patterns.lobes));
if ~isempty(lobesFiles)
    fileStruct.lobes = fullfile(patientFolder, lobesFiles(1).name);
end

% Find consolidation segmentation
consolidationFiles = files(contains({files.name}, patterns.consolidation));
if ~isempty(consolidationFiles)
    fileStruct.consolidation = fullfile(patientFolder, consolidationFiles(1).name);
end

% Find injury/high attenuation segmentation
injuryFiles = files(contains({files.name}, patterns.injury));
if ~isempty(injuryFiles)
    fileStruct.injury = fullfile(patientFolder, injuryFiles(1).name);
end

end

