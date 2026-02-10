function params = getProcessingParameters(nii)
%GETPROCESSINGPARAMETERS Get processing parameters for lung segmentation
%
%   params = getProcessingParameters(nii)
%
%   Returns a structure containing all processing parameters for the lung
%   segmentation pipeline. Parameters are adjusted based on scanner
%   manufacturer and institution-specific protocols.
%
%   Inputs:
%       nii - NIfTI structure with header containing voxel dimensions
%
%   Output:
%       params - Structure with fields organized by processing stage:
%           .voxel         - Voxel size information
%           .filtering     - Anisotropic diffusion parameters
%           .lungs         - Lung preprocessing parameters
%           .thresholds    - Hounsfield Unit thresholds
%           .air           - Air segmentation parameters
%           .injury        - Injury segmentation parameters
%           .airways       - Airway segmentation parameters
%           .vessels       - Vessel segmentation parameters
%           .fissures      - Fissure segmentation parameters
%           .distance      - Distance classification parameters
%           .parallel      - Parallel processing parameters
%
%   Example:
%       params = getProcessingParameters(nii);
%       filteredVol = anisodiff3D(volume, params.filtering.kappa, ...
%           params.filtering.gamma, params.filtering.radius, ...
%           params.filtering.iterations, params.voxel.dimensions);
%
%   See also: anisodiff3D, vesselness3D

%% Voxel Information
% Extract the 3D voxel dimensions from the NIfTI header
voxelDimensions = nii.hdr.dime.pixdim(2:4);

params.voxel.px = voxelDimensions(1);
params.voxel.py = voxelDimensions(2);
params.voxel.vz = voxelDimensions(3);
params.voxel.dimensions = voxelDimensions;
params.voxel.volume = abs(voxelDimensions(1) * voxelDimensions(2) * voxelDimensions(3));

%% Filtering Parameters (Anisotropic Diffusion)
params.filtering.kappa = 1;           % Edge preservation parameter
params.filtering.gamma = 3/44;        % Time step
params.filtering.radius = 50;         % Filter radius
params.filtering.iterations = 1;      % Number of iterations

%% Lung Preprocessing Parameters
params.lungs.erosionSize = 2;         % Erosion kernel size

%% Hounsfield Unit Thresholds
params.thresholds.air = -980;                 % Air threshold for airways/cysts
params.thresholds.consolidationMin = -300;    % Minimum HU for consolidation
params.thresholds.ggoMin = -700;              % Minimum HU for GGO

%% Air Segmentation Parameters
params.air.minSize26 = 16;            % Minimum air region size (26-conn)
params.air.minSizeInjury26 = 25;      % Min size for injury analysis (26-conn)
params.air.minSizeInjury8 = 4;        % Min size for injury analysis (8-conn)

%% Injury Segmentation Parameters
params.injury.adaptiveThresholdSep = 150;     % Adaptive threshold to separate GGO from consolidation
params.injury.kernelSize = 21;                % Gaussian kernel size for adaptive threshold
params.injury.peakOffset = 40;                % Offset from histogram peak for threshold
params.injury.vesselPeakOffset = 100;         % Offset for vessel/pathology separation
params.injury.histoMinHU = -970;              % Min HU for histogram analysis
params.injury.histoMaxHU = -600;              % Max HU for histogram analysis
params.injury.minSize6 = 35;                  % Min injury size (6-conn)
params.injury.minSize4 = 16;                  % Min injury size (4-conn)

%% Airway Segmentation Parameters
params.airways.minVolume = 10;                        % Minimum object volume (voxels)
params.airways.SmallVolumes = 1000;                   % Small volume threshold
params.airways.IslandsVolumes = 500;                  % Island volume threshold
params.airways.honeycombing.minDiameter = 3.5;        % Min diameter for honeycomb detection
params.airways.honeycombing.largeDiameter = 5;        % Large diameter for cyst detection
params.airways.skelMinBranch = 3;                     % Minimum branch length for skeleton
params.airways.anisotropyThreshold = 0.90;            % FA threshold for tubular airways
params.airways.wallMinDistance = 0;                   % Min distance for wall detection
params.airways.wallMaxDistance = 2;                   % Max distance for wall detection

%% Vessel Segmentation Parameters
% Jerman vesselness filter
params.vessels.jerman.scales = 0.5:0.5:1.5;   % Scales for Jerman filter
params.vessels.jerman.tau = 0.5;              % Jerman tau parameter
params.vessels.jerman.useBright = true;       % Use bright vessels

% Vesselness thresholds
params.vessels.threshold.initial = 0.65;      % Initial vesselness threshold
params.vessels.threshold.complete = 0.70;     % Threshold for complete vessels in healthy tissue

% Anisotropy thresholds
params.vessels.anisotropy.check = 0.92;       % FA for orientation check
params.vessels.anisotropy.final = 0.95;       % Final FA threshold

% Minimum eigenvalue for eigenvector filtering
params.vessels.minEigenvalue = 100;

% Pathological airway (honeycombing) parameters
params.vessels.aircystsVolume = 150;          % Min volume for main pathological airways
params.vessels.aircystsDistance = 4;          % Max distance from main pathological airways

% Graph-based vessel connectivity
params.vessels.graph.MinBranchLength = 15;    % Minimum branch length for vessel graph
params.vessels.graph.minVolume = 15;          % Minimum volume for vessel graph branches
params.vessels.graph.mainVolume = 300;        % Minimum volume for main vascular branches
params.vessels.graph.keepVolume = 15;         % Minimum volume for final vascular skeleton
params.vessels.graph.diameter.small = 2;      % Small vessel diameter threshold
params.vessels.graph.diameter.main = 2.5;     % Main branch diameter

% Vessel cleaning
params.vessels.cleaning.injuryRegionVolume = 100;   % Remove small vessels in injury regions
params.vessels.cleaning.finalMinVolume = 30;        % Remove small fragments in final vessels

% Large vessel detection
params.vessels.largeVessels.MinHU = -600;     % Minimum HU for large vessel detection
params.vessels.largeVessels.MaxHU = 100;      % Maximum HU for large vessel detection
params.vessels.largeVessels.diameter = 8;     % Min diameter for large vessel detection
params.vessels.largeVessels.minLength = 5;    % Min length for large vessel detection
params.vessels.largeVessels.minVol = 30;      % Min volume for large vessel detection

% Wall segmentation
params.vessels.wallMinDistance = 0;           % Min distance for wall detection
params.vessels.wallMaxDistance = 2;           % Max distance for wall detection

%% Fissure Segmentation Parameters
params.fissures.scales = [0.5, 1.0];          % Scales for vesselness filter (thin sheets)
params.fissures.tau = 0.5;                    % Vesselness tau parameter
params.fissures.vesselnessThreshold = 0.3;    % Min vesselness for fissure candidates
params.fissures.vesselnessMaxThreshold = 0.8; % Max vesselness (exclude strong tubular)
params.fissures.minHU = -950;                 % Min HU for fissures (air-like)
params.fissures.maxHU = -700;                 % Max HU for fissures
params.fissures.maxDistanceToBoundary = 5;    % Max distance to lobe boundary
params.fissures.minSize = 5000;               % Min fissure size (from ORIGINAL.m line 393)
params.fissures.dilationSize = 2;             % Dilation disk size (from ORIGINAL.m line 392)

%% Distance Classification Parameters (Close/Mid/Far from mediastinum)
params.distance.closeThreshold = 0.5;         % Threshold for close region
params.distance.closeThresholdClass = 0.6;    % Classification threshold
params.distance.midfarThreshold = 0.85;       % Midfar region threshold
params.distance.farThreshold = 0.85;          % Far region threshold

%% Parallel Processing Parameters
params.parallel.workers = 3;                  % Number of parallel workers
params.parallel.useParallel = true;           % Enable parallel processing

end
