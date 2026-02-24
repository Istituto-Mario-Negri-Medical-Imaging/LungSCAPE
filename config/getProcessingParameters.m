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
%           .airSeg        - Low-attenuation refinement parameters
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
params.thresholds.healthyMaxHU = -600;        % Max HU for healthy tissue (above → pathological)
params.thresholds.ggoMaxHU = -200;            % Max HU for GGO (above → dense)

%% Air Segmentation Parameters
params.air.minSize26 = 16;            % Minimum air region size (26-conn)
params.air.minSizeInjury26 = 25;      % Min size for injury analysis (26-conn)
params.air.minSizeInjury8 = 4;        % Min size for injury analysis (8-conn)

%% Injury Segmentation Parameters
params.injury.brightStructuresThreshold = 200; % Adaptive threshold for initial bright structures detection
params.injury.adaptiveThresholdSep = 150;     % Adaptive threshold to separate GGO from consolidation
params.injury.kernelSize = 21;                % Gaussian kernel size for adaptive threshold
params.injury.peakOffset = 40;                % Offset from histogram peak for threshold
params.injury.vesselPeakOffset = 100;         % Offset for vessel/pathology separation
params.injury.histoMinHU = -970;              % Min HU for histogram analysis
params.injury.histoMaxHU = -600;              % Max HU for histogram analysis
params.injury.minSize6 = 35;                  % Min injury size (6-conn)
params.injury.minSize4 = 16;                  % Min injury size (4-conn)
params.injury.combingRecoMinVol = 500;        % Min volume for combing reconstruction seed (8-conn)

%% Airway Segmentation Parameters
params.airways.minVolume = 10;                        % Minimum object volume (voxels)
params.airways.SmallVolumes = 1000;                   % Small volume threshold
params.airways.IslandsVolumes = 500;                  % Island volume threshold
params.airways.honeycombing.minDiameter = 3.5;        % Min diameter for honeycomb detection
params.airways.honeycombing.largeDiameter = 5;        % Large diameter for cyst detection
params.airways.skelMinBranch = 3;                     % Minimum branch length for skeleton
params.airways.anisotropyThreshold = 0.90;            % FA threshold for tubular airways
params.airways.leakage.largeDiameter = 6;             % Large leakage diameter threshold
params.airways.leakage.dilationSize = 5;              % Dilation size for far region in leakage check
params.airways.wallMinDistance = 0;                   % Min distance for wall detection
params.airways.wallMaxDistance = 2;                   % Max distance for wall detection

%% Trachea Segmentation Parameters
params.trachea.bodyMaskHU = -191;                       % HU threshold for body mask extraction
params.trachea.intensityThreshold = -850;               % HU threshold for air in trachea
params.trachea.minRoundness = 0.65;                     % Min roundness (sphericity) for trachea seeds
params.trachea.minSeedArea = 300;                       % Min area for trachea seed objects
params.trachea.airwaysMaxHU = -600;                     % Max HU for airway voxels
params.trachea.minVolumeForSkel = 200;                  % Min volume before skeletonization
params.trachea.branchMinLength = 15;                    % Min branch length for carina detection
params.trachea.minFragmentVolume = 150;                 % Fragments smaller than this go to AirSeg
params.trachea.bronchi.vesselness.scale = 1;            % Vesselness scale for bronchi check
params.trachea.bronchi.vesselness.tau = 0.7;            % Vesselness tau parameter
params.trachea.bronchi.vesselness.threshold = 0.9;      % Vesselness threshold
params.trachea.bronchi.endpointDistance = 6;             % Max distance from airway endpoints
params.trachea.bronchi.FA = 0.95;                       % FA threshold for bronchi recovery
params.trachea.bronchi.minRecoverVolume = 30;           % Min volume for recovered bronchi fragments
params.trachea.bronchi.minVolume18 = 70;                % Min volume for cleanup (18-conn)
params.trachea.bronchi.minVolume6 = 300;                % Min volume for cleanup (6-conn)
params.trachea.bronchi.skelMinLength = 70;              % Min skeleton length for length-based check
params.trachea.bronchi.maxCheckVolume = 1000;           % Max volume for second FA pass
params.trachea.bronchi.minVesselnessVol = 20;           % Min vesselness candidate volume

%% Low-Attenuation Refinement Parameters (air trapping / emphysema)
params.airSeg.maxHU = -950;                              % Remove AirSeg voxels above this HU
params.airSeg.expandHU = -900;                           % HU threshold for AirSeg expansion via injury combing
params.airSeg.healthyMinVol = 5;                         % Min volume for healthy fragments (8-conn)
params.airSeg.dilationSize = 2;                          % Dilation disk size for healthy→AirSeg recovery

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

% Multi-threshold vessel refinement in pathological zones
params.vessels.pathZone.vesselness85 = 0.85;         % Vesselness threshold for distal consolidation and honeycombing (#0, #2)
params.vessels.pathZone.vesselness80 = 0.80;         % Vesselness threshold for general consolidation (#3)
params.vessels.pathZone.midfarDilateSize = 2;        % Dilation for midfar→far reconstruction
params.vessels.pathZone.farDilateSize = 5;           % Dilation for far_LungsCenter

% Vessel cleaning
params.vessels.cleaning.injuryRegionVolume = 100;   % Remove small vessels in injury regions
params.vessels.cleaning.finalMinVolume = 30;        % Remove small fragments in final vessels

% Large vessel detection
params.vessels.largeVessels.MinHU = -600;     % Minimum HU for large vessel detection
params.vessels.largeVessels.MaxHU = 100;      % Maximum HU for large vessel detection
params.vessels.largeVessels.diameter = 8;     % Min diameter for large vessel detection
params.vessels.largeVessels.minLength = 5;    % Min length for large vessel detection
params.vessels.largeVessels.minVol = 30;      % Min volume for large vessel detection

% Large consolidation reclassification
params.vessels.largeConsol.minHU = -150;              % HU threshold for dense consolidation
params.vessels.largeConsol.initMinVol = 10;            % Initial cleanup (8-conn)
params.vessels.largeConsol.mainMinVol = 300;            % Main large consolidation min vol (26-conn)
params.vessels.largeConsol.distalMinVol = 50;           % Distal consolidation min vol (26-conn)
params.vessels.largeConsol.vesselMinVol = 30;           % Vessel-consolidated min vol (8-conn)
params.vessels.largeConsol.vesselness.scales = 0.5:0.5:1.5;  % Vesselness scales for re-filtering
params.vessels.largeConsol.vesselness.tau = 0.6;        % Vesselness tau
params.vessels.largeConsol.vesselness.threshold = 0.87; % Vesselness threshold
params.vessels.largeConsol.postMinVol = 150;            % Small post-vesselness fragment threshold

% Pre-wall reassignment
params.vessels.preWall.consolMaxVol = 30;            % Max vol for small consolidation fragments (8-conn)
params.vessels.preWall.airwaysDist = 3;              % Max distance from airways for wall reassignment
params.vessels.preWall.airwaysCloseDist = 1;         % Close distance for HU-based airway walls
params.vessels.preWall.vesselsDist = 2;              % Max distance from vessels for wall reassignment
params.vessels.preWall.airwayWallHU = -900;          % HU threshold for close airway wall detection

% Last false positives removal
params.vessels.lastFP.pass1.maxVol = 500;            % Max vol for FA check
params.vessels.lastFP.pass1.minDiam = 3;             % Min diameter seed
params.vessels.lastFP.pass1.minArea = 4;             % Min area for seed
params.vessels.lastFP.pass1.FA = 0.95;               % FA threshold
params.vessels.lastFP.pass2.maxVol = 70;             % Max vol for volume+diameter check
params.vessels.lastFP.pass2.minDiam = 3;             % Min diameter seed
params.vessels.lastFP.distal.minDiam = 3;            % Min diameter for distal segregated
params.vessels.lastFP.distal.minArea = 3;            % Min area for distal seed

% Wall segmentation
params.vessels.wallMinDistance = 0;           % Min distance for wall detection
params.vessels.wallMaxDistance = 2;           % Max distance for wall detection

% Last corrections - progressive leakage detection
params.vessels.lastCorr.fissureMaxVol = 1000;         % Max vol for fissure cleaning
params.vessels.lastCorr.outBordersDilate = 9;         % Dilation disk for distal region
params.vessels.lastCorr.outBorders2Dilate = 12;       % Dilation disk for external boundary
params.vessels.lastCorr.skelMinBranch = 5;            % Min branch length for diameter skeleton
params.vessels.lastCorr.bsapxMinVol = 7000;           % Min volume for base/apex detection
% Pass 0: small disconnected objects in distal region
params.vessels.lastCorr.pass0.maxVol = 1000;          % Max object volume
params.vessels.lastCorr.pass0.minDiam = 3;            % Min diameter for leakage seed
params.vessels.lastCorr.pass0.minArea = 5;            % Min area for leakage seed
% Pass 1: very small disconnected objects in distal region
params.vessels.lastCorr.pass1.maxVol = 300;           % Max object volume
params.vessels.lastCorr.pass1.minDiam = 3;            % Min diameter for leakage seed
% Pass 2: leakage in consolidation's distal region
params.vessels.lastCorr.pass2.minDiam = 2.5;          % Min diameter for leakage seed
params.vessels.lastCorr.pass2.minArea = 4;            % Min area for leakage seed
% Pass 3: general distal leakage
params.vessels.lastCorr.pass3.minDiam = 4;            % Min diameter for leakage seed
params.vessels.lastCorr.pass3.minArea = 5;            % Min area for leakage seed
% Pass 4: diameter-based disconnected from inner region
params.vessels.lastCorr.pass4.minDiam = 3;            % Min diameter for leakage seed
% Eigenvector orientation filter
params.vessels.lastCorr.eigvec.maxVol = 3000;         % Max vol for eigenvector check
params.vessels.lastCorr.eigvec.FA = 0.87;             % FA threshold
% Seed-based cleanup
params.vessels.lastCorr.seed.maxVol = 1000;           % Max vol for seed cleanup
% Final diameter + FA check
params.vessels.lastCorr.final.maxVol = 150;           % Max vol for final check
params.vessels.lastCorr.final.minDiam = 3.5;          % Min diameter for final check
params.vessels.lastCorr.final.FA = 0.94;              % FA threshold for final check

%% Fissure Segmentation Parameters (morphological lobe boundary approach)
params.fissures.dilationSize = 2;             % Dilation disk size 
params.fissures.minSize = 5000;               % Min fissure size 

%% Finalization Parameters
% Airway walls
params.finalize.airwayWalls.largeAirwayDiam = 4;     % Min diameter for large airways
params.finalize.airwayWalls.largeWallDist = 2.5;     % Wall distance for large airways
params.finalize.airwayWalls.normalWallDist = 2;      % Wall distance for normal airways
params.finalize.airwayWalls.tracheaWallDist = 2.5;   % Wall distance for trachea
% Large vessel walls
params.finalize.vesselWalls.largeDist = 1.5;         % Wall distance for large vessels
% Fibrotic bands recovery
params.finalize.fibrotic.C = 200;                    % Adaptive threshold constant
params.finalize.fibrotic.kernelSize = 9;             % Gaussian kernel size
params.finalize.fibrotic.minVol = 500;               % Min volume for fibrotic bands (6-conn)

%% Distance Classification Parameters (Close/Mid/Far from mediastinum)
params.distance.closeThreshold = 0.5;         % Threshold for close region
params.distance.closeThresholdClass = 0.6;    % Classification threshold
params.distance.midfarThreshold = 0.85;       % Midfar region threshold
params.distance.farThreshold = 0.85;          % Far region threshold

%% Parallel Processing Parameters
params.parallel.workers = 3;                  % Number of parallel workers
params.parallel.useParallel = true;           % Enable parallel processing

end
