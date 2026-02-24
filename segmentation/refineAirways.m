function [airwaysSegmented, airwaysWip] = refineAirways(volumes, lungsProcessed, distanceMasks, params)
%SEGMENTAIRWAYS Segment and refine airway tree
%
%   [airwaysSegmented, airwaysRemoved, airwaysWip] = segmentAirways(volumes, binary, params)
%
%   Segments airways and removes incorrectly segmented regions (honeycombing,
%   cysts, consolidation leakages) using diameter analysis and anisotropy.
%
%   Inputs:
%       volumes        - Structure with CT, initial airway/consolidation segmentation
%       lungsProcessed - Structure with the eroded lung mask (binary)
%       distanceMasks - Lung portion far from center
%       params         - Processing parameters
%
%   Output:
%       airwaysSegmented - Structure with fields:
%           .airways           - Clean airway segmentation
%           .trachea           - Trachea segmentation
%           .walls             - Airway wall segmentation
%           .spurious          - Spurious airways (honeycombing)
%       airwaysWip   - Structure with fields:
%           .airways           - Updated airway segmentation
%           .trachea           - Updated trachea segmentation
%   See also: refineAirways, computeVesselDiameter

fprintf('Segmenting airways...\n');

%% Initial cleanup
fprintf('  Removing spurious fragments - Step0 ...\n');
Airway = volumes.airways;
Airway(lungsProcessed.binary == 0) = 0;
Airway = bwareaopen(Airway, params.airways.minVolume);
airwaysWip.step0 = Airway;
AirwayInit = Airway;  % Preserve for diameter-based leakage step 

%% Detect honeycombing and cysts
fprintf('  Removing honeycombing (distal consolidation-connected airspaces)- Step1 ...\n');
% Step1 - Identify large airways (diameter > 3.5mm)
[diameterImage, skeletonImage] = computeVesselDiameter(Airway, ...
    params.airways.skelMinBranch);

diameterImage_airways = diameterImage > params.airways.honeycombing.minDiameter;

% Step1.2 - Identify distal consolidations and reconstruct connected large airspaces
honeycombing_walls = volumes.consolidation & distanceMasks.far;
honeycombing_walls = imfill(honeycombing_walls,8,'holes');
honeycombing_space = imreconstruct(honeycombing_walls, Airway, 8); 
honeycombing_space = imreconstruct(diameterImage_airways, honeycombing_space, 8);

% Step1.3 - % Remove small objects containing honeycombing markers from Step1.2
Airways_Small = Airway & ~bwareaopen(Airway,params.airways.SmallVolumes,26);
Airways_Remove = imreconstruct(honeycombing_space, Airways_Small, 26);
Airways_Small = Airways_Small & ~Airways_Remove;

% Update Airway
Airway = bwareaopen(Airway,params.airways.SmallVolumes,26) | Airways_Small;
honeycombing_space = honeycombing_space | Airways_Remove;
airwaysWip.step1 = Airway; 

% Step2
fprintf('  Removing cysts and enlargements - Step2 ...\n');
% Identify very large putative airways in distal regions
diameterImage_cysts = diameterImage > params.airways.honeycombing.largeDiameter;
diameterImage_cysts(distanceMasks.far == 0) = 0;

% Reconstruct full cysts from diameter map using distance labeling
Airways_distance = Airway & ~skeletonImage;
[~, Airways_distance_idx] = bwdist(skeletonImage);
Airways_distance_idx(~Airways_distance) = 0;

honeycombing_space_cysts = ismember(Airways_distance_idx, find(diameterImage_cysts));
honeycombing_space_cysts = honeycombing_space_cysts | diameterImage_cysts;

% Combine pathological airways
honeycombing_space = honeycombing_space | honeycombing_space_cysts;

% Update Airway
Airway = Airway & ~honeycombing_space;
airwaysWip.step2 = Airway; 

% Step3
fprintf('  Removing islands near cleaning interventions - Step3 ...\n');
% Check small airways near pathological regions
cleaningFragments = Airway & ~bwareaopen(Airway, params.airways.IslandsVolumes);

cleaningFragments = imreconstruct(imdilate(honeycombing_space, strel('sphere',1)), cleaningFragments, 26);

% Update Airway
Airway = Airway & ~cleaningFragments;
airwaysWip.step3 = Airway; 

honeycombing_space = honeycombing_space | cleaningFragments;

% Step3.5 - Diameter-based leakage removal
% Volume-independent rule: remove airways with diameter > 6 in far region
fprintf('  Removing diameter-based leakages (d>6) - Step3.5 ...\n');
air_leakage = diameterImage > params.airways.leakage.largeDiameter;
air_leakage(distanceMasks.far == 0) = 0;
AirwayCorr = AirwayInit;
AirwayCorr(imdilate(distanceMasks.far, strel('sphere', params.airways.leakage.dilationSize)) == 0) = 0;
air_leakage = imreconstruct(air_leakage, AirwayCorr, 26);
honeycombing_space = honeycombing_space | air_leakage;
Airway = Airway & ~air_leakage;
airwaysWip.step3_5 = Airway;

fprintf('  Diameter/Volume based cleaning complete\n');

% Step4
%% Filter by anisotropy (keep small structures with high FA)
fprintf('    Filtering islands by anisotropy (FA)...\n');
Airways_SmallTubular = Airway & ~bwareaopen(Airway, params.airways.IslandsVolumes, 26);
[Airways_tubular, Airways_blob] = filterByAnisotropy(Airways_SmallTubular, ...
    params.airways.anisotropyThreshold);

% Update Airway
Airway = bwareaopen(Airway, params.airways.IslandsVolumes, 26) | Airways_tubular;
airwaysWip.step4 = Airway;

honeycombing_space = honeycombing_space | Airways_blob;

%% Output
airwaysSegmented.airways = Airway;
airwaysSegmented.spurious = honeycombing_space;
airwaysSegmented.trachea = [];  % Populated by segmentTrachea
airwaysSegmented.walls = [];    % Populated by segmentWalls

fprintf('  Airway segmentation complete\n\n');

end
