function vesselsSeg = segmentVessels(volumes, lungsProcessed, airwaysSeg, injurySeg, distanceMasks, params)
%SEGMENTVESSELS Complete vessel segmentation pipeline with Jerman filter
%
%   vesselsSeg = segmentVessels(volumes, lungsProcessed, airwaysSeg, injurySeg, distanceMasks, params)
%
%   Segments pulmonary vessels using Jerman's vesselness filter (vesselness3D)
%   and sophisticated morphological post-processing. Handles both healthy and
%   pathological tissue regions differently.
%
%   Inputs:
%       volumes        - Structure with fields:
%           .ct          - Original CT volume (int16, HU values)
%           .filtered    - Anisotropic diffusion filtered volume
%           .vessels     - Combined artery+vein mask (arteries | veins)
%           .arteries    - TotalSegmentator artery mask (uint8, masked to lungs)
%           .veins       - TotalSegmentator vein mask   (uint8, masked to lungs)
%       lungsProcessed - Structure with fields:
%           .binary       - Eroded lung mask
%       airwaysSeg     - Structure with fields:
%           .airways     - Clean airway segmentation
%           .spurious    - Spurious airways
%       injurySeg      -  Structure with fields:
%           .dense       - Compromised tissue mask
%       distanceMasks  - Structure with distance-based lung masks:
%           .close       - Close to center mask
%           .midfar      - Middle distance mask
%           .far         - Far from center mask
%       params - Processing parameters
%
%   Outputs:
%       vesselsSeg - Structure with fields:
%           .final            - Final vessel segmentation (logical)
%           .main             - Main vasculature only (logical)
%           .vesselness       - Vesselness response map (double, 0-1)
%           .preliminary      - Initial vessel segmentation after reconstruction (logical)
%           .rawThreshold     - Initial vessel segmentation before reconstruction (logical)
%           .loopEndDiscarded - Vessels discarded by loop/end classification (logical)
%           .pathZoneDiscarded - Vessels discarded by pathological zone refinement (logical)
%
%   Algorithm Steps:
%       1. Apply Jerman vesselness filter to combined vessel mask (arteries | veins)
%       2. Threshold vesselness response
%       3. Reconstruct vessels from central lung regions
%       4. Recover bronchovascular bundles near airways (arterial territory only)
%       5. Remove loops and spurious extremities using graph analysis
%          - Generic (non-venous) voxels: FA-based shape filtering
%          - Venous voxels (near veins mask): permissive filter (FA bypassed)
%       6. Refine in pathological regions
%
%   See also: vesselness3D, refineVessels, classifyVesselsByDiameter, classifyVesselsByType

fprintf('  Segmenting pulmonary vessels...\n');

%% Apply Jerman vesselness filter on vessels preliminary segmentation
fprintf('    Applying Jerman vesselness filter...\n');

% Jerman filter with multiple scales to capture vessels of different sizes
% Scales: 0.5, 1.0, 1.5 mm capture small to medium vessels
% The last parameter (true) specifies the "Jerman" mode
vesselnessEnhanced = vesselness3D(volumes.vessels, ...
    params.vessels.jerman.scales, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    params.vessels.jerman.tau, ...
    params.vessels.jerman.useBright);  % true = Jerman mode

vesselsSeg.vesselness = vesselnessEnhanced;

%% Threshold vesselness response
fprintf('    Thresholding vesselness response...\n');

% Initial binary segmentation from vesselness
vesselsRawThreshold = vesselnessEnhanced > params.vessels.threshold.initial;

% Save the raw threshold before reconstruction (needed for post-trachea refinement)
vesselsSeg.rawThreshold = vesselsRawThreshold;

% Reconstruct vessels branching from the mediastinal region
% Dilate the "close" mask slightly to serve as seeds
seedMask = imdilate(distanceMasks.close, strel('sphere', 3));
vesselsInitialSeg = imreconstruct(seedMask, vesselsRawThreshold, 26);

% Store initial segmentation (after reconstruction)
vesselsSeg.preliminary = vesselsInitialSeg;

%% Extract skeleton and diameter information
fprintf('    Computing vessel diameter map...\n');

% Compute diameter map for the initial segmentation
[diameterImage, skeletonImage] = computeVesselDiameter(vesselsInitialSeg, 0);

% Identify endpoints and branchpoints
endpointsImage = bwmorph3(skeletonImage, 'endpoints');
branchpointsImage = bwmorph3(skeletonImage, 'branchpoints');

%% Segment main vasculature (diameter > 2.5mm)
fprintf('    Identifying main vascular branches...\n');

diameterImageMain = diameterImage > params.vessels.graph.diameter.main;
skeletonImageMain = bwskel(diameterImageMain, 'MinBranchLength', params.vessels.graph.MinBranchLength);
skeletonImageMain = bwareaopen(skeletonImageMain, params.vessels.graph.minVolume);

%% Graph-based connectivity reconstruction
fprintf('    Reconstructing vessel connectivity using graph analysis...\n');

% Connect fragmented main branches using shortest path analysis
connectionTracts = connectVesselBranches(skeletonImageMain, ...
    skeletonImage, ...
    params.vessels);

% Add connection tracts to main vasculature
skeletonImageMain = skeletonImageMain | connectionTracts;

% Reconstruct full vessels from main skeleton
skeletonImageMainDil = imdilate(skeletonImageMain, strel('sphere', 2));
mainVasculature = imreconstruct((distanceMasks.close | distanceMasks.midfar), ...
    vesselsInitialSeg | skeletonImageMainDil, 26);
mainVasculature = bwareaopen(mainVasculature, params.vessels.graph.mainVolume);

vesselsSeg.main = mainVasculature;

%% Remove spurious vessels and loops
fprintf('    Removing loops and spurious extremities...\n');

% Classify smaller vessels (diameter <= 2.5mm) into:
% - Extremities (ends): connected to endpoints
% - Loops: closed structures not connected to endpoints

diameterImageSmall = diameterImage <= params.vessels.graph.diameter.small & ...
    diameterImage > 0;

% Separate into loops and ends
vesselsLoops = diameterImageSmall & ...
    ~imreconstruct(endpointsImage, diameterImageSmall, 26);

vesselsEnds = imreconstruct(endpointsImage, diameterImageSmall, 26);

% Remove branchpoints from ends classification
vesselsEnds0 = vesselsEnds;
vesselsEnds0(branchpointsImage ~= 0) = 0;
vesselsEnds0 = imreconstruct(endpointsImage, vesselsEnds0, 26);

% --- Artery / vein pre-split ---
% Pulmonary veins run inter-segmentally and do not follow the bronchial
% tree. Their peripheral branches frequently appear as "loops" (both ends
% connect to existing vasculature, leaving no free endpoint) or as short
% "ends" not anchored to the main arterial tree. Applying the same FA-based
% filter to veins would cause systematic over-discarding of true venous
% inter-segmental bridges.
%
% Strategy: split loops and ends by overlap with the TS vein mask (dilated
% by 1 voxel to absorb boundary differences). Apply the standard
% filterSmallVessels to the generic (non-venous) subset and a more
% permissive filterVenousVessels to the venous subset.
veinsMaskDil = imdilate(volumes.veins, strel('sphere', params.vessels.graph.venous.maskDilSize));

loopsGeneric = vesselsLoops & ~veinsMaskDil;
loopsVenous  = vesselsLoops &  veinsMaskDil;
endsGeneric  = vesselsEnds0 & ~veinsMaskDil;
endsVenous   = vesselsEnds0 &  veinsMaskDil;

% Generic filter (FA-based, unchanged)
[loopsGenericKeep, loopsGenericDiscard] = filterSmallVessels(loopsGeneric, injurySeg.dense, params);
[endsGenericKeep,  endsGenericDiscard]  = filterSmallVessels(endsGeneric,  injurySeg.dense, params);

% Venous filter (permissive: no FA check, consolidation rule only)
[loopsVenousKeep, loopsVenousDiscard] = filterVenousVessels(loopsVenous, injurySeg.dense);
[endsVenousKeep,  endsVenousDiscard]  = filterVenousVessels(endsVenous,  injurySeg.dense);

% Extremities originating from main vasculature (always keep, valid for
% both arteries and veins that happen to be connected to the main tree)
vesselsEndsFromMain = imreconstruct(skeletonImageMainDil, vesselsEnds);

% Combine kept structures
vesselsToKeep = mainVasculature      | ...
                loopsGenericKeep     | ...
                loopsVenousKeep      | ...
                endsGenericKeep      | ...
                endsVenousKeep       | ...
                vesselsEndsFromMain;

vesselsToKeep = bwareaopen(vesselsToKeep, params.vessels.graph.keepVolume);

% Combine discarded structures (loop/end classification)
vesselsToDiscard = loopsGenericDiscard | loopsVenousDiscard | ...
                   endsGenericDiscard  | endsVenousDiscard;

%% Reconstruct full vessels from skeleton
fprintf('    Reconstructing full vessels...\n');

% Reconstruct vessels using distance-based labeling
vesselsReconstruction = single(vesselsInitialSeg);
[~, lbl_dist] = bwdist(skeletonImage);
lbl_dist = lbl_dist .* uint32(vesselsReconstruction);

reco_indexes = find(vesselsToKeep);
vesselsUpdate = ismember(lbl_dist, reco_indexes);

% Reconstruct discarded vessels (for downstream use in pathological zone and complete vessels)
discard_indexes = find(vesselsToDiscard);
loopEndDiscarded = ismember(lbl_dist, discard_indexes);
vesselsSeg.loopEndDiscarded = loopEndDiscarded;

%% Refine vessels near airways
fprintf('    Refining vessels near airways...\n');

% Recover missed vessels along airways (bronchovascular bundles).
% Constrained to arterial territory: pulmonary arteries travel alongside
% bronchi (bronchovascular bundle), while veins run inter-segmentally and
% are not spatially coupled to airways. Applying this recovery to veins
% would risk pulling in fibrotic bands or consolidated tissue near bronchi.
% A 1-voxel dilation of the TS artery mask accounts for minor boundary
% differences between the deep-learning segmentation and the Jerman output.
arteriesMask = imdilate(volumes.arteries, strel('sphere', 1));

vesselsNearAirways = recoverVesselsNearAirways(vesselsUpdate, ...
    airwaysSeg.airways, ...
    diameterImage, ...
    params);
vesselsNearAirways = vesselsNearAirways & arteriesMask;

vesselsUpdate = vesselsUpdate | vesselsNearAirways;

%% Remove vessels near pathological air cysts (honeycombing)
fprintf('    Removing fibrotic bands near honeycombing...\n');

% Vessels near pathological airways are likely fibrotic bands
if isfield(airwaysSeg, 'spurious')
    pathologicalAirwaysMain = bwareaopen(airwaysSeg.spurious, params.vessels.aircystsVolume);
    [pathAirwaysDist, ~] = bwdist(single(pathologicalAirwaysMain));
    pathAirwaysDist(skeletonImage == 0) = 0;
    vesselsNearCysts = pathAirwaysDist < params.vessels.aircystsDistance & pathAirwaysDist > 0;

    % Reconstruct full vessels near cysts
    vesselsNearCystsReco = reconstructBySeeds(vesselsNearCysts, ...
        vesselsInitialSeg);

    vesselsUpdate = vesselsUpdate & ~vesselsNearCystsReco;
end

%% Remove reticulations in distal consolidations
fprintf('    Removing reticulations in consolidated regions...\n');

vesselsUpdate = removeDistalReticulations(vesselsUpdate, ...
    injurySeg.dense, ...
    distanceMasks.far, ...
    params);

%% Multi-threshold vessel refinement in pathological zones
fprintf('    Multi-threshold vessel refinement in pathological zones...\n');

[vesselsUpdate, pathZoneDiscarded] = refineVesselsInPathologicalZones(...
    vesselsUpdate, vesselnessEnhanced, lungsProcessed.binary, ...
    injurySeg, airwaysSeg, distanceMasks, loopEndDiscarded, params);

vesselsSeg.pathZoneDiscarded = pathZoneDiscarded;

%% Final cleanup
fprintf('    Performing final cleanup...\n');

% Remove small objects in pathological regions
vesselsClean = bwareaopen(vesselsUpdate, params.vessels.cleaning.injuryRegionVolume);
vesselsCheck = vesselsUpdate & ~vesselsClean;

if isfield(injurySeg, 'highAttenuation')
    vesselsDiscard = imreconstruct(imdilate(injurySeg.highAttenuation, ...
        strel('sphere', 1)), vesselsCheck, 26);
    vesselsUpdate = vesselsUpdate & ~vesselsDiscard;
end

% Generic removal of small objects if not highly anisotropic
vesselsClean = bwareaopen(vesselsUpdate, 200);
vesselsCheck = vesselsUpdate & ~vesselsClean;
[vesselsCheckKeep, ~] = filterByAnisotropy(vesselsCheck, ...
    params.vessels.anisotropy.final);

vesselsUpdate = vesselsClean | vesselsCheckKeep;
vesselsUpdate = bwareaopen(vesselsUpdate, params.vessels.cleaning.finalMinVolume, 6);

%% Output
vesselsSeg.final = vesselsUpdate;
vesselsSeg.main = mainVasculature;

fprintf('    Vessel segmentation complete!\n');

end

%% Helper function: Filter small vessels by anisotropy (generic / arterial)
function [vesselsKeep, vesselsDiscard] = filterSmallVessels(vesselsMask, consolidationMask, params)
    % Filter small vessels based on size, location, and anisotropy.
    % Used for generic (non-venous) loops and ends.

    % Very small structures (< 7 voxels) - keep
    vesselsSmall = vesselsMask & ~bwareaopen(vesselsMask, 7);

    % Large structures in compromised tissue (>= 15 voxels) - discard
    vesselsLarge = bwareaopen(vesselsMask, 15);
    vesselsLarge(consolidationMask == 0) = 0;

    % Medium-sized - filter by anisotropy
    vesselsMedium = vesselsMask & ~(vesselsSmall | vesselsLarge);
    [vesselsMediumKeep, vesselsMediumDiscard] = filterByAnisotropy(vesselsMedium, ...
        params.vessels.anisotropy.final);

    vesselsKeep    = vesselsSmall | vesselsMediumKeep;
    vesselsDiscard = vesselsLarge | vesselsMediumDiscard;
end


%% Helper function: Filter venous loops/ends (permissive — no FA check)
function [vesselsKeep, vesselsDiscard] = filterVenousVessels(vesselsMask, consolidationMask)
    % Permissive filter for structures overlapping the TS vein mask.
    %
    % Pulmonary veins run inter-segmentally: their peripheral branches can
    % appear as topological loops (no free endpoint) or short ends not
    % connected to the arterial main tree. Both patterns produce low FA
    % at curved or branching segments, causing systematic over-discarding
    % when the standard FA-based filter is applied.
    %
    % Rules:
    %   < 7 voxels              → Keep  (same as generic filter)
    %   >= 15 voxels in consol. → Discard (consolidation false-positive
    %                             risk is unchanged regardless of vessel type)
    %   everything else         → Keep  (FA check bypassed: overlap with the
    %                             TS vein mask is sufficient evidence of
    %                             true vascular origin)

    % Very small structures (< 7 voxels) - keep
    vesselsSmall = vesselsMask & ~bwareaopen(vesselsMask, 7);

    % Large structures in compromised tissue (>= 15 voxels) - discard
    vesselsLarge = bwareaopen(vesselsMask, 15);
    vesselsLarge(consolidationMask == 0) = 0;

    % Remaining structures: keep without FA filtering
    vesselsRest = vesselsMask & ~(vesselsSmall | vesselsLarge);

    vesselsKeep    = vesselsSmall | vesselsRest;
    vesselsDiscard = vesselsLarge;
end


%% Helper function: Recover vessels near airways
function vesselsRecovered = recoverVesselsNearAirways(vessels, airways, diameterImage, params)
    % Recover vessels that run alongside airways (bronchovascular bundles)

    % Compute distance to airways
    [airwaysDist, ~] = bwdist(single(airways));

    % Only consider large vessels (diameter > 2mm)
    diameterImageMain = diameterImage > 2;
    airwaysDist(diameterImageMain == 0) = 0;

    % Vessels within 6mm of airways
    vesselsCloseToAirways = airwaysDist < 6 & airwaysDist > 0;

    % Reconstruct full vessels
    vesselsRecovered = imreconstruct(vesselsCloseToAirways, vessels, 8);
end


%% Helper function: Multi-threshold vessel refinement in pathological zones
function [vesselsUpdate, pathZoneDiscarded] = refineVesselsInPathologicalZones(...
    vesselsUpdate, vesselnessEnhanced, lungsbin, injurySeg, airwaysSeg, ...
    distanceMasks, loopEndDiscarded, params)
    % Applies increasing vesselness requirements in pathological zones
    % to remove false positive vessels (reticulations, consolidation leakage).
    % Implements 4 sub-steps with different vesselness thresholds:
    %   #0: Distal consolidation region → require vesselness > 0.85
    %   #1: Outer borders + bright structures in far region
    %   #2: Near honeycombing → require vesselness > 0.85
    %   #3: General consolidation → require vesselness > 0.80
    % Also includes loop/end discarded vessels from the graph analysis step.

    pz = params.vessels.pathZone;

    % Compute external boundary of lungs
    outBorders2 = bwmorph3(imfill(lungsbin, 8, 'holes'), 'remove');
    outBorders = distanceMasks.far & outBorders2;

    % Compute distal region excluding apex/base
    far_LungsCenter = imdilate(distanceMasks.midfar, strel('square', pz.midfarDilateSize));
    far_LungsCenter = imreconstruct(far_LungsCenter, distanceMasks.far, 8);
    far_LungsCenter = imdilate(far_LungsCenter, strel('disk', pz.farDilateSize));
    outBorders2far = outBorders2 & far_LungsCenter;

    consolidation = logical(injurySeg.dense);
    brightstructures = logical(injurySeg.brightstructures);

    % #0: Distal consolidation, require vesselness > 0.85
    discard0 = imreconstruct(uint8(outBorders2far), uint8(consolidation), 8);
    discard0(vesselnessEnhanced > pz.vesselness85) = 0;

    % #1: Outer borders + bright structures in far region
    discard1 = imreconstruct(uint8(outBorders2far), uint8(brightstructures), 8);
    discard1 = uint8(outBorders) | discard1;
    discard1(distanceMasks.far == 0) = 0;

    % #2: Near honeycombing, require vesselness > 0.85
    if isfield(airwaysSeg, 'spurious') && any(airwaysSeg.spurious(:))
        honeycombDil = imdilate(airwaysSeg.spurious, strel('disk', 1));
        discard2 = imreconstruct(uint8(honeycombDil), uint8(consolidation), 8);
        discard2(vesselnessEnhanced > pz.vesselness85) = 0;
    else
        discard2 = false(size(vesselsUpdate));
    end

    % #3: General consolidation, require vesselness > 0.80
    discard3 = consolidation;
    discard3(vesselnessEnhanced > pz.vesselness80) = 0;

    % Combine all discard zones (including loop/end discards from graph analysis)
    allDiscard = logical(discard3) | logical(discard0) | logical(discard1) | ...
        logical(discard2) | loopEndDiscarded;

    % Track what was removed for downstream reclassification
    pathZoneDiscarded = vesselsUpdate & allDiscard;

    % Update vessels
    vesselsUpdate = vesselsUpdate & ~allDiscard;
end


%% Helper function: Remove distal reticulations
function vesselsUpdate = removeDistalReticulations(vessels, consolidation, farMask, params)
    % Remove reticulations in distal consolidated regions using
    % eigenvalue/eigenvector analysis

    % Focus on distal consolidations
    vesselsCheck = vessels;
    vesselsCheck(consolidation == 0) = 0;
    vesselsCheck(farMask == 0) = 0;

    % Compute diameter
    [diameterImage, ~] = computeVesselDiameter(vesselsCheck, 5);

    % Only check larger structures (diameter > 4mm)
    diameterImageMain = diameterImage > 4;
    vesselsCheck = imreconstruct(diameterImageMain, vesselsCheck, 26);

    % Filter by eigenvector orientation (longitudinal = likely leakage)
    [vesselsDiscard, ~] = filterByEigenVectorOrientation(vesselsCheck, ...
        params.vessels.anisotropy.check, ...
        'longitudinal', ...
        params.vessels.minEigenvalue);

    % Expand to include nearby consolidation
    vesselsDiscard = imreconstruct(vesselsDiscard, logical(consolidation), 8);

    vesselsUpdate = vessels & ~vesselsDiscard;
end