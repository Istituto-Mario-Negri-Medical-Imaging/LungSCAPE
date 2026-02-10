function [fissures, fissureInfo] = segmentFissures(volumes, masks, lobes, params)
%SEGMENTFISSURES Segment lung fissures from CT volumes
%
%   [fissures, fissureInfo] = segmentFissures(volumes, masks, lobes, params)
%
%   Segments interlobar fissures by:
%   1. Extracting interfaces between lobes
%   2. Applying vesselness filter to enhance fissure-like structures
%   3. Combining with intensity-based thresholding
%   4. Morphological refinement
%
%   Fissures appear as thin, sheet-like structures separating lung lobes.
%   They are typically darker (lower HU) than surrounding parenchyma.
%
%   Inputs:
%       volumes - Structure with CT volumes
%       masks   - Structure with segmentation masks
%       lobes   - Lobe segmentation (uint8: 1-5)
%       params  - Processing parameters
%
%   Outputs:
%       fissures    - Binary fissure segmentation (logical)
%       fissureInfo - Structure with additional information:
%           .vesselness - Vesselness response emphasizing sheets
%           .boundaries - Lobe boundaries
%
%   See also: vesselness3D, segmentLobes

fprintf('  Segmenting lung fissures...\n');

% Check if lobes are available
if isempty(lobes) || all(lobes(:) == 0)
    warning('No lobe segmentation available. Skipping fissure segmentation.');
    fissures = false(size(masks.lungs));
    fissureInfo.vesselness = zeros(size(masks.lungs));
    fissureInfo.boundaries = false(size(masks.lungs));
    return;
end

%% Extract lobe boundaries
fprintf('    Extracting lobe boundaries...\n');

% Dilate each lobe and find overlaps (these are boundaries)
lobeBoundaries = false(size(lobes));

uniqueLobes = unique(lobes(:));
uniqueLobes(uniqueLobes == 0) = [];  % Remove background

for i = 1:length(uniqueLobes)
    lobeMask = lobes == uniqueLobes(i);
    lobeDilated = imdilate(lobeMask, strel('sphere', 2));

    % Overlap with other lobes defines boundaries
    for j = (i+1):length(uniqueLobes)
        otherLobeMask = lobes == uniqueLobes(j);
        overlap = lobeDilated & otherLobeMask;
        lobeBoundaries = lobeBoundaries | overlap;
    end
end

fissureInfo.boundaries = lobeBoundaries;

%% Apply vesselness filter for sheet-like structures
fprintf('    Applying vesselness filter for fissures...\n');

% Fissures are sheet-like (plate-like) structures
% Use vesselness with appropriate scale to detect them
% Smaller scale emphasizes thin fissure planes
fissureVesselness = vesselness3D(volumes.filtered, ...
    params.fissures.scales, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    params.fissures.tau, ...
    true);

fissureInfo.vesselness = fissureVesselness;

%% Threshold based on vesselness and intensity
fprintf('    Thresholding fissure candidates...\n');

% Fissures have moderate vesselness response
fissuresFromVesselness = fissureVesselness > params.fissures.vesselnessThreshold & ...
    fissureVesselness < params.fissures.vesselnessMaxThreshold;

% Fissures have relatively low HU values (darker than parenchyma)
fissuresFromIntensity = volumes.ct > params.fissures.minHU & ...
    volumes.ct < params.fissures.maxHU;

% Combine criteria
fissureCandidates = fissuresFromVesselness & fissuresFromIntensity;

% Must be near lobe boundaries
[boundaryDist, ~] = bwdist(lobeBoundaries);
fissureCandidates = fissureCandidates & (boundaryDist < params.fissures.maxDistanceToBoundary);

%% Morphological refinement
fprintf('    Refining fissure segmentation...\n');

% Keep structures near lobe boundaries
fissures = imreconstruct(lobeBoundaries, fissureCandidates, 26);

% Remove small isolated structures
fissures = bwareaopen(fissures, params.fissures.minSize);

% Thin to sheet-like structures
fissures = bwskel(fissures);

% Keep only within lung ROI
fissures(masks.binary == 0) = 0;

fprintf('    Fissure segmentation complete\n');

end
