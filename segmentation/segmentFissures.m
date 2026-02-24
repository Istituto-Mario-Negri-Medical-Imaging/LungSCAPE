function [fissures, fissureInfo] = segmentFissures(volumes, masks, lobes, params)
%SEGMENTFISSURES Segment lung fissures as morphological lobe boundaries
%
%   [fissures, fissureInfo] = segmentFissures(volumes, masks, lobes, params)
%
%   Extracts interlobar fissures by:
%   1. Computing morphological boundary of each lobe (bwmorph3 'remove')
%   2. Removing external lung boundaries
%   3. Dilating and cleaning small fragments
%
%   Inputs:
%       volumes - Structure with CT volumes (unused, kept for interface)
%       masks   - Structure with .binary field (eroded lung mask)
%       lobes   - Lobe segmentation (uint8: 1-5)
%       params  - Processing parameters (.fissures.dilationSize, .minSize)
%
%   Outputs:
%       fissures    - Binary fissure segmentation (logical)
%       fissureInfo - Structure with:
%           .externalBoundary - External lung boundary mask
%
%   See also: refineLobeSegmentation

fprintf('  Segmenting lung fissures...\n');

% Check if lobes are available
if isempty(lobes) || all(lobes(:) == 0)
    warning('No lobe segmentation available. Skipping fissure segmentation.');
    fissures = false(size(masks.binary));
    fissureInfo.externalBoundary = false(size(masks.binary));
    return;
end

%% Extract fissures as morphological boundaries of each lobe
fprintf('    Extracting lobe boundaries via bwmorph3...\n');
fissures = false(size(masks.binary));
lobeLabels = unique(lobes(:));
lobeLabels(lobeLabels == 0) = [];

for i = 1:length(lobeLabels)
    lobeMask = lobes == lobeLabels(i);
    lobeBorder = bwmorph3(lobeMask, 'remove');
    fissures = fissures | lobeBorder;
end

%% Remove external lung boundaries
fprintf('    Removing external lung boundaries...\n');
externalBoundary = imfill(masks.binary, 8, 'holes');
externalBoundary = bwmorph3(externalBoundary, 'remove');
fissures = fissures & ~externalBoundary;

fissureInfo.externalBoundary = externalBoundary;

%% Morphological refinement
fprintf('    Refining fissure segmentation...\n');
fissures = imdilate(fissures, strel('disk', params.fissures.dilationSize));
fissures = bwareaopen(fissures, params.fissures.minSize, 26);

fprintf('    Fissure segmentation complete\n');

end
