function [airwaysSeg, lowAttenuation] = segmentTrachea(volumes, lungsProcessed, airwaysSeg, lowAttenuation, params)
%SEGMENTTRACHEA Segment trachea and finalize bronchi
%
%   [airwaysSeg, lowAttenuation] = segmentTrachea(volumes, lungsProcessed,
%       airwaysSeg, lowAttenuation, params)
%
%   Performs three-phase trachea processing:
%     Phase 1: Body mask extraction and trachea seeding (2D slice-by-slice)
%     Phase 2: Trachea reconstruction and branchpoint-based separation
%     Phase 3: Bronchi finalization using vesselness and FA filtering
%
%   Inputs:
%       volumes        - Structure with .ct field
%       lungsProcessed - Structure with .binary, .convexhull, .lungsBinFill
%       airwaysSeg     - Structure with .airways, .spurious
%       lowAttenuation - Low-attenuation segmentation mask (air trapping)
%       params         - Processing parameters (requires .trachea, .lungs, .voxel)
%
%   Outputs:
%       airwaysSeg     - Updated: .airways = full tree, .trachea = above carina
%       lowAttenuation - Updated with small airway fragments
%
%   See also: refineAirways, computeFractionalAnisotropy, filterByAnisotropy

fprintf('  Segmenting trachea and finalizing bronchi...\n');

ct = volumes.ct;
lungsbin = lungsProcessed.binary;
lungsBinFill = lungsProcessed.lungsBinFill;
convexhull = lungsProcessed.convexhull;
volSize = size(ct);
tp = params.trachea;

%% Phase 1: Body mask and trachea seeding
fprintf('    Phase 1: Body mask extraction and trachea seeding...\n');
[tracheaSeeds, cavitiesStack] = seedTrachea(ct, lungsbin, lungsBinFill, ...
    volSize, params.lungs.erosionSize, tp);

%% Phase 2: Trachea/bronchi separation at branchpoint
fprintf('    Phase 2: Trachea-bronchi separation at carina...\n');

% Reconstruct trachea from extra-pulmonary seeds
seeds = tracheaSeeds .* 255;
seedsInLungs = imreconstruct(uint8(lungsbin), seeds, 8);
seeds = uint8(seeds & ~seedsInLungs);
tracheaReco = imreconstruct(seeds, cavitiesStack);

% Combine with existing airways, fill holes, mask to ROI
combined = tracheaReco | airwaysSeg.airways;
combined = uint8(imfill(combined, 8, 'holes'));
combined = uint8(imfill(combined, 26, 'holes'));
combined(convexhull == 0) = 0;

% Reconstruct trachea within combined tree (only seed-connected parts)
trachea = combined;
trachea = imreconstruct(seeds, trachea, 26);
trachea(imdilate(lungsBinFill, strel('disk', 1)) ~= 0) = 0;

% Find carina: flip z so search goes from apex downward
tracheaFlip = flip(trachea, 3);
tracheaReduced = bwareaopen(tracheaFlip, tp.minVolumeForSkel, 8);
tracheaReduced = imclose(tracheaReduced, strel('sphere', 2));
skel = bwskel(tracheaReduced, 'MinBranchLength', tp.branchMinLength);
branchPts = bwmorph3(skel, 'branchpoints');

splitSlice = volSize(3);
for k = 1:volSize(3)
    if length(unique(branchPts(:,:,k))) > 1
        splitSlice = k;
        break
    end
end
tracheaFlip(:,:,splitSlice:volSize(3)) = 0;

% Apply split: trachea = above carina, bronchi = below
trachea = combined;
trachea(flip(tracheaFlip, 3) == 0) = 0;
trachea = imreconstruct(seeds, trachea, 26);

bronchi = combined & ~trachea;

% Morphological close-open on bronchi
se = strel('sphere', 2);
bronchiClosed = imdilate(bronchi, se);
bronchiClosed = imerode(bronchiClosed, se);
bronchi = bronchiClosed | bronchi;
bronchi(ct > tp.airwaysMaxHU) = 0;

% Re-combine and mask to ROI
fullTree = trachea | bronchi;
fullTree(convexhull == 0) = 0;

% Remove small fragments from tree → lowAttenuation
smallFragments = fullTree & ~bwareaopen(fullTree, tp.minFragmentVolume);
fullTree = fullTree & ~smallFragments;
lowAttenuation = lowAttenuation | smallFragments;

% HU filter on full tree
fullTree(ct > tp.airwaysMaxHU) = 0;

% Ensure lowAttenuation does not overlap airways
lowAttenuation(fullTree ~= 0) = 0;

%% Phase 3: Bronchi finalization with FA
fprintf('    Phase 3: Bronchi finalization (vesselness + FA)...\n');
bp = tp.bronchi;
voxDim = [params.voxel.px; params.voxel.px; params.voxel.vz];

% --- Pass A: Recover tubular AirSeg fragments near airway endpoints ---
% Compute endpoint proximity map
skelAirways = bwskel(fullTree, 'MinBranchLength', 2);
endpointsImg = bwmorph3(skelAirways, 'endpoints');
endpointsDist = bwdist(endpointsImg, 'quasi-euclidean');
nearEndpoints = endpointsDist > 0 & endpointsDist < bp.endpointDistance;

% Vesselness on lowAttenuation to find tubular candidates
airVesselness = vesselness3D(lowAttenuation, bp.vesselness.scale, ...
    voxDim, bp.vesselness.tau, true);
airCandidates = airVesselness > bp.vesselness.threshold;
airCandidates = bwareaopen(airCandidates, bp.minVesselnessVol, 26);

% Close small gaps
airCandidatesClosed = imdilate(airCandidates, strel('sphere', 1));
airCandidatesClosed = imerode(airCandidatesClosed, strel('sphere', 1));
airCandidates = airCandidatesClosed | airCandidates;

% Keep only candidates near airway endpoints
airCandidates = imreconstruct(nearEndpoints, airCandidates, 26);

% FA filter: keep anisotropic structures (FA > threshold)
[add2Airways, ~] = filterByAnisotropy(airCandidates, bp.FA);
add2Airways = bwareaopen(add2Airways, bp.minRecoverVolume);

% Update segmentations
lowAttenuation = lowAttenuation & ~add2Airways;
fullTree = fullTree | add2Airways;

% Morphological cleanup
fullTreeClosed = imdilate(fullTree, strel('sphere', 1));
fullTreeClosed = imerode(fullTreeClosed, strel('sphere', 1));
fullTree = fullTreeClosed | fullTree;
fullTree = bwareaopen(fullTree, bp.minVolume18, 18);
fullTree = bwareaopen(fullTree, bp.minVolume6, 6);

% --- Pass B: Remove non-anisotropic small airways ---
% Check small fragments (< maxCheckVolume) that have short skeletons
airwayCheck = fullTree & ~bwareaopen(fullTree, bp.maxCheckVolume);
lengthImg = bwskel(airwayCheck);
lengthImg = bwareaopen(lengthImg, bp.skelMinLength);
checkOK = imreconstruct(lengthImg, airwayCheck);
airwayCheck = airwayCheck & ~checkOK;   % Only short-skeleton fragments

% FA filter: keep anisotropic, discard the rest → lowAttenuation
[~, add2AirSeg] = filterByAnisotropy(airwayCheck, bp.FA);
fullTree = fullTree & ~add2AirSeg;
lowAttenuation = lowAttenuation | add2AirSeg;

%% Output
airwaysSeg.airways = fullTree;
airwaysSeg.trachea = trachea;

fprintf('  Trachea segmentation complete\n\n');

end

%% ======================== Helper Functions ========================

function [tracheaSeeds, cavitiesStack] = seedTrachea(ct, lungsbin, ...
    lungsBinFill, volSize, erosionSize, tp)
%SEEDTRACHEA Extract body mask and identify trachea seeds slice-by-slice
%
%   For each axial slice:
%     1. Body mask via HU thresholding + morphological cleanup
%     2. Trachea seeding via air threshold + roundness metric
%     3. Cavity extraction outside lung parenchyma

cavitiesStack = uint8(zeros(volSize));
tracheaSeeds = uint8(zeros(volSize));

for k = 1:volSize(3)
    slice = ct(:,:,k);
    lungsFilled_k = lungsBinFill(:,:,k);
    lungsbin_k = lungsbin(:,:,k);

    %% Body mask extraction
    bodyBin = double(slice >= tp.bodyMaskHU);
    bodyBin = imbinarize(bodyBin);
    se = strel('line', 10, 110);
    bodyBin = imdilate(bodyBin, se);
    se = strel('line', 10, 70);
    bodyBin = imerode(bodyBin, se);
    bodyBin = imdilate(bodyBin, se);

    % Find largest white object (body)
    [labelMap, nObj] = bwlabel(bodyBin);

    if nObj ~= 0
        stats = regionprops(labelMap, 'Area');
        [~, maxIdx] = max([stats.Area]);
        bodyObj = (labelMap == maxIdx);
        bodyObj(:,1) = 0;
        bodyObj(:,end) = 0;

        % Find largest black region (background); complement = filled body
        [bgLabel, nBg] = bwlabel(~bodyObj);
        statsBg = regionprops(bgLabel, 'Area');
        [~, bgIdx] = max([statsBg.Area]);
        mask = ~(bgLabel == bgIdx);
    else
        mask = bodyBin;
    end

    % Smooth body mask with horizontal line
    se = strel('line', 40, 0);
    mask = imerode(mask, se);
    mask = imdilate(mask, se);

    %% Trachea seeding: intensity + sphericity
    cavities = uint8(mask);
    cavities(slice < tp.intensityThreshold & cavities ~= 0) = 255;
    cavities = imbinarize(cavities);
    cavities = imfill(cavities, 'holes');
    [~, L] = bwboundaries(cavities);
    statsL = regionprops(L, 'Area', 'Perimeter');
    trc = uint8(zeros(size(slice)));

    % Only seed trachea in upper half of volume (trachea is superior)
    if k > ceil(volSize(3) / 2)
        for i = 1:length(statsL)
            metric = 4 * pi * statsL(i).Area / statsL(i).Perimeter^2;
            if metric < tp.minRoundness
                L(L == i) = 0;
            end
            if statsL(i).Area < tp.minSeedArea
                L(L == i) = 0;
            end
            trc = uint8(imbinarize(L, 0));
        end
    end

    %% Separate extra-pulmonary cavities from intra-pulmonary
    cavities0 = cavities;
    cavities(lungsbin_k ~= 0) = 0;
    lungsFilled_dilated = imdilate(lungsFilled_k, ...
        strel('disk', erosionSize + 1));
    cavities0(lungsFilled_dilated ~= 0) = 0;
    cavities = imreconstruct(cavities0, cavities, 8);

    tracheaSeeds(:,:,k) = trc;
    cavitiesStack(:,:,k) = cavities;
end

end
