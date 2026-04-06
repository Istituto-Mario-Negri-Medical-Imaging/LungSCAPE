function vesselTypeMap = classifyVesselsByType(vessels, arteries, veins, params, typeMapJerman)
%CLASSIFYVESSELSBYTYPE Classify vessel voxels as arteries, veins, or undetermined
%
%   vesselTypeMap = classifyVesselsByType(vessels, arteries, veins, params)
%
%   Classifies each voxel of the refined vessel segmentation as artery (1),
%   vein (2), or undetermined (3) using a skeleton-based two-phase approach:
%
%   Phase 1 — Skeleton-based propagation (primary criterion)
%       The vessel binary is reduced to its medial-axis skeleton (bwskel).
%       Labels are propagated from artery/vein anchor nodes (skeleton voxels
%       that directly overlap the TotalSegmentator masks) through the skeleton
%       graph via morphological reconstruction. Because the skeleton is
%       1-voxel wide, two vessel branches that merely touch in the full binary
%       volume — but are anatomically separate — are correctly split at the
%       skeleton level, whereas the previous volumetric imreconstruct would
%       propagate both labels through the entire connected region.
%
%       Skeleton nodes reachable from both anchor types (residual conflicts at
%       true bifurcation junctions) are resolved by shortest-path distance on
%       the skeleton graph using binaryImageGraph3Weighted with two virtual
%       super-source nodes. This respects vascular tree topology rather than
%       Euclidean geometry.
%
%       Skeleton labels are projected to the full vessel volume by assigning
%       each vessel voxel the label of its nearest skeleton voxel (bwdist).
%
%   Phase 2 — Euclidean distance fallback (isolated and unreached voxels)
%       Voxels not resolved by Phase 1 (disconnected micro-components, blobs
%       too small for skeletonisation, skeleton nodes with no reachable anchor)
%       are classified by nearest-neighbour Euclidean distance to the TS masks,
%       within maxDistVox. Voxels beyond maxDistVox from both masks are marked
%       undetermined (3). This threshold is kept conservative intentionally:
%       peripheral unreached voxels are the most topologically ambiguous and
%       expanding their coverage risks contaminating the artery/vein classes.
%
%   The TotalLabelMap (label 5 = all vessels) is left unchanged for
%   backward compatibility. This function produces a separate VesselsTypeMap.
%
%   Inputs:
%       vessels  - Binary refined vessel segmentation (logical)
%       arteries - TS artery mask (uint8/logical), volumes.arteries
%       veins    - TS vein mask   (uint8/logical), volumes.veins
%       params   - Processing parameters
%
%   Output:
%       vesselTypeMap - uint8 volume:
%           0 - not a vessel
%           1 - artery
%           2 - vein
%           3 - undetermined vessel
%
%   See also: classifyVesselsByDiameter, binaryImageGraph3Weighted, saveResults

fprintf('  Classifying vessels by type (artery / vein)...\n');

% typeMapJerman is optional (Strategy A prior from segmentVessels).
% When absent the function behaves identically to the pre-Strategy-A version.
if nargin < 5
    typeMapJerman = [];
end

maxDistVox = params.vessels.typeClassification.maxDistVox;
vessMask   = logical(vessels);

%% Precompute Euclidean distances from TS masks
% Computed early: needed for both Phase 2 fallback and anchor snap in Phase 1.
[distArt,  ~] = bwdist(logical(arteries));
[distVein, ~] = bwdist(logical(veins));

%% Phase 1: Skeleton-based label propagation

% --- Skeleton ---
fprintf('    Computing vessel skeleton...\n');
skel = bwskel(vessMask);

% Degenerate case: skeleton empty (vessMask too thin/small for bwskel).
% Degrade gracefully to the full-volume Euclidean fallback.
if ~any(skel(:))
    fprintf('    [Warning] Vessel skeleton is empty — using distance-only classification.\n');
    vesselTypeMap = distanceFallback(vessMask, distArt, distVein, maxDistVox);
    printSummary(vesselTypeMap, 0);
    fprintf('    Vessel type classification complete\n');
    return;
end

% Nearest-skeleton index map.
% idxNearestSkel(i) = linear index of the nearest skeleton voxel to voxel i.
% Computed once; reused for both anchor snap (below) and volume projection.
% Note: bwdist(BW) gives distance to nearest nonzero voxel in BW, and its
% second output gives the linear index of that nearest voxel.
[~, idxNearestSkel] = bwdist(skel);

% --- Anchor nodes on skeleton ---
% Build anchor source masks in volume space, then project onto skeleton.
%
% When typeMapJerman is available (Strategy A), voxels labeled in the
% preliminary binary take precedence over TS masks: they carry type
% information derived from the Jerman filter response on each individual
% vascular system rather than from TS mask spatial proximity alone.
% TS masks supplement the Jerman prior for voxels added to the final binary
% by downstream steps (adaptive thresholding, large vessel recovery, etc.)
% that were not present in the preliminary binary.
%
% Without typeMapJerman the anchor source falls back to TS masks only,
% preserving the pre-Strategy-A behaviour.
if ~isempty(typeMapJerman)
    % Primary: Jerman-labeled voxels present in the final binary
    artPrimary  = vessMask & (typeMapJerman == 1);
    veinPrimary = vessMask & (typeMapJerman == 2);

    % Supplementary: TS masks restricted to voxels NOT covered by the prior
    notInJerman = vessMask & (typeMapJerman == 0);
    artSuppl    = logical(arteries) & notInJerman;
    veinSuppl   = logical(veins)    & notInJerman;

    artSource  = artPrimary  | artSuppl;
    veinSource = veinPrimary | veinSuppl;
else
    artSource  = vessMask & logical(arteries);
    veinSource = vessMask & logical(veins);
end

artSkelAnchors  = skel & artSource;
veinSkelAnchors = skel & veinSource;

% Snap fallback: if no source voxel lies on the skeleton, map each source
% voxel to its nearest skeleton voxel.
artSkelAnchors  = snapAnchorsToSkeleton(artSkelAnchors,  artSource,  skel, idxNearestSkel);
veinSkelAnchors = snapAnchorsToSkeleton(veinSkelAnchors, veinSource, skel, idxNearestSkel);

% Resolve overlap: if a skeleton node is claimed by both types, assign it
% to the Euclidean-closer TS mask.
both = artSkelAnchors & veinSkelAnchors;
if any(both(:))
    artSkelAnchors( both & (distVein <= distArt)) = false;
    veinSkelAnchors(both & (distArt  <  distVein)) = false;
end

% --- Label propagation through the 1-voxel-wide skeleton ---
% imreconstruct on the sparse skeleton binary is far less susceptible to
% label leakage than on the full volume: branches that are anatomically
% separate but touch in the vessel binary typically have disjoint skeletons.
artSkelReach  = false(size(skel));
veinSkelReach = false(size(skel));

if any(artSkelAnchors(:))
    artSkelReach  = logical(imreconstruct(artSkelAnchors,  skel, 26));
end
if any(veinSkelAnchors(:))
    veinSkelReach = logical(imreconstruct(veinSkelAnchors, skel, 26));
end

% --- Classify skeleton nodes ---
onlyArtSkel    =  artSkelReach & ~veinSkelReach;   % unambiguous artery
onlyVeinSkel   = veinSkelReach & ~artSkelReach;    % unambiguous vein
conflictedSkel =  artSkelReach &  veinSkelReach;   % true junction / merge artefact
% unreached skeleton nodes: skelLabel stays 0 → handled by Phase 2

nConflSkel = sum(conflictedSkel(:));

% --- Path-distance tiebreaker for conflicted skeleton nodes ---
% Build the skeleton graph once and resolve all conflicts vectorially.
if nConflSkel > 0
    fprintf('    Resolving %d conflicted skeleton nodes via shortest-path distance...\n', ...
        nConflSkel);
    [onlyArtSkel, onlyVeinSkel] = resolveConflicts( ...
        onlyArtSkel, onlyVeinSkel, conflictedSkel, ...
        artSkelAnchors, veinSkelAnchors, skel);
end

% --- Skeleton label image: 1=artery, 2=vein, 0=unreached ---
skelLabel = zeros(size(skel), 'uint8');
skelLabel(onlyArtSkel)  = 1;
skelLabel(onlyVeinSkel) = 2;

% --- Project skeleton labels to full vessel volume ---
% Each voxel in vessMask inherits the label of its nearest skeleton voxel.
vesselTypeMap  = zeros(size(vessels), 'uint8');
vessMaskIdx    = find(vessMask);
nearestSkelIdx = idxNearestSkel(vessMaskIdx);       % linear indices into volume
vesselTypeMap(vessMaskIdx) = skelLabel(nearestSkelIdx);

% Voxels projected onto an unreached skeleton node (skelLabel = 0) remain
% zero here and are resolved by the Phase 2 Euclidean fallback below.
unreached = vessMask & (vesselTypeMap == 0);

%% Phase 2: Euclidean distance fallback for unreached voxels

unreachArtery = unreached & (distArt  <= distVein) & (distArt  <= maxDistVox);
unreachVein   = unreached & (distVein <  distArt)  & (distVein <= maxDistVox);
unreachUndet  = unreached & ~unreachArtery & ~unreachVein;

vesselTypeMap(unreachArtery) = 1;
vesselTypeMap(unreachVein)   = 2;
vesselTypeMap(unreachUndet)  = 3;

%% Summary
printSummary(vesselTypeMap, nConflSkel);
fprintf('    Vessel type classification complete\n');

end


%% -----------------------------------------------------------------------
%  Helper: path-distance tiebreaker for conflicted skeleton nodes
%  -----------------------------------------------------------------------
function [onlyArtSkel, onlyVeinSkel] = resolveConflicts( ...
    onlyArtSkel, onlyVeinSkel, conflictedSkel, ...
    artSkelAnchors, veinSkelAnchors, skel)
%RESOLVECONFLICTS Resolve skeleton nodes reachable from both artery and vein.
%
%   Constructs the weighted skeleton graph via binaryImageGraph3Weighted and
%   adds two virtual super-source nodes (one per type) with zero-cost edges
%   to the respective anchor nodes. A single Dijkstra call per super-source
%   yields the shortest-path distance from each anchor set to every skeleton
%   node. Each conflicted node is assigned to whichever type is closer along
%   the skeleton tree.
%
%   Edge weights in binaryImageGraph3Weighted equal sqrt(dx²+dy²+dz²) + 1
%   (the +1 comes from the binary map value at the neighbour). This constant
%   per-hop offset does not affect the relative ranking between artery and
%   vein path distances.

    [gSkel, skelNodeNums] = binaryImageGraph3Weighted(skel, 26);
    nSkelNodes = numnodes(gSkel);

    artAnchorIDs  = unique(nonzeros(skelNodeNums( artSkelAnchors)));
    veinAnchorIDs = unique(nonzeros(skelNodeNums(veinSkelAnchors)));

    % Safety: if both anchor sets are empty the graph has no labelled seeds;
    % leave all conflicted nodes as unreached (skelLabel = 0, Phase 2 fallback).
    if isempty(artAnchorIDs) && isempty(veinAnchorIDs)
        return;
    end

    % Add two virtual super-source nodes (indices nSkelNodes+1, nSkelNodes+2)
    % connected to their respective anchor nodes with zero-cost edges.
    gAug      = addnode(gSkel, 2);
    superArt  = nSkelNodes + 1;
    superVein = nSkelNodes + 2;

    if ~isempty(artAnchorIDs)
        gAug = addedge(gAug, ...
            repmat(superArt,  numel(artAnchorIDs),  1), artAnchorIDs, ...
            zeros(numel(artAnchorIDs),  1));
    end
    if ~isempty(veinAnchorIDs)
        gAug = addedge(gAug, ...
            repmat(superVein, numel(veinAnchorIDs), 1), veinAnchorIDs, ...
            zeros(numel(veinAnchorIDs), 1));
    end

    % Shortest-path distances from each super-source to all skeleton nodes.
    % distances() returns a 1 × numnodes(gAug) row vector.
    dFromArt  = distances(gAug, superArt);
    dFromVein = distances(gAug, superVein);
    dFromArt  = dFromArt( 1:nSkelNodes);   % discard super-source row entries
    dFromVein = dFromVein(1:nSkelNodes);

    % Vectorised resolution: for each conflicted voxel, read its graph node
    % ID and compare the two shortest-path distances.
    conflictedIdx = find(conflictedSkel);
    nodeIDs       = skelNodeNums(conflictedIdx);
    valid         = nodeIDs > 0;              % safety guard (should always be true)

    dA = dFromArt( nodeIDs(valid));
    dV = dFromVein(nodeIDs(valid));
    confValidIdx = conflictedIdx(valid);

    onlyArtSkel( confValidIdx(dA <= dV)) = true;
    onlyVeinSkel(confValidIdx(dA >  dV)) = true;
end


%% -----------------------------------------------------------------------
%  Helper: snap TS anchors to skeleton when no direct overlap exists
%  -----------------------------------------------------------------------
function anchors = snapAnchorsToSkeleton(anchors, tsMask, skel, idxNearestSkel)
%SNAPANCHORSTOSKELETON Map TS anchor voxels to their nearest skeleton voxel.
%
%   Called only when no TS voxel coincides with the skeleton (anchors is
%   all-false). Maps each TS mask voxel to the nearest skeleton voxel via the
%   precomputed idxNearestSkel index map, then restricts the result to actual
%   skeleton voxels (skel & snapMask).

    if any(anchors(:)) || ~any(tsMask(:))
        return;   % already have direct overlap, or no TS mask at all
    end

    tsMaskIdx      = find(tsMask);
    nearestSkelIdx = idxNearestSkel(tsMaskIdx);   % nearest skeleton voxel per TS voxel
    snapMask       = false(size(skel));
    snapMask(nearestSkelIdx) = true;
    anchors = skel & snapMask;
end


%% -----------------------------------------------------------------------
%  Helper: full Euclidean-distance classification (skeleton-empty fallback)
%  -----------------------------------------------------------------------
function typeMap = distanceFallback(vessMask, distArt, distVein, maxDistVox)
%DISTANCEFALLBACK Pure distance-based fallback used when bwskel returns empty.
%
%   All vessel voxels are classified by nearest-neighbour distance to the TS
%   masks, within maxDistVox. Voxels beyond that threshold from both masks are
%   marked undetermined (3). This matches the Phase 2 logic applied to the
%   full vessel mask.

    typeMap       = zeros(size(vessMask), 'uint8');
    unreachArtery = vessMask & (distArt  <= distVein) & (distArt  <= maxDistVox);
    unreachVein   = vessMask & (distVein <  distArt)  & (distVein <= maxDistVox);
    unreachUndet  = vessMask & ~unreachArtery & ~unreachVein;
    typeMap(unreachArtery) = 1;
    typeMap(unreachVein)   = 2;
    typeMap(unreachUndet)  = 3;
end


%% -----------------------------------------------------------------------
%  Helper: print classification statistics
%  -----------------------------------------------------------------------
function printSummary(vesselTypeMap, nConflSkel)
%PRINTSUMMARY Log per-class voxel counts and conflict resolution note.

    nArt   = sum(vesselTypeMap(:) == 1);
    nVein  = sum(vesselTypeMap(:) == 2);
    nUndet = sum(vesselTypeMap(:) == 3);
    nTot   = nArt + nVein + nUndet;

    fprintf('    Arteries:      %d voxels (%.1f%%)\n', nArt,   100*nArt   / max(nTot,1));
    fprintf('    Veins:         %d voxels (%.1f%%)\n', nVein,  100*nVein  / max(nTot,1));
    fprintf('    Undetermined:  %d voxels (%.1f%%)\n', nUndet, 100*nUndet / max(nTot,1));
    if nConflSkel > 0
        fprintf('    [Note] %d conflicted skeleton nodes resolved by path distance\n', nConflSkel);
    end
end
