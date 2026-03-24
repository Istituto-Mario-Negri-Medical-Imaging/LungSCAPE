function vesselTypeMap = classifyVesselsByType(vessels, arteries, veins, params)
%CLASSIFYVESSELSBYTYPE Classify vessel voxels as arteries, veins, or undetermined
%
%   vesselTypeMap = classifyVesselsByType(vessels, arteries, veins, params)
%
%   Classifies each voxel of the refined vessel segmentation as artery (1),
%   vein (2), or undetermined (3) using a two-phase approach:
%
%   Phase 1 — Topological propagation (primary criterion)
%       Vessel voxels that directly overlap the TotalSegmentator artery or
%       vein mask serve as "anchors" of known type. Labels are then
%       propagated from these anchors through the connected vessel lumen
%       using morphological reconstruction (imreconstruct). This correctly
%       classifies peripheral voxels that lie outside the TS mask boundary
%       by following the vessel tree topology rather than geometric distance
%       — a voxel near an artery but connected via vessel lumen to a vein
%       is correctly classified as vein.
%
%   Phase 2 — Distance tiebreaker (conflicts and isolated components)
%       Voxels reachable from both artery and vein anchors (segmentation
%       merge artefact) are resolved by nearest-neighbour distance to the
%       respective TS masks. Voxels not reachable from any anchor and
%       within maxDistVox of a TS mask get a distance-based label. Voxels
%       beyond maxDistVox from both masks are marked undetermined (3).
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
%   See also: classifyVesselsByDiameter, saveResults

fprintf('  Classifying vessels by type (artery / vein)...\n');

maxDistVox = params.vessels.typeClassification.maxDistVox;
vessMask   = logical(vessels);

%% Phase 1: Topological propagation from TS anchor voxels

% Anchors: vessel voxels that directly overlap the TS masks.
% These are the high-confidence seeds from which labels propagate.
arteryAnchors = vessMask & logical(arteries);
veinAnchors   = vessMask & logical(veins);

% Propagate each label through the connected vessel lumen.
% imreconstruct expands the seed into any connected region of the mask,
% following 26-connectivity through the vessel binary volume.
arteryReach = false(size(vessMask));
veinReach   = false(size(vessMask));

if any(arteryAnchors(:))
    arteryReach = logical(imreconstruct(arteryAnchors, vessMask, 26));
end
if any(veinAnchors(:))
    veinReach = logical(imreconstruct(veinAnchors, vessMask, 26));
end

% First-pass assignment from topology alone
onlyArtery = arteryReach & ~veinReach;   % unambiguous artery
onlyVein   = veinReach   & ~arteryReach; % unambiguous vein
conflicted = arteryReach & veinReach;    % reachable from both (merge artefact)
unreached  = vessMask & ~arteryReach & ~veinReach; % disconnected from all anchors

%% Phase 2: Distance tiebreaker for conflicts and unreached voxels

[distArt,  ~] = bwdist(logical(arteries));
[distVein, ~] = bwdist(logical(veins));

% Conflicted voxels: nearest TS mask wins
conflictArtery = conflicted & (distArt <= distVein);
conflictVein   = conflicted & (distVein <  distArt);

% Unreached voxels: distance-based if within maxDistVox, else undetermined
unreachArtery = unreached & (distArt  <= distVein) & (distArt  <= maxDistVox);
unreachVein   = unreached & (distVein <  distArt)  & (distVein <= maxDistVox);
unreachUndet  = unreached & ~unreachArtery & ~unreachVein;

%% Assemble output map
vesselTypeMap = zeros(size(vessels), 'uint8');
vesselTypeMap(onlyArtery)    = 1;
vesselTypeMap(onlyVein)      = 2;
vesselTypeMap(conflictArtery)= 1;
vesselTypeMap(conflictVein)  = 2;
vesselTypeMap(unreachArtery) = 1;
vesselTypeMap(unreachVein)   = 2;
vesselTypeMap(unreachUndet)  = 3;

%% Summary
nArt   = sum(vesselTypeMap(:) == 1);
nVein  = sum(vesselTypeMap(:) == 2);
nUndet = sum(vesselTypeMap(:) == 3);
nConfl = sum(conflicted(:));
nTot   = nArt + nVein + nUndet;

fprintf('    Arteries:      %d voxels (%.1f%%)\n', nArt,   100*nArt/max(nTot,1));
fprintf('    Veins:         %d voxels (%.1f%%)\n', nVein,  100*nVein/max(nTot,1));
fprintf('    Undetermined:  %d voxels (%.1f%%)\n', nUndet, 100*nUndet/max(nTot,1));
if nConfl > 0
    fprintf('    [Note] %d conflict voxels resolved by distance (likely segmentation merges)\n', nConfl);
end
fprintf('    Vessel type classification complete\n');

end
