function vesselsRefined = postTracheaVesselRefinement(vesselsRefined, vesselsSeg, ...
    volumes, lungsProcessed, airwaysSeg, injurySeg, params)
%POSTTRACHEAVESSELREFINEMENT Second pass vessel recovery near airways (post-trachea)
%
%   vesselsRefined = postTracheaVesselRefinement(vesselsRefined, vesselsSeg,
%       volumes, lungsProcessed, airwaysSeg, injurySeg, params)
%
%   After trachea finalization, performs a second pass of vessel recovery
%   near the airways and cyst removal. This step uses the initial Jerman
%   output (rawThreshold) which includes vessels that were discarded during
%   the graph-based loop/end classification, but may still be valid
%   bronchovascular structures.
%
%   Steps:
%       1. Recover vessels along airways using rawThreshold segmentation
%          (constrained to arterial territory — bronchovascular bundle anatomy)
%       2. Remove vessels near pathological airways (honeycombing cysts)
%       3. Remove small vessels near AirSeg-derived cysts in compromised tissue
%
%   Inputs:
%       vesselsRefined  - Current refined vessel segmentation (logical)
%       vesselsSeg      - Structure with .rawThreshold field from segmentVessels
%       volumes         - Structure with .ct, .arteries fields
%       lungsProcessed  - Structure with .binary field
%       airwaysSeg      - Structure with .airways, .spurious fields (post-trachea)
%       injurySeg       - Structure with .dense, .lowAttenuation fields
%       params          - Processing parameters
%
%   Output:
%       vesselsRefined  - Updated vessel segmentation (logical)
%
%   See also: segmentVessels, refineVessels

fprintf('  Post-trachea vessel refinement...\n');

pt = params.vessels.postTrachea;

%% Get the raw Jerman threshold (before loop/end reconstruction)
rawThreshold = vesselsSeg.rawThreshold;

%% Compute diameter and skeleton on raw threshold
diameter = 2 * bwdist(~rawThreshold);
skeletonImage = bwskel(rawThreshold);
diameterImage = diameter .* double(skeletonImage);
diameterImageMain = diameterImage >= pt.mainDiamThr;

%% Recover vessels along airways (bronchovascular bundles)
fprintf('    Recovering vessels along airways...\n');

[airwaysDist, ~] = bwdist(single(airwaysSeg.airways));
airwaysDist(diameterImageMain == 0) = 0;

% Vessels within max distance but beyond min distance from airways
vesselsCloseToAirways = airwaysDist < pt.airwaysMaxDist & ...
    airwaysDist > pt.airwaysMinDist;

% Reconstruct full vessels from the rawThreshold segmentation.
% Constrained to arterial territory: same rationale as in segmentVessels —
% the bronchovascular bundle relationship holds only for arteries.
arteriesMask = imdilate(volumes.arteries, strel('sphere', 1));
vesselsCloseToAirways = imreconstruct(vesselsCloseToAirways, rawThreshold, 8);
vesselsCloseToAirways = vesselsCloseToAirways & arteriesMask;

%% Remove vessels near pathological airways (honeycombing cysts)
fprintf('    Removing vessels near honeycombing cysts...\n');

vesselsNearCysts = false(size(rawThreshold));
if isfield(airwaysSeg, 'spurious') && any(airwaysSeg.spurious(:))
    % Remove trachea airways from spurious before distance computation
    pathologicalAirways = airwaysSeg.spurious & ~airwaysSeg.airways;
    pathologicalAirwaysMain = bwareaopen(pathologicalAirways, pt.cystsMinVol);
    [pathAirwaysDist, ~] = bwdist(single(pathologicalAirwaysMain));

    vesselsNearCysts = pathAirwaysDist < pt.cystsMaxDist & pathAirwaysDist > 0;
end

%% Remove small vessels near AirSeg-derived cysts in compromised tissue
fprintf('    Removing vessels near AirSeg cysts...\n');

vesselsNearCysts2 = false(size(rawThreshold));
if isfield(injurySeg, 'lowAttenuation') && any(injurySeg.lowAttenuation(:))
    % Use AirSeg (low-attenuation) to find additional cysts
    airSegMain = bwareaopen(injurySeg.lowAttenuation, pt.cysts2MinAirSegVol, 4);
    [airSegDist, ~] = bwdist(single(airSegMain));
    cysts2proximity = airSegDist < pt.cysts2MaxDist & airSegDist > 0;

    % Seed from consolidation
    if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
        cysts2proximity = imreconstruct(logical(injurySeg.dense), cysts2proximity);
    end

    % Only small vessel fragments in compromised tissue
    checkSmall = vesselsRefined & ~bwareaopen(vesselsRefined, pt.cysts2SmallVesselMaxVol, 6);
    checkSmall2 = bwareaopen(checkSmall & cysts2proximity, pt.cysts2MinOverlapVol);
    vesselsNearCysts2 = imreconstruct(checkSmall2, checkSmall);

    % Restrict to compromised tissue
    if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
        vesselsNearCysts2 = vesselsNearCysts2 & logical(injurySeg.dense);
    end
end

%% Apply corrections
% Remove cyst-adjacent vessels first, so airways recovery doesn't re-add them
vesselsCloseToAirways = vesselsCloseToAirways & ~vesselsNearCysts;

% Update vessels
vesselsRefined = vesselsRefined | vesselsCloseToAirways;
vesselsRefined = vesselsRefined & ~vesselsNearCysts;
vesselsRefined = vesselsRefined & ~vesselsNearCysts2;

fprintf('    Post-trachea vessel refinement complete\n');

end
