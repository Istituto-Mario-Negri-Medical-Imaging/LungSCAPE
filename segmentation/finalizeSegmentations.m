function segs = finalizeSegmentations(vessels, largeVessels, airways, trachea, ...
    dense, ggo, airSeg, lungsbin, distanceMasks, fissures, volumes, params)
%FINALIZESEGMENTATIONS Cross-segment wall interactions and final cleanup
%
%   segs = finalizeSegmentations(vessels, largeVessels, airways, trachea,
%       dense, ggo, airSeg, lungsbin, distanceMasks, fissures, volumes, params)
%
%   Performs post-segmentation finalization:
%   1. Airway walls computation with LargeAirways merge and validation
%   2. Trachea walls computation with reconstruction
%   3. First masking: remove airway walls from consolidation/GGO
%   4. Large vessel lumen cleanup and masking of consolidation/GGO
%   5. Fibrotic bands recovery via adaptive thresholding
%
%   Inputs:
%       vessels       - Refined vessel segmentation (small/medium)
%       largeVessels  - Large vessel segmentation (from addLargeVessels)
%       airways       - Full airway tree segmentation
%       trachea       - Trachea segmentation (above carina)
%       dense         - Consolidation/dense opacities
%       ggo           - Ground glass opacity
%       airSeg        - Low-attenuation (air trapping) segmentation
%       lungsbin      - Binary lung mask (eroded)
%       distanceMasks - Structure with .close, .far fields
%       fissures      - Binary fissure segmentation
%       volumes       - Structure with .filtered field
%       params        - Processing parameters
%
%   Output:
%       segs - Structure with finalized segmentations:
%           .vessels, .airways, .trachea, .dense, .ggo, .airSeg,
%           .airwayWalls, .tracheaWalls
%
%   See also: refineVessels, createLabelMap

fprintf('  Finalizing segmentations (wall interactions)...\n');

fp = params.finalize;

%% 1. Airway walls
fprintf('    Computing airway walls...\n');

% Extract LargeAirways (diameter > 4) for wider wall computation
[diamAW, ~] = bwdist(single(~airways));
diamAW = 2 * diamAW;
skelAW = bwskel(airways);
diamImgAW = diamAW .* double(skelAW);
largeAirwaysSeed = diamImgAW > fp.airwayWalls.largeAirwayDiam;
largeAirwaysSeed = imdilate(largeAirwaysSeed, strel('sphere', 1));
largeAirwaysReco = imreconstruct(largeAirwaysSeed, airways, 8);
largeAirwaysReco = imreconstruct(distanceMasks.close, largeAirwaysReco, 26);

% Compute distance maps
[airwayDist, ~] = bwdist(single(airways));
[tracheaDist, ~] = bwdist(single(trachea));
[largeAirwayDist, ~] = bwdist(single(largeAirwaysReco));

% Airway walls: normal airways + large airways (wider walls)
airwayWalls = airwayDist < fp.airwayWalls.normalWallDist & airwayDist > 0;
largeAirwayWalls = largeAirwayDist < fp.airwayWalls.largeWallDist & largeAirwayDist > 0;
airwayWalls = airwayWalls | largeAirwayWalls;

% Trachea walls: merge with connected airway walls, fill, reconstruct from trachea seed
tracheaWalls = tracheaDist < fp.airwayWalls.tracheaWallDist & tracheaDist > 0;
tracheaWalls = tracheaWalls | imreconstruct(tracheaWalls, airwayWalls, 8);
tracheaWalls = imfill(tracheaWalls, 8, 'holes');
tracheaWalls = imreconstruct(uint8(trachea), uint8(tracheaWalls), 8);
tracheaWalls = logical(tracheaWalls) & ~trachea;

% Remove trachea walls from airway walls (no overlap)
airwayWalls = airwayWalls & ~tracheaWalls;

% Validate airway walls: must be connected to airway lumen
airwayWallsFilled = imfill(airwayWalls, 8, 'holes');
airwayWallsFilled = imreconstruct(uint8(airways), uint8(airwayWallsFilled), 8);
airwayWalls = airwayWalls & logical(airwayWallsFilled);
airwayWalls = airwayWalls & ~airways;

%% 2. First masking: consolidation/GGO -= airway walls
fprintf('    First masking: removing airway walls from pathology...\n');
dense = dense & ~airwayWalls;
ggo = ggo | dense;
ggo = ggo & ~airwayWalls;

%% 3. Large vessel lumen cleanup and masking
fprintf('    Masking pathology with large vessel lumen...\n');

% Remove large vessels from airway walls, airways, far region
largeVessels = largeVessels & ~airwayWalls;
largeVessels = largeVessels & ~airways;
largeVessels = largeVessels & ~distanceMasks.far;
largeVessels = bwareaopen(largeVessels, 30, 6);

dense = dense & ~largeVessels;
ggo   = ggo   & ~largeVessels;

%% 5. Reclassify healthy tissue with HU above threshold
fprintf('    Reclassifying high-HU healthy tissue as pathological...\n');

% Compute healthy tissue (everything not otherwise classified)
healthy = zeros(size(lungsbin), 'uint8');
healthy(lungsbin ~= 0) = 1;
healthy(ggo ~= 0) = 0;
healthy(dense ~= 0) = 0;
healthy(airSeg ~= 0) = 0;
healthy(vessels ~= 0) = 0;
healthy(largeVessels ~= 0) = 0;
healthy(airways ~= 0) = 0;
healthy(airwayWalls ~= 0) = 0;
if any(trachea(:))
    healthy(trachea ~= 0) = 0;
end
healthy(tracheaWalls ~= 0) = 0;

% Healthy voxels above healthyMaxHU → pathological
% Exclude voxels within 1 voxel of vessel borders (partial volume guard)
ct = volumes.ct;
[vesselBorderDist, ~] = bwdist(vessels | largeVessels);
healthyToInjury = logical(healthy) & (ct > params.thresholds.healthyMaxHU) ...
    & (vesselBorderDist > fp.vesselGuardDist);
clear vesselBorderDist;

% Split reclassified voxels: GGO (≤ ggoMaxHU) vs Dense (> ggoMaxHU)
ggo = ggo | (healthyToInjury & (ct <= params.thresholds.ggoMaxHU));
dense = dense | (healthyToInjury & (ct > params.thresholds.ggoMaxHU));
healthy(healthyToInjury) = 0;

% Enforce GGO upper HU limit on ALL GGO voxels (including pre-existing)
ggoToDense = (ggo ~= 0) & (ct > params.thresholds.ggoMaxHU);
dense = dense | ggoToDense;
ggo = ggo & ~ggoToDense;

%% 6. Fibrotic bands recovery
fprintf('    Recovering fibrotic bands in healthy tissue...\n');
fb = fp.fibrotic;
sigma = 0.3 * ((fb.kernelSize - 1) * 0.5 - 1) + 0.8;

% Adaptive thresholding on healthy tissue
roi = double(lungsbin > 0);
imSrc = double(imcomplement(volumes.filtered)) .* roi;

lungblur = imgaussfilt(imSrc, sigma, 'FilterSize', fb.kernelSize);
roilungs = imgaussfilt(roi, sigma, 'FilterSize', fb.kernelSize);
smooth = lungblur ./ roilungs;
fibroticBands = 255 .* (imSrc <= (smooth - fb.C));
fibroticBands = uint8(fibroticBands .* roi);
fibroticBands = logical(fibroticBands) & logical(healthy);
fibroticBands = bwareaopen(fibroticBands, fb.minVol, 6);

% Exclude vessel regions (small/medium + large)
fibroticBands = fibroticBands & ~vessels;
fibroticBands = fibroticBands & ~largeVessels;

% Exclude fissure regions
if any(fissures(:))
    fibroticBands = fibroticBands & ~fissures;
end

% Add fibrotic bands to consolidation
dense = dense | fibroticBands;

%% 7. Final cleanup: size filter + morphological opening (GGO only)
fprintf('    Final cleanup: size filter and GGO opening...\n');
dense = bwareaopen(dense, fp.finalMinVol.dense, 6);
ggo   = bwareaopen(ggo,   fp.finalMinVol.ggo,   6);

% Morphological opening on GGO to remove thin shells and isolated voxels.
% Not applied to dense: reticulations are thin structures that opening would damage.
ggo = imopen(ggo, strel('sphere', 1));

%% Output
segs.vessels = vessels | largeVessels;
segs.airways = airways;
segs.trachea = trachea;
segs.dense = dense;
segs.ggo = ggo;
segs.airSeg = airSeg;
segs.airwayWalls = airwayWalls;
segs.tracheaWalls = tracheaWalls;

fprintf('    Finalization complete\n');

end
