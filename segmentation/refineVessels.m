function [vesselsRefined, largeVessels, addConsolidation] = refineVessels(vesselsSeg, volumes, ...
    lungsProcessed, airwaysSeg, injurySeg, distanceMasks, fissures, params)
%REFINEVESSELS Advanced vessel post-processing and refinement
%
%   [vesselsRefined, largeVessels, addConsolidation] = refineVessels(vesselsSeg, volumes,
%       lungsProcessed, airwaysSeg, injurySeg, distanceMasks, fissures, params)
%
%   Performs vessel refinement:
%   1. Adaptive thresholding to recover vessels in healthy tissue
%   2. Fissure cleaning (remove small vessels in fissure regions)
%   3. Progressive leakage detection (5-pass distal region filtering)
%   4. Eigenvector orientation filter (horizontal leakage removal)
%   5. Seed-based cleanup of residual leakages
%   6. Large consolidation reclassification (dense tissue misclassified as vessels)
%   7. Large vessel reconstruction (TS artery+vein masks via morphological reconstruction)
%
%   Inputs:
%       vesselsSeg     - Structure with .final, .vesselness fields
%       volumes        - Structure with .filtered, .ct, .arteries, .veins fields
%       lungsProcessed - Structure with .binary, .original fields
%       airwaysSeg     - Structure with .airways field
%       injurySeg      - Structure with .dense, .ggo fields
%       distanceMasks  - Structure with .close, .far fields
%       fissures       - Binary fissure segmentation (logical)
%       params         - Processing parameters
%
%   Outputs:
%       vesselsRefined  - Final refined vessel segmentation (logical)
%       largeVessels    - Large vessels recovered from intensity (logical)
%       addConsolidation - Dense tissue reclassified from vessels (logical)
%
%   See also: segmentVessels, classifyVesselsByDiameter

fprintf('  Post-processing vessel segmentation...\n');

lungsbin = lungsProcessed.binary;
lc = params.vessels.lastCorr;

%% Initialize working copy
vessels = vesselsSeg.final;

%% Honeycombing walls masking (narrow bright structures mask to honeycombing region)
fprintf('    Masking honeycombing walls...\n');
cv = params.vessels.completeVessels;

honeycombingWallsMask = false(size(lungsbin));
if isfield(injurySeg, 'lowAttenuation') && any(injurySeg.lowAttenuation(:))
    honeycombingLargeAirSeg = bwareaopen(injurySeg.lowAttenuation, cv.honeycomb.minAirSegVol);
    honeycombingLargeAirSeg = bwareaopen(honeycombingLargeAirSeg, cv.honeycomb.minAreaVol, 8);
    [airSegDist, ~] = bwdist(single(honeycombingLargeAirSeg));
    honeycombingWallsMask = airSegDist < cv.honeycomb.maxDist & airSegDist > 0;
end

% Narrow bright structures to honeycombing walls region
brightStructuresMasked = false(size(lungsbin));
if isfield(injurySeg, 'brightstructures') && any(injurySeg.brightstructures(:))
    brightStructuresMasked = logical(injurySeg.brightstructures) & honeycombingWallsMask;
    brightStructuresMasked = bwareaopen(brightStructuresMasked, cv.honeycomb.minOutThrVol);
end

%% Adaptive thresholding for complete vessels in healthy tissue
fprintf('    Applying adaptive thresholding in healthy tissue...\n');

masks.lungs = lungsProcessed.original;
masks.lungsBinary = lungsbin;
masks.airways = airwaysSeg.airways;
masks.far = distanceMasks.far;

completeVessels = adaptiveVesselThreshold(volumes.filtered, masks, params);

% Fibrotic check: near-consolidation portion of adaptive threshold
fibroticCheck = false(size(lungsbin));
if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
    fibroticCheck = imreconstruct(uint8(injurySeg.dense), ...
        uint8(completeVessels), 26);
end

% Remove pathological vessels and consolidation from complete vessels
completeVessels = imbinarize(completeVessels);
if isfield(vesselsSeg, 'pathZoneDiscarded')
    completeVessels = completeVessels & ~vesselsSeg.pathZoneDiscarded;
end
if isfield(injurySeg, 'dense')
    completeVessels = completeVessels & ~logical(injurySeg.dense);
end

% Remove near-airway region
[airwaysDist_cv, ~] = bwdist(single(airwaysSeg.airways));
outROI = airwaysDist_cv < 3 & airwaysDist_cv > 0;
completeVessels = completeVessels & ~outROI;

% Remove fissures
if any(fissures(:))
    completeVessels = completeVessels & ~fissures;
end

% Border cleanup
outBorders2_cv = bwmorph3(imfill(lungsbin, 8, 'holes'), 'remove');
lungs_R_cv = lungsProcessed.original == 1; lungs_R_cv(lungsbin == 0) = 0;
lungs_L_cv = lungsProcessed.original == 2; lungs_L_cv(lungsbin == 0) = 0;
lungs_R_bsapx_cv = lungs_R_cv & ~bwareaopen(lungs_R_cv, cv.bsapxMinVol, 8);
lungs_L_bsapx_cv = lungs_L_cv & ~bwareaopen(lungs_L_cv, cv.bsapxMinVol, 8);
lungs_bsapx_cv = lungs_R_bsapx_cv | lungs_L_bsapx_cv;
cleanBorders0_cv = imdilate(outBorders2_cv, strel('disk', cv.borderDilate0));
cleanBorders1_cv = imdilate(outBorders2_cv, strel('disk', cv.borderDilate1));
cleanBorders0_cv(lungs_bsapx_cv ~= 0) = 0;
cleanBorders1_cv(lungs_bsapx_cv == 0) = 0;
cleanBorders_cv = cleanBorders0_cv | cleanBorders1_cv;
cleanBorders_cv(lungsbin == 0) = 0;
completeVessels = completeVessels & ~cleanBorders_cv;

% Keep fibrotic check only within complete vessels, then remove from complete
fibroticCheck = logical(fibroticCheck) & completeVessels;
completeVessels = completeVessels & ~fibroticCheck;
completeVessels = bwareaopen(completeVessels, cv.minVol);

%% Vesselness on filtered volume (double threshold)
fprintf('    Applying vesselness on filtered volume...\n');

% Vesselness on filtered volume with adjusted z-spacing (vz*2)
completeVesselsFrangi0 = vesselness3D(volumes.filtered, cv.frangi.scale, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz * cv.frangi.spacingZmult], ...
    params.vessels.jerman.tau, true);

% High threshold reconstructed from mediastinal region
completeVesselsFrangi1 = completeVesselsFrangi0 > cv.frangi.highThr;
completeVesselsFrangi1 = imreconstruct((distanceMasks.close | distanceMasks.midfar), ...
    completeVesselsFrangi1, 26);

% Low threshold in non-consolidation
completeVesselsFrangi0(logical(injurySeg.dense)) = 0;
completeVesselsFrangi0 = completeVesselsFrangi0 > cv.frangi.lowThr;

completeVesselsFrangi = completeVesselsFrangi0 | completeVesselsFrangi1;
completeVesselsFrangi = completeVesselsFrangi & ~cleanBorders_cv;
if any(fissures(:))
    completeVesselsFrangi = completeVesselsFrangi & ~fissures;
end
if isfield(vesselsSeg, 'loopEndDiscarded')
    completeVesselsFrangi = completeVesselsFrangi & ~vesselsSeg.loopEndDiscarded;
end

%% Separate fibrotic check vesselness
fprintf('    Applying vesselness on fibrotic check...\n');
if isfield(injurySeg, 'dense')
    fibroticCheck(logical(injurySeg.dense)) = 0;
end
fibroticCheckFrangi = vesselness3D(fibroticCheck, cv.fibroticFrangi.scale, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    params.vessels.jerman.tau, true);
fibroticCheck = fibroticCheckFrangi > cv.fibroticFrangi.thr;
fibroticCheck = fibroticCheck & ~vessels;
fibroticCheck = bwareaopen(fibroticCheck, cv.minVolFibrotic);

%% Combine all vessel sources
completeVessels = completeVessels | fibroticCheck;
completeVessels = completeVessels | completeVesselsFrangi;
completeVessels = completeVessels | vessels;
completeVessels = completeVessels & ~brightStructuresMasked;

% Jerman filter on the combined binary to connect fragmented tracts
fprintf('    Connecting fragmented tracts with Jerman filter...\n');
completeVesselsFrangiConn = vesselness3D(completeVessels, cv.binaryFrangi.scales, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    params.vessels.jerman.tau, true);
completeVesselsFinal = completeVesselsFrangiConn > cv.binaryFrangi.thr;

vessels = completeVesselsFinal | vessels;
vessels = bwareaopen(vessels, cv.volThr);

%% FA check #1: Small objects with large diameter and low anisotropy
fprintf('    FA check on small objects after integration...\n');
cleanFA1 = bwareaopen(vessels, cv.faCheck1.maxVol);
checkFA1 = vessels & ~cleanFA1;
diamFA1 = 2 * bwdist(~checkFA1);
skelFA1 = bwskel(checkFA1);
bpFA1 = bwmorph3(skelFA1, 'branchpoints');
skelFA1(bpFA1 ~= 0) = 0;
diamSkelFA1 = diamFA1 .* double(skelFA1);
seedFA1 = diamSkelFA1 > cv.faCheck1.minDiam;

checkLblFA1 = bwlabeln(seedFA1);
checkEigFA1 = regionprops3(checkLblFA1, 'EigenValues');
idxKeepFA1 = false(1, size(checkEigFA1, 1));
for i = 1:size(checkEigFA1, 1)
    eigenvalues = checkEigFA1.EigenValues{i};
    FAval = computeFractionalAnisotropy(eigenvalues);
    if FAval > cv.faCheck1.FA
        idxKeepFA1(i) = true;
    end
end
discardFA1 = seedFA1 & ~ismember(checkLblFA1, find(idxKeepFA1));
discardFA1 = imreconstruct(discardFA1, checkFA1, 26);
vessels = vessels & ~discardFA1;

%% FA check #2: Border+far objects with low anisotropy (inverted)
fprintf('    FA check on border objects (inverted)...\n');
checkFA2 = vessels;
checkFA2_small = checkFA2 & ~bwareaopen(checkFA2, cv.faCheck2.maxVol);
checkFA2(cleanBorders_cv == 0) = 0;
checkFA2(distanceMasks.far == 0) = 0;
checkFA2 = bwareaopen(checkFA2, cv.faCheck2.minVol);
checkFA2 = imreconstruct(checkFA2, checkFA2_small, 26);
discardFA2 = checkFA2;

checkLblFA2 = bwlabeln(discardFA2);
checkEigFA2 = regionprops3(checkLblFA2, 'EigenValues');
idxKeepFA2 = false(1, size(checkEigFA2, 1));
for i = 1:size(checkEigFA2, 1)
    eigenvalues = checkEigFA2.EigenValues{i};
    FAval = computeFractionalAnisotropy(eigenvalues);
    if FAval < cv.faCheck2.FA
        idxKeepFA2(i) = true;
    end
end
discardFA2 = ismember(checkLblFA2, find(idxKeepFA2));
vessels = vessels & ~discardFA2;

%% Pre-wall reassignment
fprintf('    Pre-wall reassignment of small consolidation fragments...\n');
pw = params.vessels.preWall;

% Small consolidation fragments (< 30 voxels)
wallFragments = logical(injurySeg.dense) & ...
    ~bwareaopen(logical(injurySeg.dense), pw.consolMaxVol, 8);

% Distance from airways and vessels
[airwaysDist, ~] = bwdist(single(airwaysSeg.airways));
[vesselsDist, ~] = bwdist(single(vessels));

% Identify wall fragments near each structure
wallsNearAirways = wallFragments & (airwaysDist < pw.airwaysDist & airwaysDist > 0);
wallsNearAirwaysClose = (volumes.filtered < pw.airwayWallHU) & ...
    (airwaysDist <= pw.airwaysCloseDist & airwaysDist > 0);
wallsNearAirways = wallsNearAirways | wallsNearAirwaysClose;
wallsNearVessels = wallFragments & (vesselsDist < pw.vesselsDist & vesselsDist > 0);

% Reassign: small consol near vessels → add to vessels
vessels = (vessels | wallsNearVessels) & ~wallsNearAirways;

%% Fissure cleaning
fprintf('    Fissure cleaning...\n');
vessels = bwareaopen(vessels, 16, 26);

if any(fissures(:))
    vesselsFissure = vessels & ~bwareaopen(vessels, lc.fissureMaxVol, 26);
    vesselsFissure(fissures == 0) = 0;
    vessels = vessels & ~vesselsFissure;
    vessels(lungsbin == 0) = 0;
end

%% Last corrections pt1: Progressive leakage detection
fprintf('    Progressive leakage detection (5 passes)...\n');

% Compute OutBorders = far region on external boundary
outBorders2 = bwmorph3(imfill(lungsbin, 8, 'holes'), 'remove');
outBorders = distanceMasks.far & outBorders2;

% Compute base/apex mask
lungs_R = lungsProcessed.original == 1; lungs_R(lungsbin == 0) = 0;
lungs_L = lungsProcessed.original == 2; lungs_L(lungsbin == 0) = 0;
lungs_R_bsapx = lungs_R & ~bwareaopen(lungs_R, lc.bsapxMinVol, 8);
lungs_L_bsapx = lungs_L & ~bwareaopen(lungs_L, lc.bsapxMinVol, 8);
lungsBsApx = lungs_R_bsapx | lungs_L_bsapx;

% Dilated distal region (excluding base/apex)
outBordersDilate = imdilate(outBorders, strel('disk', lc.outBordersDilate));
outBordersDilate(lungsBsApx ~= 0) = 0;
outBordersDilate(lungsbin == 0) = 0;

% Consolidation in distal region
outBordersDilateConsol = logical(injurySeg.dense);
outBordersDilateConsol(outBordersDilate == 0) = 0;

% Inner region (complement of dilated distal)
inBorders = lungsbin & ~outBordersDilate;

% Diameter skeleton for leakage detection
diameter0 = 2 * bwdist(~vessels);
skel0 = bwskel(vessels, 'MinBranchLength', lc.skelMinBranch);
diamImg0 = diameter0 .* double(skel0);

vesselsWork = vessels;

% --- Pass 0: Small disconnected objects (<1000) with diameter>3 in distal ---
small0 = vesselsWork & ~bwareaopen(vesselsWork, lc.pass0.maxVol);
leak0 = diamImg0;
leak0(small0 == 0) = 0;
leak0(outBordersDilate == 0) = 0;
leak0 = leak0 > lc.pass0.minDiam;
leak0 = bwareaopen(leak0, lc.pass0.minArea);
leakReco0 = imreconstruct(leak0, small0, 26);

% --- Pass 1: Very small objects (<300) with diameter>3 in distal ---
small1 = vesselsWork & ~bwareaopen(vesselsWork, lc.pass1.maxVol);
leak1 = diamImg0;
leak1(small1 == 0) = 0;
leak1(outBordersDilate == 0) = 0;
leak1 = leak1 > lc.pass1.minDiam;
leakReco1 = imreconstruct(leak1, leakReco0, 26);

% --- Pass 2: Leakage in consolidation's distal region, diameter>2.5 ---
leak2 = diamImg0;
leak2(outBordersDilateConsol == 0) = 0;
leak2 = leak2 > lc.pass2.minDiam;
leak2 = bwareaopen(leak2, lc.pass2.minArea);
leakReco2 = imreconstruct(leak2, vesselsWork, 8);

% --- Pass 3: General distal leakage, diameter>4 ---
leak3 = diamImg0;
leak3(outBordersDilate == 0) = 0;
leak3 = leak3 > lc.pass3.minDiam;
leak3 = bwareaopen(leak3, lc.pass3.minArea);
leakReco3 = imreconstruct(leak3, vesselsWork, 8);

% --- Pass 4: Diameter>3, disconnected from inner region ---
leak4 = diamImg0 > lc.pass4.minDiam;
leakReco4 = imreconstruct(leak4, vesselsWork, 26);
keep4 = imreconstruct(inBorders, leakReco4, 26);
remove4 = leakReco4 & ~keep4;
leakReco4 = imreconstruct(remove4, leakReco4, 26);

% Combine all leakage passes
allLeakage = leakReco0 | leakReco1 | leakReco2 | leakReco3 | leakReco4;
outBordersLastCorr = leakReco2 | leakReco3;  % Used for seed-based cleanup

vessels = vessels & ~allLeakage;

%% Last corrections pt2: Eigenvector orientation filter
fprintf('    Eigenvector orientation filter...\n');

outsideBsApx = lungsbin & ~lungsBsApx;
vesselsCheck = vessels & ~bwareaopen(vessels, lc.eigvec.maxVol);
vesselsCheckInner = imreconstruct(outsideBsApx, vesselsCheck, 26);
vesselsCheck = vesselsCheck & ~vesselsCheckInner;  % Only base/apex objects

vesselsCheckLbl = bwlabeln(vesselsCheck);
checkEig = regionprops3(vesselsCheckLbl, 'EigenValues', 'EigenVectors');
idxDiscard = false(1, size(checkEig, 1));

for i = 1:size(checkEig, 1)
    eigenvalues = checkEig.EigenValues{i};
    eigenvectors = checkEig.EigenVectors{i};
    FAval = (1/sqrt(2)) * sqrt((eigenvalues(1)-eigenvalues(2))^2 ...
        + (eigenvalues(2)-eigenvalues(3))^2 ...
        + (eigenvalues(1)-eigenvalues(3))^2) / ...
        sqrt(eigenvalues(1)^2 + eigenvalues(2)^2 + eigenvalues(3)^2);
    if FAval > lc.eigvec.FA
        if abs(eigenvectors(3,1)) < abs(eigenvectors(2,1)) && ...
                abs(eigenvectors(3,1)) < abs(eigenvectors(1,1))
            idxDiscard(i) = true;
        end
    end
end

vesselsDiscard = ismember(vesselsCheckLbl, find(idxDiscard));
vesselsDiscard = imreconstruct(vesselsDiscard, vesselsCheck, 26);
vessels = vessels & ~vesselsDiscard;

%% Last corrections pt3: Seed-based cleanup
fprintf('    Seed-based cleanup of residual leakages...\n');

% D1: Small objects connected to previous leakage corrections → remove
cleanMain = bwareaopen(vessels, lc.seed.maxVol);
checkSmall = vessels & ~cleanMain;
seedDiscard = imreconstruct(imdilate(outBordersLastCorr, strel('sphere', 1)), ...
    checkSmall, 26);
vessels = vessels & ~seedDiscard;

% D2: Remaining small objects (<150): diameter+FA check
cleanMain2 = bwareaopen(vessels, lc.final.maxVol);
checkSmall2 = vessels & ~cleanMain2;
diam2 = 2 * bwdist(~checkSmall2);
skel2 = bwskel(checkSmall2);
bp2 = bwmorph3(skel2, 'branchpoints');
skel2(bp2 ~= 0) = 0;
diamSkel2 = diam2 .* double(skel2);
largeSkel2 = diamSkel2 > lc.final.minDiam;

checkLbl2 = bwlabeln(largeSkel2);
checkEig2 = regionprops3(checkLbl2, 'EigenValues');
idxKeep2 = false(1, size(checkEig2, 1));

for i = 1:size(checkEig2, 1)
    eigenvalues = checkEig2.EigenValues{i};
    FAval = (1/sqrt(2)) * sqrt((eigenvalues(1)-eigenvalues(2))^2 ...
        + (eigenvalues(2)-eigenvalues(3))^2 ...
        + (eigenvalues(1)-eigenvalues(3))^2) / ...
        sqrt(eigenvalues(1)^2 + eigenvalues(2)^2 + eigenvalues(3)^2);
    if FAval > lc.final.FA
        idxKeep2(i) = true;
    end
end

checkKeep2 = ismember(checkLbl2, find(idxKeep2));
finalDiscard = largeSkel2 & ~checkKeep2;
finalDiscard = imreconstruct(finalDiscard, checkSmall2, 26);
vessels = vessels & ~finalDiscard;

%% Large consolidation reclassification
fprintf('    Reclassifying large consolidations misidentified as vessels...\n');
lcp = params.vessels.largeConsol;

% Detect dense tissue (HU > -150)
largeConsol = volumes.ct > lcp.minHU;
largeConsol = bwareaopen(largeConsol, lcp.initMinVol, 8);
largeConsol(lungsbin == 0) = 0;

% Two paths: distal consolidation region + general far region
largeConsol2 = largeConsol & outBordersDilateConsol;
largeConsol(distanceMasks.far == 0) = 0;
largeConsol = largeConsol | largeConsol2;
largeConsol = bwareaopen(largeConsol, lcp.mainMinVol, 26);
largeConsol2 = bwareaopen(largeConsol2, lcp.distalMinVol, 26);

% Border reconstruction: find consolidations touching pleural surface
cleanBorders = bwmorph3(lungsbin, 'remove');
largeConsol1 = imreconstruct(cleanBorders, largeConsol, 8);
largeConsol1 = largeConsol1 | imreconstruct(largeConsol1, vessels, 8);

% Identify vessels that overlap with border-connected consolidation
vesselsConsolidated = vessels & largeConsol1;
vesselsConsolidated = bwareaopen(vesselsConsolidated, lcp.vesselMinVol, 8);

% Combine all consolidation candidates and expand into vessels
addConsolidation = vesselsConsolidated | largeConsol2;
addConsolidation = addConsolidation | imreconstruct(addConsolidation, vessels, 8);

% Remove consolidation from vessels
vessels = vessels & ~addConsolidation;

% Re-filter remaining vessels through vesselness to recover true tubular structures
vesselsRefiltered = vesselness3D(vessels, lcp.vesselness.scales, ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    lcp.vesselness.tau, true);
vessels = vesselsRefiltered > lcp.vesselness.threshold;

% Small post-vesselness fragments connected to consolidation → reclassify
largeConsol2post = vessels & ~bwareaopen(vessels, lcp.postMinVol);
largeConsol2post = imreconstruct(largeConsol1, largeConsol2post, 26);
vessels = vessels & ~largeConsol2post;

% Finalize consolidation output
addConsolidation = addConsolidation | largeConsol2post;
addConsolidation = addConsolidation & ~vessels;

% Also add the dense tissue mask itself (minus final vessels)
largeConsol0 = largeConsol & ~vessels;
addConsolidation = addConsolidation | largeConsol0;

%% Add large vessels that may have been missed by the Jerman filter
fprintf('    Segmenting large vessels from intensity...\n');
[largeVessels, add2Final] = addLargeVessels(volumes, lungsProcessed, vessels, ...
    airwaysSeg, injurySeg, distanceMasks, params);

% Add2Final goes into the main vessel segmentation
vessels = vessels | add2Final;

%% Last removal of false positives
fprintf('    Last removal of false positives...\n');
lfp = params.vessels.lastFP;

% --- Pass E1: FA check on small objects (<500) with large diameter ---
checkE1 = vessels & ~bwareaopen(vessels, lfp.pass1.maxVol);
diamE1 = 2 * bwdist(~checkE1);
skelE1 = bwskel(checkE1);
diamSkelE1 = diamE1 .* double(skelE1);
seedE1 = diamSkelE1 > lfp.pass1.minDiam;
seedE1 = bwareaopen(seedE1, lfp.pass1.minArea);
checkE1 = imreconstruct(seedE1, vessels, 26);

checkLblE1 = bwlabeln(checkE1);
checkEigE1 = regionprops3(checkLblE1, 'EigenValues');
idxKeepE1 = false(1, size(checkEigE1, 1));
for i = 1:size(checkEigE1, 1)
    eigenvalues = checkEigE1.EigenValues{i};
    FAval = (1/sqrt(2)) * sqrt((eigenvalues(1)-eigenvalues(2))^2 ...
        + (eigenvalues(2)-eigenvalues(3))^2 ...
        + (eigenvalues(1)-eigenvalues(3))^2) / ...
        sqrt(eigenvalues(1)^2 + eigenvalues(2)^2 + eigenvalues(3)^2);
    if FAval > lfp.pass1.FA
        idxKeepE1(i) = true;
    end
end
discardE1 = checkE1 & ~ismember(checkLblE1, find(idxKeepE1));

% --- Pass E2: Volume+diameter check on very small objects (<70) ---
checkE2 = vessels & ~bwareaopen(vessels, lfp.pass2.maxVol);
diamE2 = 2 * bwdist(~checkE2);
skelE2 = bwskel(checkE2);
diamSkelE2 = diamE2 .* double(skelE2);
seedE2 = diamSkelE2 > lfp.pass2.minDiam;
discardE2 = imreconstruct(seedE2, vessels, 6);

discardE = discardE1 | discardE2;
vessels = vessels & ~discardE;
addConsolidation = addConsolidation | discardE;

% --- Pass E3: Distally segregated objects with large diameter ---
% InBorders2 = wider inner region (dilation 12 instead of 9)
outBordersDilate2 = imdilate(outBorders2, strel('disk', lc.outBorders2Dilate));
outBordersDilate2(lungsbin == 0) = 0;
outBordersDilate2(lungsBsApx ~= 0) = 0;
inBorders2 = lungsbin & ~(outBordersDilate2 | lungsBsApx);

lastSeg = imreconstruct(inBorders2, vessels, 26);
lastSeg = vessels & ~lastSeg;  % Objects NOT connected to inner region
diamLS = 2 * bwdist(~lastSeg);
skelLS = bwskel(lastSeg);
diamSkelLS = diamLS .* double(skelLS);
seedLS = diamSkelLS >= lfp.distal.minDiam;
seedLS = bwareaopen(seedLS, lfp.distal.minArea, 26);
lastSeg = imreconstruct(seedLS, lastSeg, 26);

% Remove objects in base/apex
lastSegRemove = imreconstruct(lungsBsApx, lastSeg, 26);
lastSeg = lastSeg & ~lastSegRemove;

% Remove objects connected to consolidation
lastSeg = imreconstruct(logical(injurySeg.dense) | addConsolidation, lastSeg, 26);
vessels = vessels & ~lastSeg;

%% Final size filtering
vessels = bwareaopen(vessels, params.vessels.cleaning.finalMinVolume, 6);

%% Output
vesselsRefined = vessels;

fprintf('    Vessel post-processing complete\n');

end


%% Helper: Adaptive vessel thresholding in healthy tissue
function completeVessels = adaptiveVesselThreshold(volumeFiltered, masks, ~)
    C = 180;
    ksize = 15;
    sigma = 0.3 * ((ksize - 1) * 0.5 - 1) + 0.8;

    ROI = masks.lungsBinary;
    completeVessels = uint8(zeros(size(masks.lungsBinary)));

    for k = 1:size(volumeFiltered, 3)
        ROIslice = double(ROI(:, :, k));
        src = double(imcomplement(volumeFiltered(:, :, k))) .* (ROIslice > 0);
        roi = double(ones(size(ROIslice)));
        roi = roi .* (ROIslice > 0);

        blur = imgaussfilt(src, sigma, 'FilterSize', ksize);
        roi = imgaussfilt(roi, sigma, 'FilterSize', ksize);
        smooth = blur ./ roi;
        output = 255 .* (src <= (smooth - C));
        completeVessels(:, :, k) = uint8(output .* (ROIslice > 0));
    end
end

%% Helper: Recover large vessel lumen from HU-based solid-tissue candidates
function [largeVessels, add2Final] = addLargeVessels(volumes, lungsProcessed, ...
    vesselsUpdate, airwaysSeg, injurySeg, distanceMasks, params)
%
% Builds a solid-tissue candidate pool from HU thresholding (HU > huMin1),
% hole-filling, and pruning (HU > huMin2), then constrains it by distance
% from Jerman skeleton-derived large-diameter seeds. This compensates for
% the attenuated Jerman response at the centre of very large vessels and
% for any gaps in the main filter-based segmentation.
%
% Two-level seeding (seedLow / seedHigh) anchors the reconstruction to
% confirmed large-calibre regions and prevents expansion into small vessels.
% Three additional sources supplement the core:
%   - 4-conn expansion: contiguous lumen beyond the distance gate, filtered
%     against large consolidations to limit leakage.
%   - Bronchial seed (largeVesselsAdd0): pulmonary arteries run alongside
%     large bronchi (bronchovascular bundle). A corona around large-calibre
%     airways (diameter > bronchDiamThr, distance < bronchCrownDist) seeds
%     an 8-conn growth in the HU pool, recovering vessel lumen that the
%     Jerman skeleton seeds alone may miss.
%   - LargeVessels1 wall-zone addition: a permissive intermediate
%     segmentation (seedLow only, no distance gate) is computed and its
%     intersection with a narrow corona around the already-confirmed
%     largeVessels is added back. This recovers partial-volume lumen at the
%     vessel boundary that the stricter distance gate excluded.
%
% Outputs:
%   largeVessels - Hilar/central large vessel lumen for wall computation in
%                  finalizeSegmentations. NOT added to the main segmentation.
%   add2Final    - Mid-to-large vessel lumen added to the main segmentation.
%                  Anchored directly to the mediastinal region (close mask).

    lv = params.vessels.largeVessels;
    lungsbin = lungsProcessed.binary;

    %% Diameter seeds from Jerman-refined vessels
    % seedLow : diam > diamSeedLow (3.5 mm), retained only where connected to seedHigh
    % seedHigh: diam > diamSeedHigh (4.5 mm) — used to constrain seedLow
    diamVess    = 2 * bwdist(~vesselsUpdate);
    skelVess    = bwskel(vesselsUpdate);
    diamImgVess = diamVess .* double(skelVess);

    seedLow  = imdilate(diamImgVess > lv.diamSeedLow,  strel('sphere', 1));
    seedHigh = imdilate(diamImgVess > lv.diamSeedHigh, strel('sphere', 1));
    seedLow  = imreconstruct(seedHigh, seedLow, 26);  % keep only parts connected to high-diam
    clear seedHigh;
    seedLow  = bwareaopen(seedLow, lv.diamSeedMinVol, 6);

    %% HU-based solid-tissue pool
    huLumen0 = volumes.ct > lv.huMin1;
    huLumen0 = huLumen0 & lungsbin;
    huLumen0 = imfill(huLumen0, 8, 'holes');
    huLumen0(volumes.ct < lv.huMin2) = 0;

    %% LargeVessels1: permissive intermediate (seedLow only, no distance gate)
    % Grows from seedLow into the full HU pool without the distance constraint
    % applied to the main largeVessels. Used later to recover partial-volume
    % lumen within the wall zone of the confirmed vessel region.
    largeVessels1 = imreconstruct(seedLow, huLumen0, 8);
    largeVessels1 = imreconstruct(distanceMasks.close, largeVessels1, 26);

    %% Distance gate from seedLow
    [seedLowDist, ~] = bwdist(single(seedLow));
    huGate_large = seedLowDist <= lv.huSeedDistLarge & seedLowDist > 0;
    huGate_add   = seedLowDist <  lv.huSeedDistAdd   & seedLowDist > 0;

    %% largeVessels: grow from seedLow within the 3.5 mm-gated HU pool,
    %  then keep only components connected to the mediastinal region.
    %  8-connectivity (consistent with ORIGINAL) limits diagonal growth,
    %  keeping the core region topologically close to the confirmed seeds.
    huCand_large = huLumen0 & (huGate_large | seedLow);
    largeVessels = imreconstruct(seedLow, huCand_large, 8);
    largeVessels = imreconstruct(distanceMasks.close, largeVessels, 26);

    %% Narrow 4-conn expansion into the full HU pool + consolidation filtering
    % The expansion captures contiguous lumen beyond the distance gate while
    % the 4-connectivity keeps it topologically close to the confirmed region.
    % Only the expanded portion (not the core) is filtered against consolidation.
    huLargeExp     = imreconstruct(largeVessels, huLumen0, 4);
    huLargeExpOnly = huLargeExp & ~largeVessels;

    if any(injurySeg.dense(:))
        consolBig      = bwareaopen(logical(injurySeg.dense), lv.huConsolMinVol, 6);
        huLargeRem     = imreconstruct(consolBig, huLargeExpOnly, 6);
        huLargeExpOnly = huLargeExpOnly & ~huLargeRem;
    end
    largeVessels = largeVessels | huLargeExpOnly;

    %% Bronchial seed (largeVesselsAdd0)
    % Pulmonary arteries travel alongside large bronchi (bronchovascular bundle).
    % A corona around large-calibre airways seeds an 8-conn growth in the HU
    % pool, recovering vessel lumen adjacent to the bronchial tree that is not
    % reachable from the Jerman diameter seeds alone.
    diamAw    = 2 * bwdist(~logical(airwaysSeg.airways));
    skelAw    = bwskel(logical(airwaysSeg.airways));
    diamImgAw = diamAw .* double(skelAw);

    largeAw = diamImgAw > lv.bronchDiamThr;
    largeAw = imreconstruct(largeAw, logical(airwaysSeg.airways), 8);
    largeAw = imreconstruct(distanceMasks.close, largeAw, 26);

    [awCrownDist, ~] = bwdist(single(largeAw));
    bronchialCrown = awCrownDist < lv.bronchCrownDist & awCrownDist > 0;
    largeVesselsAdd0 = imreconstruct(bronchialCrown, huLumen0, 8);
    largeVessels = largeVessels | largeVesselsAdd0;
    clear diamAw skelAw diamImgAw largeAw awCrownDist bronchialCrown largeVesselsAdd0;

    %% LargeVessels1 wall-zone addition
    % The intersection of the permissive intermediate with a narrow corona
    % around the confirmed largeVessels recovers partial-volume lumen at the
    % vessel boundary that the stricter distance gate had excluded.
    [lv1Dist, ~] = bwdist(single(largeVessels));
    lv1WallZone = lv1Dist < lv.lv1WallDist & lv1Dist > 0;
    largeVessels = largeVessels | (largeVessels1 & lv1WallZone);
    clear lv1Dist lv1WallZone largeVessels1;

    %% add2Final: HU pool within 3 mm of seedLow, anchored directly to mediastinum
    % More conservative than largeVessels (tighter gate, no seed-grow step).
    huCand_add = huLumen0 & (huGate_add | seedLow);
    add2Final  = imreconstruct(distanceMasks.close, huCand_add, 26);
    add2Final  = add2Final & lungsbin;

    %% Airway exclusion (approximate wall proxy; precise exclusion in finalizeSegmentations)
    [awDist, ~] = bwdist(single(airwaysSeg.airways));
    awWallApprox = awDist < lv.awWallDist & awDist > 0;
    largeVessels = largeVessels & ~awWallApprox & ~airwaysSeg.airways;
    add2Final    = add2Final    & ~airwaysSeg.airways;

    %% Far region and size cleanup
    largeVessels = largeVessels & ~distanceMasks.far;
    largeVessels = bwareaopen(largeVessels, 4,  4);   % remove isolated voxels (4-conn)
    largeVessels = bwareaopen(largeVessels, lv.minVol, 6);

    %% Mutual exclusivity: largeVessels holds only what is NOT in add2Final or Jerman vessels
    largeVessels = largeVessels & ~add2Final;
    largeVessels = largeVessels & ~vesselsUpdate;
end
