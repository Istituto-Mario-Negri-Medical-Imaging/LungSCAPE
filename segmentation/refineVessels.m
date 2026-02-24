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
%   7. Large vessel reconstruction for structures missed by Jerman filter
%
%   Inputs:
%       vesselsSeg     - Structure with .final, .vesselness fields
%       volumes        - Structure with .filtered, .ct fields
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

%% Adaptive thresholding for complete vessels in healthy tissue
fprintf('    Applying adaptive thresholding in healthy tissue...\n');

masks.lungs = lungsProcessed.original;
masks.lungsBinary = lungsbin;
masks.airways = airwaysSeg.airways;
masks.far = distanceMasks.far;

completeVessels = adaptiveVesselThreshold(volumes.filtered, masks, params);

% Vesselness filtering on adaptive threshold result
fprintf('    Applying vesselness filter to complete vessels...\n');
completeVesselsFrangi = vesselness3D(completeVessels, ...
    [0.5, 1], ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    0.5, true);

completeVessels = completeVesselsFrangi > params.vessels.threshold.complete;

% Check for fibrotic bands (near consolidation)
if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
    fibroticCheck = imreconstruct(uint8(injurySeg.dense), ...
        uint8(completeVessels), 26);
    fibroticCheck(vesselsSeg.vesselness > 0.8) = 0;
    completeVessels = completeVessels & ~fibroticCheck;
end

% Combine with vessels from segmentVessels
vessels = vessels | completeVessels;

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
largeVessels = addLargeVessels(volumes, lungsProcessed, vessels, ...
    airwaysSeg, distanceMasks, params);

% Combine large vessels with refined vessels
vessels = vessels | largeVessels;

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
function completeVessels = adaptiveVesselThreshold(volumeFiltered, masks, params)
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

    % Clean up borders
    if isfield(masks, 'far')
        lungs_R = masks.lungs == 1;
        lungs_L = masks.lungs == 2;

        lungs_R_bsapx = lungs_R & ~bwareaopen(lungs_R, 7000, 8);
        lungs_L_bsapx = lungs_L & ~bwareaopen(lungs_L, 7000, 8);
        lungs_bsapx = lungs_R_bsapx | lungs_L_bsapx;

        cleanBorders0 = imdilate(masks.far, strel('disk', 5));
        cleanBorders1 = imdilate(masks.far, strel('disk', 2));
        cleanBorders0(lungs_bsapx ~= 0) = 0;
        cleanBorders1(lungs_bsapx == 0) = 0;
        cleanBorders = cleanBorders0 | cleanBorders1;
        cleanBorders(masks.lungs == 0) = 0;

        completeVessels(cleanBorders ~= 0) = 0;
    end

    completeVessels = logical(completeVessels);
end

%% Helper: Add large vessels
function largeVessels = addLargeVessels(volumes, lungsProcessed, vesselsUpdate, ...
    airwaysSeg, distanceMasks, params)

    largeVessels0 = volumes.filtered > params.vessels.largeVessels.MinHU & ...
                    volumes.filtered < params.vessels.largeVessels.MaxHU;
    largeVessels0(lungsProcessed.binary == 0) = 0;
    largeVessels0 = largeVessels0 & ~airwaysSeg.airways;

    [diameterImage, ~] = computeVesselDiameter(largeVessels0, 0);

    diameterImageLarge = diameterImage > params.vessels.largeVessels.diameter;
    diameterImageLarge = bwareaopen(diameterImageLarge, ...
        params.vessels.largeVessels.minLength, 26);

    [largeDist, ~] = bwdist(single(diameterImageLarge));
    largeDistMask = largeDist < 3 & largeDist > 0;
    largeVessels = largeVessels0 & (largeDistMask | diameterImageLarge);
    largeVessels = imreconstruct(diameterImageLarge, largeVessels, 8);
    largeVessels = imreconstruct(distanceMasks.close, largeVessels, 26);

    largeVessels = largeVessels & ~airwaysSeg.airways;
    largeVessels = bwareaopen(largeVessels, params.vessels.largeVessels.minVol, 6);
    largeVessels = largeVessels & ~vesselsUpdate;
end
