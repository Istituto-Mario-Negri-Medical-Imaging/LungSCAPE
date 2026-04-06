function vesselsClassified = classifyVesselsByDiameter(vessels, params)
%CLASSIFYVESSELSBYDIAMETER Classify vessels into small/medium/large categories
%
%   vesselsClassified = classifyVesselsByDiameter(vessels, params)
%
%   Classifies vessels into three diameter categories using skeleton-based
%   diameter estimation in mm and nearest-neighbor propagation for
%   non-skeleton voxels. Uses parallel processing for efficiency.
%
%   Classification (BV5/BV10 framework, Estépar et al. AJRCCM 2013):
%       Label 1: Large vessels  — BV10  (diameter > 3.6 mm, CSA > 10 mm²)
%       Label 2: Medium vessels — BV5-10 (2.5 mm < diameter <= 3.6 mm)
%       Label 3: Small vessels  — BV5   (diameter <= 2.5 mm, CSA <= 5 mm²)
%
%   Diameters are measured in mm using the mean in-plane voxel spacing.
%
%   Algorithm:
%       1. Compute diameter map in mm from distance transform × voxel spacing
%       2. Extract skeleton and measure diameter at skeleton points
%       3. Classify skeleton points by diameter
%       4. Handle boundary cases (reclassify small fragments)
%       5. Propagate labels to non-skeleton voxels using nearest neighbor
%
%   Inputs:
%       vessels - Binary vessel segmentation (logical)
%       params  - Processing parameters (see getProcessingParameters)
%
%   Outputs:
%       vesselsClassified - Labeled volume (uint8):
%           0 = background
%           1 = large vessels  (BV10,  > 3.6 mm)
%           2 = medium vessels (BV5-10, 2.5–3.6 mm)
%           3 = small vessels  (BV5,   <= 2.5 mm)
%
%   See also: computeVesselDiameter, propagateLabelsNN

fprintf('  Classifying vessels by diameter...\n');

% Diameter thresholds (mm)
dp = params.vessels.diameterClassification;
smallThr = dp.smallDiameterThreshold;   % mm
largeThr = dp.largeDiameterThreshold;   % mm

%% Compute diameter map in mm
fprintf('    Computing vessel diameter map (mm)...\n');

% bwdist returns distances in voxel units. Convert to mm by scaling with
% the mean in-plane spacing. On the skeleton the dominant cross-section is
% transverse, so in-plane spacing is the appropriate metric. For nearly
% isotropic in-plane voxels (px ≈ py) this is accurate; the z-spacing
% contribution is negligible because skeleton points rarely have their
% nearest border along z.
voxelSpacing_mm = mean([params.voxel.px, params.voxel.py]);
diameter = 2 * bwdist(~vessels) * voxelSpacing_mm;   % mm
skeletonImage = bwskel(vessels);
diameterImage = diameter .* double(skeletonImage);

fprintf('    Voxel spacing used for diameter: %.3f mm (mean in-plane)\n', voxelSpacing_mm);

%% Initial classification by diameter
fprintf('    Classifying skeleton points by diameter...\n');

diameterImageSmall = diameterImage <= smallThr & diameterImage > 0;
diameterImageMid   = diameterImage <= largeThr & diameterImage > smallThr;
diameterImageLarge = diameterImage > largeThr;

%% Reclassify medium-sized fragments (span of tolerance)
fprintf('    Reclassifying medium-sized fragments...\n');

% Medium fragments < 27 voxels may need reclassification
mid2change = diameterImageMid & ~bwareaopen(diameterImageMid, 27, 26);
diameterImageMid = diameterImageMid & ~mid2change;

% Create label image at skeleton endpoints
% This will be used to propagate labels to reclassified regions
diameterImageLbl = uint8(zeros(size(vessels)));
diameterImageLbl(diameterImageSmall > 0) = 3;
diameterImageLbl(diameterImageMid > 0) = 0;  % Will be filled
diameterImageLbl(diameterImageLarge > 0) = 1;

% Extract endpoints only (for nearest neighbor search)
diameterImageLbl = diameterImageLbl .* uint8(bwmorph3(...
    imbinarize(diameterImageLbl, 0), 'endpoints'));

% Propagate labels to medium-sized fragments
[diameterImageMid, diameterImageLarge, diameterImageSmall] = ...
    propagateLabelsToFragments(mid2change, diameterImageLbl, ...
    diameterImageMid, diameterImageLarge, diameterImageSmall, params);

%% Reclassify large-sized fragments (span of tolerance)
fprintf('    Reclassifying large-sized fragments...\n');

% Large fragments < 17 voxels may need reclassification
large2change = diameterImageLarge & ~bwareaopen(diameterImageLarge, 17, 26);
diameterImageLarge = diameterImageLarge & ~large2change;

% Create label image at skeleton endpoints
diameterImageLbl = uint8(zeros(size(vessels)));
diameterImageLbl(diameterImageSmall > 0) = 0;  % Will be filled
diameterImageLbl(diameterImageMid > 0) = 2;
diameterImageLbl(diameterImageLarge > 0) = 0;  % Will be filled

% Extract endpoints only
diameterImageLbl = diameterImageLbl .* uint8(bwmorph3(...
    imbinarize(diameterImageLbl, 0), 'endpoints'));

% Propagate labels to large-sized fragments
[diameterImageMid, ~, ~] = propagateLabelsToFragments(large2change, ...
    diameterImageLbl, diameterImageMid, diameterImageLarge, ...
    diameterImageSmall, params);

%% Create final skeleton label map
fprintf('    Creating final skeleton classification...\n');

diameterImageLbl = uint8(zeros(size(vessels)));
diameterImageLbl(diameterImageSmall > 0) = 3;
diameterImageLbl(diameterImageMid > 0) = 2;
diameterImageLbl(diameterImageLarge > 0) = 1;

%% Propagate labels to all vessel voxels (not just skeleton)
fprintf('    Propagating labels to all vessel voxels...\n');

% Identify non-skeleton vessel voxels
classify = vessels;
classify(skeletonImage ~= 0) = 0;

% Use nearest neighbor to assign labels
vesselsClassified = propagateLabelsNN(classify, diameterImageLbl, params);

% Add skeleton labels
vesselsClassified = vesselsClassified + diameterImageLbl;

fprintf('    Vessel classification complete\n');

end


%% Helper: Propagate labels to fragments
function [mid, large, small] = propagateLabelsToFragments(fragments, labelImg, ...
    mid, large, small, params)
    % Propagate labels from endpoints to unclassified fragments

    % Find nearest labeled endpoint for each fragment voxel
    [~, lbl_dist] = bwdist(imbinarize(labelImg, 0));
    lbl_dist = lbl_dist .* uint32(fragments);

    % Get unique fragment indices
    lbl_idx = unique(lbl_dist);
    lbl_idx(1) = [];  % Remove background
    lbl_idx_new = zeros(size(lbl_idx));

    % For each fragment, find nearest endpoint and get its label
    parfor (lbl_idx_k = 1:length(lbl_idx), params.parallel.workers)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(fragments), val2change);
        lbl_idx_new(lbl_idx_k) = labelImg(val_row, val_col, val_z);
    end

    % Assign fragments to appropriate size category
    fill3 = lbl_idx(lbl_idx_new == 3);
    fill3 = (ismember(lbl_dist, fill3)) .* 3;
    fill2 = lbl_idx(lbl_idx_new == 2);
    fill2 = (ismember(lbl_dist, fill2)) .* 2;
    fill1 = lbl_idx(lbl_idx_new == 1);
    fill1 = ismember(lbl_dist, fill1);

    large = large | fill1;
    mid = mid | fill2;
    small = small | fill3;
end


%% Helper: Propagate labels using nearest neighbor
function classifiedImg = propagateLabelsNN(toClassify, labelImg, params)
    % Propagate labels to unlabeled voxels using nearest neighbor

    [~, lbl_dist] = bwdist(imbinarize(labelImg, 0));
    lbl_dist = lbl_dist .* uint32(toClassify);

    % Get unique indices
    lbl_idx = unique(lbl_dist);
    lbl_idx(1) = [];  % Remove background
    lbl_idx_new = zeros(size(lbl_idx));

    % For each voxel, find nearest labeled voxel
    parfor (lbl_idx_k = 1:length(lbl_idx), params.parallel.workers)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(toClassify), val2change);
        lbl_idx_new(lbl_idx_k) = labelImg(val_row, val_col, val_z);
    end

    % Create classified image
    fill3 = lbl_idx(lbl_idx_new == 3);
    fill3 = (ismember(lbl_dist, fill3)) .* 3;
    fill2 = lbl_idx(lbl_idx_new == 2);
    fill2 = (ismember(lbl_dist, fill2)) .* 2;
    fill1 = lbl_idx(lbl_idx_new == 1);
    fill1 = ismember(lbl_dist, fill1);

    classifiedImg = uint8(fill1) + uint8(fill2) + uint8(fill3);
end
