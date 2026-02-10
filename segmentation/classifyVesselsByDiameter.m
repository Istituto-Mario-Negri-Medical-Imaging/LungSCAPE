function vesselsClassified = classifyVesselsByDiameter(vessels, params)
%CLASSIFYVESSELSBYDIAMETER Classify vessels into small/medium/large categories
%
%   vesselsClassified = classifyVesselsByDiameter(vessels, params)
%
%   Classifies vessels into three diameter categories using sophisticated
%   skeleton-based diameter estimation and nearest-neighbor propagation for
%   non-skeleton voxels. Uses parallel processing for efficiency.
%
%   Classification:
%       Label 1: Large vessels (diameter > 6mm)
%       Label 2: Medium vessels (2mm < diameter <= 6mm)
%       Label 3: Small vessels (diameter <= 2mm)
%
%   Algorithm:
%       1. Compute diameter map from distance transform
%       2. Extract skeleton and measure diameter at skeleton points
%       3. Classify skeleton points by diameter
%       4. Handle boundary cases (reclassify small fragments)
%       5. Propagate labels to non-skeleton voxels using nearest neighbor
%
%   Inputs:
%       vessels - Binary vessel segmentation (logical)
%       params  - Processing parameters with fields:
%           .smallDiameterThreshold  - Threshold for small vessels (default: 2mm)
%           .largeDiameterThreshold  - Threshold for large vessels (default: 6mm)
%           .numWorkers             - Number of parallel workers (default: 3)
%
%   Outputs:
%       vesselsClassified - Labeled volume (uint8):
%           0 = background
%           1 = large vessels (> 6mm)
%           2 = medium vessels (2-6mm)
%           3 = small vessels (<= 2mm)
%
%   Example:
%       params.smallDiameterThreshold = 2;
%       params.largeDiameterThreshold = 6;
%       classified = classifyVesselsByDiameter(vessels, params);
%       largeVessels = classified == 1;
%       smallVessels = classified == 3;
%
%   See also: computeVesselDiameter, propagateLabelsNN

fprintf('  Classifying vessels by diameter...\n');

% Set default parameters
if ~isfield(params, 'smallDiameterThreshold')
    params.smallDiameterThreshold = 2;
end
if ~isfield(params, 'largeDiameterThreshold')
    params.largeDiameterThreshold = 6;
end
if ~isfield(params, 'numWorkers')
    params.numWorkers = 3;
end

%% Compute diameter map
fprintf('    Computing vessel diameter map...\n');

diameter = 2 * bwdist(~vessels);
skeletonImage = bwskel(vessels);
diameterImage = diameter .* double(skeletonImage);

%% Initial classification by diameter
fprintf('    Classifying skeleton points by diameter...\n');

diameterImageSmall = diameterImage <= params.smallDiameterThreshold & diameterImage > 0;
diameterImageMid = diameterImage <= params.largeDiameterThreshold & ...
    diameterImage > params.smallDiameterThreshold;
diameterImageLarge = diameterImage > params.largeDiameterThreshold;

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

diameterImageLarge = diameterImageLarge & ~large2change;

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

    % Use parallel processing for speed
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        if params.numWorkers > 0
            parpool(params.numWorkers);
        end
    end

    % For each fragment, find nearest endpoint and get its label
    parfor (lbl_idx_k = 1:length(lbl_idx), params.numWorkers)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(fragments), val2change);
        lbl_idx_new(lbl_idx_k) = labelImg(val_row, val_col, val_z);
    end

    % Clean up parallel pool
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
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

    % Use parallel processing
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        if params.numWorkers > 0
            parpool(params.numWorkers);
        end
    end

    % For each voxel, find nearest labeled voxel
    parfor (lbl_idx_k = 1:length(lbl_idx), params.numWorkers)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(toClassify), val2change);
        lbl_idx_new(lbl_idx_k) = labelImg(val_row, val_col, val_z);
    end

    % Clean up parallel pool
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
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
