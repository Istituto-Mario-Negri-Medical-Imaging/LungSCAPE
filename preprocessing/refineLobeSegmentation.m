function lobesRefined = refineLobeSegmentation(lobes, lungsbinary, params)
%REFINELOBESEGMENTATION Refine lobe segmentation to fill gaps
%
%   lobesRefined = refineLobeSegmentation(lobes, lungsbinary, params)
%
%   Fills gaps in lobe segmentation by assigning unlabeled voxels within
%   the lung mask to the nearest lobe based on distance transform.
%
%   Inputs:
%       lobes        - Initial lobe segmentation (uint8: 0-5)
%       lungsbinary - Binary filled lung mask
%       params       - Processing parameters
%
%   Output:
%       lobesRefined - Refined lobe segmentation with gaps filled
%
%   Example:
%       lobesRefined = refineLobeSegmentation(lobes, lungsMask, params);
%
%   See also: preprocessLungs, bwdist

if isempty(lobes) || all(lobes(:) == 0)
    lobesRefined = [];
    return;
end

fprintf('  Refining lobe segmentation...\n');

% Mask lobes to lung region
lobes = lobes .* uint8(lungsbinary);
lobesbin = imbinarize(lobes, 0);

% Find voxels that need lobe assignment
delta = lungsbinary;
delta(lobesbin ~= 0) = 0;

classify = delta;
labelimg = lobes;

% Use distance transform to find nearest lobe for each unlabeled voxel
[~, lbl_dist] = bwdist(lobesbin);
lbl_dist = lbl_dist .* uint32(classify);
lbl_idx = unique(lbl_dist);
lbl_idx(1) = [];  % Remove background
lbl_idx_new = zeros(size(lbl_idx));

% Parallel processing for speed
fprintf('    Assigning unlabeled voxels to nearest lobe...\n');
if params.parallel.useParallel && params.parallel.workers > 0
    parfor (lbl_idx_k = 1:length(lbl_idx), params.parallel.workers)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(classify), val2change);
        lbl_idx_new(lbl_idx_k) = labelimg(val_row, val_col, val_z);
    end
else
    % Serial processing
    for lbl_idx_k = 1:length(lbl_idx)
        val2change = lbl_idx(lbl_idx_k);
        [val_row, val_col, val_z] = ind2sub(size(classify), val2change);
        lbl_idx_new(lbl_idx_k) = labelimg(val_row, val_col, val_z);
    end
end

% Assign voxels to lobes
fill5 = lbl_idx(lbl_idx_new == 5);
fill5 = ismember(lbl_dist, fill5) * 5;

fill4 = lbl_idx(lbl_idx_new == 4);
fill4 = ismember(lbl_dist, fill4) * 4;

fill3 = lbl_idx(lbl_idx_new == 3);
fill3 = ismember(lbl_dist, fill3) * 3;

fill2 = lbl_idx(lbl_idx_new == 2);
fill2 = ismember(lbl_dist, fill2) * 2;

fill1 = lbl_idx(lbl_idx_new == 1);
fill1 = ismember(lbl_dist, fill1);

% Combine all lobes
lobesRefined = labelimg + uint8(fill1) + uint8(fill2) + ...
    uint8(fill3) + uint8(fill4) + uint8(fill5);

fprintf('    Lobe refinement complete\n');

end
