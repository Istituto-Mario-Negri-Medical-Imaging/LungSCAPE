function [lungsProcessed, cardinalPoints] = preprocessLungs(volumes, params)
%PREPROCESSLUNGS Preprocess lung masks and extract lung regions
%
%   [lungsProcessed, cardinalPoints] = preprocessLungs(volumes, params)
%
%   Performs initial preprocessing of lung masks including:
%   - Margin erosion to remove partial volume effects
%   - Separation into left and right lungs
%   - Calculation of cardinal points for regional analysis
%   - Creation of a lung ROI for processing
%
%   Inputs:
%       volumes - Structure with fields:
%           .lungs - Lung segmentation (1=right, 2=left)
%           .ct    - CT volume
%       params  - Processing parameters structure
%
%   Outputs:
%       lungsProcessed - Structure with fields:
%           .original         - Original lung segmentation 
%           .binary           - Eroded lung mask (binary)
%           .right            - Right lung mask
%           .left             - Left lung mask
%           .convexhull       - Convex Hull of the lungs
%       cardinalPoints - Structure with cardinal points for regional division
%
%   See also: refineLobeSeg

fprintf('Preprocessing lung masks...\n');

% If input `volumes.lungs` is binary (0/1), split into left/right labels
uvals = unique(volumes.lungs(:));
if numel(uvals) == 2 && all(ismember(uvals, [0, 1]))
    fprintf('  Detected binary lung mask — splitting into left/right labels...\n');

    % Connected components in 3D (26-connectivity)
    cc = bwconncomp(volumes.lungs > 0, 26);
    newLungs = zeros(size(volumes.lungs), 'like', volumes.lungs);
    mid_col = size(volumes.lungs, 2) / 2;

    for i = 1:cc.NumObjects
        vox = cc.PixelIdxList{i};
        [r, c, z] = ind2sub(size(volumes.lungs), vox);
        mean_c = mean(c);

        % Assign label: 1 = right (mean column > mid), 2 = left (<= mid)
        if mean_c > mid_col
            newLungs(vox) = 1;
        else
            newLungs(vox) = 2;
        end
    end

    volumes.lungs = newLungs;
else
    fprintf('  Lung mask already labeled; skipping split.\n');
end

% Store (possibly relabeled) original lung mask
lungsProcessed.original = volumes.lungs;

%% Erode lung margins
fprintf('  Eroding lung margins (size=%d)...\n', params.lungs.erosionSize);
lungsBinFill = imfill(imbinarize(volumes.lungs, 0), 26, 'holes');
se = strel('sphere', params.lungs.erosionSize);
lungsEroded = imerode(lungsBinFill, se);

lungsProcessed.binary = lungsEroded;
lungsProcessed.lungsBinFill = lungsBinFill;

%% Separate left and right lungs
fprintf('  Separating left and right lungs...\n');
lungsProcessed.right = (lungsProcessed.original == 1);
lungsProcessed.left = (lungsProcessed.original == 2);
lungsProcessed.right(lungsProcessed.binary == 0) = 0;
lungsProcessed.left(lungsProcessed.binary == 0) = 0;

%% Extract cardinal points for regional analysis
fprintf('  Computing cardinal points...\n');
cardinalPoints = extractCardinalPoints(lungsProcessed.right, lungsProcessed.left);

%% Perform convex hull connecting lungs
fprintf('  Convex Hull connecting lungs...\n');
lungsProcessed.convexhull = convexHulling(lungsProcessed.right, ...
    lungsProcessed.left, size(volumes.ct));

fprintf('  Lung preprocessing complete\n\n');

end

%% Helper Functions

function cardinalPoints = extractCardinalPoints(lungs_R, lungs_L)
%EXTRACTCARDINALPOINTS Extract cardinal points for regional division

% Right lung
lin_indx = find(lungs_R);
[val_row, val_col, val_z] = ind2sub(size(lungs_R), lin_indx);

cardinalPoints.right.SI_min = min(val_z);
cardinalPoints.right.SI_max = max(val_z);
cardinalPoints.right.PA_min = min(val_row);
cardinalPoints.right.PA_max = max(val_row);
cardinalPoints.right.LR_min = min(val_col);
cardinalPoints.right.LR_max = max(val_col);

% Left lung
lin_indx = find(lungs_L);
[val_row, val_col, val_z] = ind2sub(size(lungs_L), lin_indx);

cardinalPoints.left.SI_min = min(val_z);
cardinalPoints.left.SI_max = max(val_z);
cardinalPoints.left.PA_min = min(val_row);
cardinalPoints.left.PA_max = max(val_row);
cardinalPoints.left.LR_min = min(val_col);
cardinalPoints.left.LR_max = max(val_col);

% Compute threshold values for regional division
% Right lung
cardinalPoints.right.LR_middle = ceil(((cardinalPoints.right.LR_max - cardinalPoints.right.LR_min)/2) + cardinalPoints.right.LR_min);
cardinalPoints.right.PA_middle = ceil(((cardinalPoints.right.PA_max - cardinalPoints.right.PA_min)/2) + cardinalPoints.right.PA_min);
cardinalPoints.right.SI_middle1 = ceil(((cardinalPoints.right.SI_max - cardinalPoints.right.SI_min)/3) + cardinalPoints.right.SI_min);
cardinalPoints.right.SI_middle2 = ceil(((cardinalPoints.right.SI_max - cardinalPoints.right.SI_min)/3) + cardinalPoints.right.SI_middle1);
cardinalPoints.right.SI_middle = ceil(((cardinalPoints.right.SI_max - cardinalPoints.right.SI_min)/2) + cardinalPoints.right.SI_min);

% Left lung
cardinalPoints.left.LR_middle = ceil(((cardinalPoints.left.LR_max - cardinalPoints.left.LR_min)/2) + cardinalPoints.left.LR_min);
cardinalPoints.left.PA_middle = ceil(((cardinalPoints.left.PA_max - cardinalPoints.left.PA_min)/2) + cardinalPoints.left.PA_min);
cardinalPoints.left.SI_middle1 = ceil(((cardinalPoints.left.SI_max - cardinalPoints.left.SI_min)/3) + cardinalPoints.left.SI_min);
cardinalPoints.left.SI_middle2 = ceil(((cardinalPoints.left.SI_max - cardinalPoints.left.SI_min)/3) + cardinalPoints.left.SI_middle1);
cardinalPoints.left.SI_middle = ceil(((cardinalPoints.left.SI_max - cardinalPoints.left.SI_min)/2) + cardinalPoints.left.SI_min);

% Combined SI middle threshold
cardinalPoints.SI_middle = ceil((cardinalPoints.right.SI_middle + cardinalPoints.left.SI_middle) / 2);

end

function convexhull = convexHulling(lungs_R, lungs_L, volSize)
%CONVEXHULLING Create a virtual volume connecting lungs

fprintf('    Performing convex hull...\n');

convexhull = zeros(volSize);
se = strel('disk', 4);

for k = 1:volSize(3)
    lungs_R_k = lungs_R(:, :, k);
    lungs_L_k = lungs_L(:, :, k);
    lungsbinarized0 = lungs_R_k | lungs_L_k;

    % Remove small regions
    lungs_R_k = bwareaopen(lungs_R_k, 10000, 8);
    lungs_L_k = bwareaopen(lungs_L_k, 10000, 8);

    % Check if both lungs are present
    go_on1 = size(unique(lungs_L_k));
    go_on2 = size(unique(lungs_R_k));

    if go_on1(1) == 2 && go_on2(1) == 2
        % Keep only largest component for each lung
        [~, n_objs_L] = bwlabel(lungs_L_k);
        [~, n_objs_R] = bwlabel(lungs_R_k);

        if n_objs_L > 1
            [objs_L, ~] = bwlabel(lungs_L_k);
            area_L = regionprops(objs_L, 'Area');
            [~, area_idx] = max([area_L.Area]);
            lungs_L_k = (objs_L == area_idx);
        end

        if n_objs_R > 1
            [objs_R, ~] = bwlabel(lungs_R_k);
            area_R = regionprops(objs_R, 'Area');
            [~, area_idx] = max([area_R.Area]);
            lungs_R_k = (objs_R == area_idx);
        end

        lungsbinarized = lungs_R_k | lungs_L_k;

        % Compute convex hull of the combined lungs to create a virtual ROI
        % Use 'union' to obtain a single hull covering both lungs efficiently
        filled = bwconvhull(lungsbinarized, 'union');

        % Smooth/close small gaps and artifacts between lungs
        filled = imclose(filled, se);

        % Fill holes inside the convex hull
        filled = imfill(filled, 'holes');
    else
        filled = lungsbinarized0;
        filled = imfill(filled, 'holes');
    end

    convexhull(:, :, k) = filled;
end

end
