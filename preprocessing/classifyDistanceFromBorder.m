function [distanceLabelmap, distanceMasks] = classifyDistanceFromBorder(convexhull, binary, params)
%CLASSIFYDISTANCEFROMBORDER Classify lung parenchyma by distance from pleural surface
%
%   [distanceLabelmap, distanceMasks] = classifyDistanceFromBorder(convexhull, binary, params)
%
%   Classifies lung parenchyma into concentric zones based on distance from
%   the pleural surface (outer borders). This is used for region-specific
%   processing and analysis.
%
%   Inputs:
%       convexhull  - convex hull connecting lungs
%       binary      - eroded lung mask (logical)
%       params      - Processing parameters with distance thresholds
%
%   Outputs:
%       distanceLabelmap - Classified volume (uint8: 0=bg, 1=close, 2=midfar, 3=far)
%       distanceMasks - Structure with fields:
%           .close    - Close to center (inner region)
%           .midfar   - Middle distance region
%           .far      - Far from center (peripheral/subpleural)
%           .normalized - Normalized distance map (0-1)
%
%
%   Classification for analysis:
%       close:  < 0.6
%       midfar: 0.6-0.85
%       far:    >= 0.85
%
%   See also: preprocessLungs, bwdist

fprintf('  Classifying distance from pleural surface...\n');

% Compute distance from the virtual lung ROI center
lungs_center_dist = bwdist(~convexhull);
lungs_center_dist = 1 - rescale(lungs_center_dist);
lungs_center_dist(binary == 0) = 0;

% Get distance range
limits0 = unique(lungs_center_dist);
limits0(1) = [];  % Remove background
limits0_min = min(limits0);
limits0_max = max(limits0);

% Rescale to range 0.1-1.0 (0 reserved for background)
lungs_center_dist = (0.9 / (limits0_max - limits0_min)) .* ...
    (lungs_center_dist - limits0_min) + 0.1;

distanceMasks.normalized = lungs_center_dist;

%% Classification for processing
% These are used during the segmentation pipeline
distanceMasks.close = lungs_center_dist < params.distance.closeThreshold;
distanceMasks.close(binary == 0) = 0;

distanceMasks.midfar = lungs_center_dist > params.distance.closeThreshold & ...
    lungs_center_dist < params.distance.midfarThreshold;

distanceMasks.far = lungs_center_dist > params.distance.farThreshold;

%% Classification for regional analysis
% Slightly different thresholds for final analysis
close_class = lungs_center_dist < params.distance.closeThresholdClass;
close_class(binary == 0) = 0;

midfar_class = lungs_center_dist >= params.distance.closeThresholdClass & ...
    lungs_center_dist < params.distance.midfarThreshold;

far_class = lungs_center_dist >= params.distance.farThreshold;

% Create labeled volume (1=close, 2=midfar, 3=far)
distanceLabelmap = uint8(close_class);
distanceLabelmap = distanceLabelmap + uint8(midfar_class) * 2;
distanceLabelmap = distanceLabelmap + uint8(far_class) * 3;

fprintf('    Distance classification complete\n');

end
