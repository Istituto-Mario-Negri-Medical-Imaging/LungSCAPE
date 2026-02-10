function [kept, discarded] = filterByEigenVectorOrientation(binaryMask, FA_threshold, ...
    orientation, minEigenvalue)
%FILTERBYEIGENVECTORORIENTATION Filter components by primary eigenvector orientation
%
%   [kept, discarded] = filterByEigenVectorOrientation(binaryMask, FA_threshold, ...
%       orientation, minEigenvalue)
%
%   Filters 3D binary components based on both fractional anisotropy and the
%   orientation of the primary eigenvector. Useful for detecting vertically or
%   horizontally oriented structures (e.g., reticulations vs vessels).
%
%   Inputs:
%       binaryMask    - 3D logical or uint8 binary volume
%       FA_threshold  - Minimum FA value for anisotropic structures (0-1)
%       orientation   - Target orientation: 'vertical', 'horizontal', or
%                       'longitudinal' (along z-axis)
%       minEigenvalue - Minimum value of largest eigenvalue (structure size)
%
%   Outputs:
%       kept      - Binary mask of components matching criteria
%       discarded - Binary mask of components not matching criteria
%
%   Example:
%       % Find vertically-oriented reticulations (along z-axis)
%       [reticulations, ~] = filterByEigenVectorOrientation(candidates, ...
%           0.92, 'longitudinal', 100);
%
%       % Find horizontally-oriented structures (in xy-plane)
%       [horizontal, ~] = filterByEigenVectorOrientation(candidates, ...
%           0.87, 'horizontal', 50);
%
%   See also: computeFractionalAnisotropy, filterByAnisotropy, regionprops3

% Validate inputs
if ~islogical(binaryMask) && ~isa(binaryMask, 'uint8')
    error('filterByEigenVectorOrientation:InvalidInput', ...
          'binaryMask must be logical or uint8');
end

validOrientations = {'vertical', 'horizontal', 'longitudinal'};
if ~ismember(lower(orientation), validOrientations)
    error('filterByEigenVectorOrientation:InvalidOrientation', ...
          'orientation must be ''vertical'', ''horizontal'', or ''longitudinal''');
end

% Convert to logical if needed
binaryMask = logical(binaryMask);

% Label connected components
labeled = bwlabeln(binaryMask);
numComponents = max(labeled(:));

% Handle empty mask
if numComponents == 0
    kept = false(size(binaryMask));
    discarded = false(size(binaryMask));
    return;
end

% Compute eigenvalues and eigenvectors for each component
props = regionprops3(labeled, 'EigenValues', 'EigenVectors');

% Filter components
keep_idx = false(numComponents, 1);
for i = 1:numComponents
    eigenvalues = props.EigenValues{i};
    eigenvectors = props.EigenVectors{i};

    % Compute FA
    FA = computeFractionalAnisotropy(eigenvalues);

    % Check FA and eigenvalue thresholds
    if FA < FA_threshold || eigenvalues(1) < minEigenvalue
        continue;
    end

    % Check orientation of primary eigenvector (corresponding to largest eigenvalue)
    primaryEigenvector = eigenvectors(:, 1);

    % Get absolute components for orientation comparison
    absVec = abs(primaryEigenvector);

    switch lower(orientation)
        case 'longitudinal'
            % Primary eigenvector should be along z-axis (3rd component largest)
            keep_idx(i) = absVec(3) > absVec(2) && absVec(3) > absVec(1);

        case 'horizontal'
            % Primary eigenvector should be in xy-plane (z-component smallest)
            keep_idx(i) = absVec(3) < absVec(2) && absVec(3) < absVec(1);

        case 'vertical'
            % For vertical, check if z-component is dominant
            keep_idx(i) = absVec(3) > max(absVec(1), absVec(2));
    end
end

% Create output masks
kept = ismember(labeled, find(keep_idx));
discarded = binaryMask & ~kept;

end
