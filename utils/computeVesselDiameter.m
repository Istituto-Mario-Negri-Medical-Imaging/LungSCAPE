function [diameterImage, skeletonImage] = computeVesselDiameter(binaryMask, minBranchLength)
%COMPUTEVESSELDIAMETER Compute diameter map of tubular structures
%
%   [diameterImage, skeletonImage] = computeVesselDiameter(binaryMask, minBranchLength)
%
%   Computes the local diameter of tubular structures (vessels, airways) by:
%   1. Skeletonizing the binary mask
%   2. Computing distance transform
%   3. Multiplying skeleton by distance (diameter = 2 * distance)
%
%   Inputs:
%       binaryMask     - 3D binary mask of tubular structures
%       minBranchLength - Minimum branch length to retain in skeleton (default: 0)
%
%   Outputs:
%       diameterImage  - Double array with diameter values at skeleton points
%       skeletonImage  - Binary skeleton of the input mask
%
%   Example:
%       % Classify vessels by diameter
%       [diameter, skel] = computeVesselDiameter(vessels, 5);
%       largeVessels = diameter > 6;
%       smallVessels = diameter > 0 & diameter <= 2;
%
%   See also: bwskel, bwdist, bwmorph3

% Set default minimum branch length
if nargin < 2
    minBranchLength = 0;
end

% Validate input
if ~islogical(binaryMask) && ~isa(binaryMask, 'uint8')
    binaryMask = logical(binaryMask);
end

% Compute distance transform (distance to nearest background pixel)
distanceTransform = bwdist(~binaryMask);

% Skeletonize with optional minimum branch length
if minBranchLength > 0
    skeletonImage = bwskel(binaryMask, 'MinBranchLength', minBranchLength);
else
    skeletonImage = bwskel(binaryMask);
end

% Compute diameter (diameter = 2 * radius = 2 * distance)
diameterImage = 2 * distanceTransform .* double(skeletonImage);

end
