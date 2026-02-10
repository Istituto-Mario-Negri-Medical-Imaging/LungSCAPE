function [kept, discarded] = filterByAnisotropy(binaryMask, FA_threshold, varargin)
%FILTERBYANISOTROPY Filter connected components based on fractional anisotropy
%
%   [kept, discarded] = filterByAnisotropy(binaryMask, FA_threshold)
%   [kept, discarded] = filterByAnisotropy(binaryMask, FA_threshold, 'invert', true)
%
%   Filters 3D binary connected components based on their fractional anisotropy (FA).
%   Components with FA above the threshold are kept (or discarded if inverted).
%
%   Inputs:
%       binaryMask   - 3D logical or uint8 binary volume
%       FA_threshold - Scalar threshold value (0-1)
%
%   Optional Name-Value Pairs:
%       'invert'     - If true, keep components with FA < threshold (default: false)
%       'minVolume'  - Minimum volume (voxels) for components to consider (default: 0)
%
%   Outputs:
%       kept      - Binary mask of kept components
%       discarded - Binary mask of discarded components
%
%   Example:
%       % Keep only highly anisotropic (tubular) structures
%       vessels = segmentVessels(volume);
%       [tubelike, bloblike] = filterByAnisotropy(vessels, 0.90);
%
%       % Keep only isotropic (blob-like) structures
%       lesions = segmentLesions(volume);
%       [blobs, tubes] = filterByAnisotropy(lesions, 0.70, 'invert', true);
%
%   See also: computeFractionalAnisotropy, regionprops3, bwlabeln

% Parse optional inputs
p = inputParser;
addParameter(p, 'invert', false, @islogical);
addParameter(p, 'minVolume', 0, @isnumeric);
parse(p, varargin{:});

invert = p.Results.invert;
minVolume = p.Results.minVolume;

% Validate inputs
if ~islogical(binaryMask) && ~isa(binaryMask, 'uint8')
    error('filterByAnisotropy:InvalidInput', 'binaryMask must be logical or uint8');
end

if FA_threshold < 0 || FA_threshold > 1
    error('filterByAnisotropy:InvalidThreshold', 'FA_threshold must be between 0 and 1');
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

% Compute eigenvalues for each component
props = regionprops3(labeled, 'EigenValues', 'Volume');

% Filter components based on FA
keep_idx = false(numComponents, 1);
for i = 1:numComponents
    % Skip if below minimum volume
    if props.Volume(i) < minVolume
        continue;
    end

    eigenvalues = props.EigenValues{i};
    FA = computeFractionalAnisotropy(eigenvalues);

    if invert
        keep_idx(i) = FA < FA_threshold;
    else
        keep_idx(i) = FA > FA_threshold;
    end
end

% Create output masks
kept = ismember(labeled, find(keep_idx));
discarded = binaryMask & ~kept;

end
