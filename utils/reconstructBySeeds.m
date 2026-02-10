function reconstructed = reconstructBySeeds(seeds, targetMask, connectivity)
%RECONSTRUCTBYSEEDS Morphological reconstruction from seed markers
%
%   reconstructed = reconstructBySeeds(seeds, targetMask, connectivity)
%
%   Performs morphological reconstruction to extract regions in targetMask that
%   are connected to seeds. This is a wrapper around imreconstruct for 3D data
%   with clearer semantics.
%
%   Inputs:
%       seeds        - Binary mask of seed points/regions
%       targetMask   - Binary mask to reconstruct from
%       connectivity - Connectivity (4, 6, 8, 18, or 26), default: 26
%
%   Output:
%       reconstructed - Binary mask of reconstructed regions
%
%   Example:
%       % Extract only vessels connected to main vasculature
%       mainVessels = bwareaopen(vessels, 1000);
%       allVessels = reconstructBySeeds(mainVessels, vessels, 26);
%
%       % Extract pathology connected to consolidations
%       pathology = reconstructBySeeds(consolidations, allAbnormalities, 8);
%
%   See also: imreconstruct, bwareaopen, bwlabeln

% Set default connectivity
if nargin < 3
    connectivity = 26;
end

% Validate inputs
if ~islogical(seeds) && ~isa(seeds, 'uint8')
    seeds = logical(seeds);
end

if ~islogical(targetMask) && ~isa(targetMask, 'uint8')
    targetMask = logical(targetMask);
end

% Ensure seeds are subset of target
seeds = seeds & targetMask;

% Perform reconstruction based on connectivity
switch connectivity
    case {4, 8}
        % 2D connectivity (slice-by-slice)
        reconstructed = imreconstruct(seeds, targetMask, connectivity);
    case {6, 18, 26}
        % 3D connectivity
        reconstructed = imreconstruct(uint8(seeds), uint8(targetMask), connectivity);
        reconstructed = logical(reconstructed);
    otherwise
        error('reconstructBySeeds:InvalidConnectivity', ...
              'Connectivity must be 4, 6, 8, 18, or 26');
end

end
