function labelMap = createLabelMap(segmentations, lungsbin)
%CREATELABELMAP Create final label map from all segmentations
%
%   labelMap = createLabelMap(segmentations, lungsbin)
%
%   Combines all segmentation results into a single label map with priority
%   ordering to handle overlaps.
%
%   Inputs:
%       segmentations - Structure with all segmentation masks
%       lungsbin      - Binary lung mask
%
%   Output:
%       labelMap - uint8 volume with labels:
%           0 - Background (outside lungs)
%           1 - Healthy parenchyma
%           2 - Ground glass opacity (GGO) /Compromised tissue
%           3 - Consolidation
%           4 - Air trapping/cysts
%           5 - Small/medium vessels
%           6 - Airways (bronchi)
%           7 - Trachea
%           8 - Large vessels
%           9 - Airway walls
%          10 - Trachea walls
%
%   Priority (later labels override earlier):
%       Background < Healthy < GGO < Consolidation < Air < Vessels < Airways
%
%   See also: saveResults

fprintf('Creating final label map...\n');

% Initialize with background
labelMap = zeros(size(lungsbin), 'uint8');

% Label 1: Healthy parenchyma (everything in lungs not otherwise labeled)
labelMap(lungsbin ~= 0) = 1;

% Label 2: GGO/Compromised tissue
if isfield(segmentations, 'ggo')
    labelMap(segmentations.ggo ~= 0) = 2;
end

% Label 3: Consolidation
if isfield(segmentations, 'consolidation')
    labelMap(segmentations.consolidation ~= 0) = 3;
end

% Label 4: Air trapping/cysts
if isfield(segmentations, 'air')
    labelMap(segmentations.air ~= 0) = 4;
end

% Label 5: Small/medium vessels
if isfield(segmentations, 'vessels')
    labelMap(segmentations.vessels ~= 0) = 5;
end

% Label 8: Large vessels
if isfield(segmentations, 'largeVessels')
    labelMap(segmentations.largeVessels ~= 0) = 8;
end

% Label 6: Airways
if isfield(segmentations, 'airways')
    labelMap(segmentations.airways ~= 0) = 6;
end

% Label 7: Trachea
if isfield(segmentations, 'trachea')
    labelMap(segmentations.trachea ~= 0) = 7;
end

% Label 9: Airway walls
if isfield(segmentations, 'airwayWalls')
    labelMap(segmentations.airwayWalls ~= 0) = 9;
end

% Label 10: Trachea walls
if isfield(segmentations, 'tracheaWalls')
    labelMap(segmentations.tracheaWalls ~= 0) = 10;
end

% Ensure background stays background
labelMap(lungsbin == 0) = 0;

fprintf('  Label map created\n\n');

end
