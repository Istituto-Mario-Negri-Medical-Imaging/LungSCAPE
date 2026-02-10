function walls = segmentWalls(structures, structureType, params)
%SEGMENTWALLS Segment walls of tubular structures (airways/vessels)
%
%   walls = segmentWalls(structures, structureType, params)
%
%   Segments the walls of tubular structures (airways or vessels) using
%   distance-based partial volume estimation. Walls appear as a thin shell
%   around the lumen.
%
%   Method:
%       1. Compute distance transform from structure
%       2. Extract voxels within wall thickness range
%       3. Morphological refinement
%
%   This captures partial volume effects at the boundary between the
%   structure lumen and surrounding tissue.
%
%   Inputs:
%       structures   - Binary segmentation of structures (airways or vessels)
%       structureType - 'airway' or 'vessel' (determines wall thickness)
%       params       - Processing parameters with wall thickness thresholds
%
%   Outputs:
%       walls - Binary segmentation of structure walls (logical)
%
%   Example:
%       airwayWalls = segmentWalls(airways, 'airway', params);
%       vesselWalls = segmentWalls(vessels, 'vessel', params);
%
%   See also: bwdist, segmentAirways, segmentVessels

fprintf('    Segmenting %s walls...\n', structureType);

%% Determine wall thickness based on structure type
switch lower(structureType)
    case 'airway'
        minDist = params.airways.wallMinDistance;
        maxDist = params.airways.wallMaxDistance;
    case 'vessel'
        minDist = params.vessels.wallMinDistance;
        maxDist = params.vessels.wallMaxDistance;
    case 'large_vessel'
        % Large vessels may have thicker walls
        minDist = 0;
        maxDist = 1.5;
    otherwise
        error('Unknown structure type: %s. Use ''airway'' or ''vessel''.', structureType);
end

%% Compute distance transform
% Distance from the structure boundary
[distanceTransform, ~] = bwdist(single(structures));

%% Extract wall region
% Walls are voxels within the specified distance range from the structure
walls = distanceTransform < maxDist & distanceTransform > minDist;

%% Morphological refinement
% Remove small isolated fragments
if isfield(params, 'wallMinSize')
    walls = bwareaopen(walls, params.wallMinSize, 6);
end

% Optional: ensure walls are connected to the structure
if isfield(params, 'ensureConnected') && params.ensureConnected
    structuresDilated = imdilate(structures, strel('sphere', 1));
    walls = imreconstruct(structuresDilated, walls, 26);
end

fprintf('      %s wall segmentation complete\n', structureType);

end


function [airwayWalls, vesselWalls, largeVesselWalls, tracheaWalls] = segmentAllWalls(...
    airways, vessels, largeVessels, trachea, params)
%SEGMENTALLWALLS Segment walls for all tubular structures
%
%   [airwayWalls, vesselWalls, largeVesselWalls, tracheaWalls] = ...
%       segmentAllWalls(airways, vessels, largeVessels, trachea, params)
%
%   Convenience function to segment walls for all tubular structures in
%   one call. Handles overlap resolution between different wall types.
%
%   Inputs:
%       airways      - Airway segmentation
%       vessels      - Vessel segmentation
%       largeVessels - Large vessel segmentation
%       trachea      - Trachea segmentation
%       params       - Processing parameters
%
%   Outputs:
%       airwayWalls      - Airway wall segmentation
%       vesselWalls      - Vessel wall segmentation (small/medium)
%       largeVesselWalls - Large vessel wall segmentation
%       tracheaWalls     - Trachea wall segmentation
%
%   See also: segmentWalls

fprintf('  Segmenting structure walls...\n');

%% Segment walls for each structure type
airwayWalls = segmentWalls(airways, 'airway', params);
tracheaWalls = segmentWalls(trachea, 'airway', params);
vesselWalls = segmentWalls(vessels, 'vessel', params);
largeVesselWalls = segmentWalls(largeVessels, 'large_vessel', params);

%% Resolve overlaps
fprintf('    Resolving wall overlaps...\n');

% Priority: airways > trachea > large vessels > small vessels
% Airways walls take precedence
tracheaWalls(airwayWalls ~= 0) = 0;
largeVesselWalls(airwayWalls ~= 0) = 0;
largeVesselWalls(tracheaWalls ~= 0) = 0;
vesselWalls(airwayWalls ~= 0) = 0;
vesselWalls(tracheaWalls ~= 0) = 0;
vesselWalls(largeVesselWalls ~= 0) = 0;

fprintf('    Wall segmentation complete\n');

end
