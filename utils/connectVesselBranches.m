function connectionTracts = connectVesselBranches(skeletonMain, skeletonAll, params)
%CONNECTVESSELBRANCHES Connect fragmented vessel branches using graph analysis
%
%   connectionTracts = connectVesselBranches(skeletonMain, skeletonAll, params)
%
%   Uses graph-based shortest path analysis to reconnect fragmented main
%   vascular branches. This sophisticated algorithm:
%   1. Creates graphs from skeleton structures
%   2. Identifies endpoints of main branches
%   3. Finds nearest neighbors between disconnected branches
%   4. Computes shortest paths to reconnect them
%   5. Filters connections by anisotropy (FA > threshold)
%
%   This is one of the most advanced features of the pipeline, using graph
%   theory to intelligently reconnect vessels that were fragmented due to
%   pathology, partial volume effects, or segmentation challenges.
%
%   Inputs:
%       skeletonMain - Binary skeleton of main vessels (diameter > 2.5mm)
%       skeletonAll  - Binary skeleton of all vessels
%       params       - Processing parameters with fields:
%           .maxConnectionDistance  - Max distance to connect (default: 25)
%           .maxPathLength          - Max shortest path length (default: 30)
%           .connectionFAThreshold  - FA threshold for connections (default: 0.87)
%
%   Outputs:
%       connectionTracts - Binary volume with reconnection paths (logical)
%
%   Algorithm:
%       1. Label connected components of main branches
%       2. Create graph G1 from main skeleton
%       3. Create graph Gout from delta skeleton (main endpoints + small vessels)
%       4. For each disconnected branch:
%          - Find nearest neighbor in other branches
%          - Compute shortest path in Gout
%          - Keep if path is short and anisotropic (vessel-like)
%       5. Filter paths by fractional anisotropy (FA > 0.87)
%
%   See also: binaryImageGraph3, shortestpath, filterByAnisotropy

fprintf('      Creating vessel connectivity graphs...\n');

% Set default parameters if not provided
if ~isfield(params, 'maxConnectionDistance')
    params.maxConnectionDistance = 25;  % mm
end
if ~isfield(params, 'maxPathLength')
    params.maxPathLength = 30;  % nodes
end
if ~isfield(params, 'connectionFAThreshold')
    params.connectionFAThreshold = 0.87;
end

%% Label main branches
[mainBranches, ~] = bwlabeln(skeletonMain);

%% Create delta skeleton (small vessels + main endpoints)
% Delta skeleton includes connections between main branches
skeletonDelta = skeletonAll & ~skeletonMain;
endpointsImageMain = bwmorph3(skeletonMain, 'endpoints');
skeletonDelta = skeletonDelta | endpointsImageMain;

%% Create graphs
% G1: Graph of main vessels only
G1 = binaryImageGraph3(skeletonMain);

% Gout: Graph of delta skeleton (used for shortest path finding)
Gout = binaryImageGraph3(skeletonDelta);

%% Find endpoints of main branches
nodeDegree = degree(G1);
endNodes = find(nodeDegree == 1);

% Get coordinates of endpoints
G1_coor = [G1.Nodes.x(endNodes), G1.Nodes.y(endNodes), G1.Nodes.z(endNodes)];

% Get branch labels for each endpoint
G1_coor_lbl = zeros(size(G1_coor, 1), 1);
for i = 1:size(G1_coor, 1)
    G1_coor_lbl(i) = mainBranches(sub2ind(size(mainBranches), ...
        G1_coor(i, 2), G1_coor(i, 1), G1_coor(i, 3)));
end

G1_coor_lbl = [G1_coor, G1_coor_lbl];

%% Find nearest neighbors between disconnected branches
fprintf('      Finding nearest neighbors between branches...\n');

distances = [];
pair1 = [];
pair2 = [];

uniqueLabels = unique(G1_coor_lbl(:, 4));

for i = 1:length(uniqueLabels)
    % Get coordinates for current branch
    currentBranchCoords = G1_coor_lbl(G1_coor_lbl(:, 4) == uniqueLabels(i), 1:3);

    % Get coordinates for all other branches
    otherBranchesCoords = G1_coor_lbl(G1_coor_lbl(:, 4) ~= uniqueLabels(i), 1:3);

    % Find nearest neighbor for each endpoint in current branch
    for ii = 1:size(currentBranchCoords, 1)
        [nearest, distance] = dsearchn(otherBranchesCoords, currentBranchCoords(ii, :));
        nearestNeighbor = otherBranchesCoords(nearest, :);

        pair1 = cat(1, pair1, currentBranchCoords(ii, :));
        pair2 = cat(1, pair2, nearestNeighbor);
        distances = cat(2, distances, distance);
    end
end

%% Filter pairs by distance threshold
endnodesConnectIdx = distances <= params.maxConnectionDistance;
inConnect = pair1(endnodesConnectIdx, :);
outConnect = pair2(endnodesConnectIdx, :);

% Remove duplicate pairs
sortPairs = [inConnect, outConnect];
pairs = [inConnect, outConnect];
[~, pairsIdx] = unique(sort(sortPairs, 2), 'rows');
pairs = pairs(pairsIdx, :);
inConnect = pairs(:, 1:3);
outConnect = pairs(:, 4:end);

nodesNum = size(inConnect, 1);

fprintf('      Found %d potential connections\n', nodesNum);

%% Compute shortest paths
fprintf('      Computing shortest paths...\n');

nodesIn = [];
nodesOut = [];

% Find node indices in Gout for each coordinate pair
for nodeIdx = 1:nodesNum
    nodeInCoor = inConnect(nodeIdx, :);
    nodeOutCoor = outConnect(nodeIdx, :);

    % Find nodes in Gout that match coordinates
    nodeIndexIn = intersect(find(Gout.Nodes.x == nodeInCoor(1)), ...
        intersect(find(Gout.Nodes.y == nodeInCoor(2)), ...
        find(Gout.Nodes.z == nodeInCoor(3))));

    nodeIndexOut = intersect(find(Gout.Nodes.x == nodeOutCoor(1)), ...
        intersect(find(Gout.Nodes.y == nodeOutCoor(2)), ...
        find(Gout.Nodes.z == nodeOutCoor(3))));

    if length(nodeIndexIn) == 1 && length(nodeIndexOut) == 1
        nodesIn = cat(1, nodesIn, nodeIndexIn);
        nodesOut = cat(1, nodesOut, nodeIndexOut);
    end
end

% Compute shortest paths and create connection tracts
connectionTracts = zeros(size(skeletonMain));

validPaths = 0;
for i = 1:length(nodesOut)
    % Find shortest path in Gout
    spath = shortestpath(Gout, nodesIn(i), nodesOut(i));
    pathLength = length(spath);

    % Keep path if it's short and not too much longer than direct distance
    if pathLength <= params.maxPathLength && pathLength < (distances(i) + 10)
        connectionTracts(sub2ind(size(skeletonMain), ...
            Gout.Nodes.y(spath), Gout.Nodes.x(spath), Gout.Nodes.z(spath))) = 1;
        validPaths = validPaths + 1;
    end
end

fprintf('      Created %d valid connection paths\n', validPaths);

%% Filter connections by anisotropy (FA)
fprintf('      Filtering connections by anisotropy...\n');

% Label each connection tract
connectionTractsLbl = bwlabeln(connectionTracts);
connectionTractsEig = regionprops3(connectionTractsLbl, 'EigenValues');

idxKeep = zeros(1, size(connectionTractsEig, 1));

for i = 1:size(connectionTractsEig, 1)
    eigenvalues = connectionTractsEig.EigenValues{i};

    % Compute fractional anisotropy
    FA = computeFractionalAnisotropy(eigenvalues);

    % Keep only highly anisotropic structures (vessel-like)
    if FA > params.connectionFAThreshold
        idxKeep(i) = 1;
    end
end

connectionTractsKeep = ismember(connectionTractsLbl, find(idxKeep));

fprintf('      Kept %d/%d connections (FA > %.2f)\n', ...
    sum(idxKeep), length(idxKeep), params.connectionFAThreshold);

connectionTracts = logical(connectionTractsKeep);

end


%% Helper function: Create binary image graph
function G = binaryImageGraph3(binaryImage)
%BINARYIMAGEGRAPH3 Create graph from 3D binary skeleton
%
%   This function creates a graph where each skeleton voxel is a node
%   and edges connect 26-connected neighbors.
%
%   Note: This requires the Image Processing Toolbox. If not available,
%   you'll need to implement a custom graph construction.

    % Get coordinates of all skeleton points
    [y, x, z] = ind2sub(size(binaryImage), find(binaryImage));

    % Create node table
    nodeTable = table(x, y, z, 'VariableNames', {'x', 'y', 'z'});

    % Build edge list by checking 26-connectivity
    numNodes = length(x);
    edges = [];

    % This is simplified - in practice, you'd use optimized connectivity checking
    % For each voxel, check 26 neighbors
    for i = 1:numNodes
        % Check 26-neighborhood
        for dx = -1:1
            for dy = -1:1
                for dz = -1:1
                    if dx == 0 && dy == 0 && dz == 0
                        continue;
                    end

                    nx = x(i) + dx;
                    ny = y(i) + dy;
                    nz = z(i) + dz;

                    % Check if neighbor exists in skeleton
                    if nx >= 1 && nx <= size(binaryImage, 2) && ...
                       ny >= 1 && ny <= size(binaryImage, 1) && ...
                       nz >= 1 && nz <= size(binaryImage, 3)

                        if binaryImage(ny, nx, nz)
                            % Find neighbor node index
                            neighborIdx = find(x == nx & y == ny & z == nz, 1);
                            if ~isempty(neighborIdx) && neighborIdx > i
                                edges = [edges; i, neighborIdx];
                            end
                        end
                    end
                end
            end
        end
    end

    % Create graph
    if isempty(edges)
        G = graph([], nodeTable);
    else
        G = graph(edges(:, 1), edges(:, 2), [], nodeTable);
    end
end
