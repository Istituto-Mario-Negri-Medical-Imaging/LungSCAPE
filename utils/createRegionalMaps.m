function regionalMaps = createRegionalMaps(lungs_R, lungs_L, lungsbin, cardinalPoints)
%CREATEREGIONALMAPS Create regional division maps for lung analysis
%
%   regionalMaps = createRegionalMaps(lungs_R, lungs_L, lungsbin, cardinalPoints)
%
%   Creates anatomical regional divisions of the lungs in three orientations:
%   - SI: Superior-Inferior (cranio-caudal)
%   - PA: Posterior-Anterior
%   - LR: Left-Right (medial-lateral)
%
%   Inputs:
%       lungs_R        - Right lung binary mask
%       lungs_L        - Left lung binary mask
%       lungsbin       - Combined binary lung mask
%       cardinalPoints - Structure with lung boundaries and midpoints
%
%   Output:
%       regionalMaps - Structure with fields:
%           .SI - Superior-Inferior map (uint8: 1-6)
%           .PA - Posterior-Anterior map (uint8: 1-4)
%           .LR - Left-Right map (uint8: 1-4)
%
%   SI Labels: 1-3 = Left (inf-mid-sup), 4-6 = Right (inf-mid-sup)
%   PA Labels: 1-2 = Left (ant-post), 3-4 = Right (ant-post)
%   LR Labels: 1-2 = Left (lat-med), 3-4 = Right (lat-med)
%
%   See also: preprocessLungs, extractCardinalPoints

fprintf('Creating regional maps...\n');

lungs_airwaybin = lungsbin;

%% PA (Posterior-Anterior) Map
lSegmentsMapPA = uint8(zeros(size(lungs_airwaybin)));
rSegmentsMapPA = uint8(zeros(size(lungs_airwaybin)));

lSegmentsMapPA(cardinalPoints.left.PA_min:cardinalPoints.left.PA_middle, :, :) = ...
    lungs_airwaybin(cardinalPoints.left.PA_min:cardinalPoints.left.PA_middle, :, :) * 1;
lSegmentsMapPA(cardinalPoints.left.PA_middle+1:cardinalPoints.left.PA_max, :, :) = ...
    lungs_airwaybin(cardinalPoints.left.PA_middle+1:cardinalPoints.left.PA_max, :, :) * 2;

rSegmentsMapPA(cardinalPoints.right.PA_min:cardinalPoints.right.PA_middle, :, :) = ...
    lungs_airwaybin(cardinalPoints.right.PA_min:cardinalPoints.right.PA_middle, :, :) * 3;
rSegmentsMapPA(cardinalPoints.right.PA_middle+1:cardinalPoints.right.PA_max, :, :) = ...
    lungs_airwaybin(cardinalPoints.right.PA_middle+1:cardinalPoints.right.PA_max, :, :) * 4;

regionalMaps.PA = (lSegmentsMapPA .* uint8(lungs_L)) + (rSegmentsMapPA .* uint8(lungs_R));

%% LR (Left-Right / Medial-Lateral) Map
lSegmentsMapLR = uint8(zeros(size(lungs_airwaybin)));
rSegmentsMapLR = uint8(zeros(size(lungs_airwaybin)));

lSegmentsMapLR(:, cardinalPoints.left.LR_min:cardinalPoints.left.LR_middle, :) = ...
    lungs_airwaybin(:, cardinalPoints.left.LR_min:cardinalPoints.left.LR_middle, :) * 1;
lSegmentsMapLR(:, cardinalPoints.left.LR_middle+1:cardinalPoints.left.LR_max, :) = ...
    lungs_airwaybin(:, cardinalPoints.left.LR_middle+1:cardinalPoints.left.LR_max, :) * 2;

rSegmentsMapLR(:, cardinalPoints.right.LR_min:cardinalPoints.right.LR_middle, :) = ...
    lungs_airwaybin(:, cardinalPoints.right.LR_min:cardinalPoints.right.LR_middle, :) * 3;
rSegmentsMapLR(:, cardinalPoints.right.LR_middle+1:cardinalPoints.right.LR_max, :) = ...
    lungs_airwaybin(:, cardinalPoints.right.LR_middle+1:cardinalPoints.right.LR_max, :) * 4;

regionalMaps.LR = (lSegmentsMapLR .* uint8(lungs_L)) + (rSegmentsMapLR .* uint8(lungs_R));

%% SI (Superior-Inferior / Cranio-Caudal) Map
lSegmentsMapSI = uint8(zeros(size(lungs_airwaybin)));
rSegmentsMapSI = uint8(zeros(size(lungs_airwaybin)));

lSegmentsMapSI(:, :, cardinalPoints.left.SI_min:cardinalPoints.left.SI_middle1) = ...
    lungs_airwaybin(:, :, cardinalPoints.left.SI_min:cardinalPoints.left.SI_middle1) * 1;
lSegmentsMapSI(:, :, cardinalPoints.left.SI_middle1:cardinalPoints.left.SI_middle2) = ...
    lungs_airwaybin(:, :, cardinalPoints.left.SI_middle1:cardinalPoints.left.SI_middle2) * 2;
lSegmentsMapSI(:, :, cardinalPoints.left.SI_middle2+1:cardinalPoints.left.SI_max) = ...
    lungs_airwaybin(:, :, cardinalPoints.left.SI_middle2+1:cardinalPoints.left.SI_max) * 3;

rSegmentsMapSI(:, :, cardinalPoints.right.SI_min:cardinalPoints.right.SI_middle1) = ...
    lungs_airwaybin(:, :, cardinalPoints.right.SI_min:cardinalPoints.right.SI_middle1) * 4;
rSegmentsMapSI(:, :, cardinalPoints.right.SI_middle1+1:cardinalPoints.right.SI_middle2) = ...
    lungs_airwaybin(:, :, cardinalPoints.right.SI_middle1+1:cardinalPoints.right.SI_middle2) * 5;
rSegmentsMapSI(:, :, cardinalPoints.right.SI_middle2+1:cardinalPoints.right.SI_max) = ...
    lungs_airwaybin(:, :, cardinalPoints.right.SI_middle2+1:cardinalPoints.right.SI_max) * 6;

regionalMaps.SI = (lSegmentsMapSI .* uint8(lungs_L)) + (rSegmentsMapSI .* uint8(lungs_R));

fprintf('  Regional maps created\n\n');

end
