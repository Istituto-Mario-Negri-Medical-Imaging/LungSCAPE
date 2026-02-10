function [vesselsRefined, largeVessels] = refineVessels(vesselsSeg, volumes, lungsProcessed, airwaysSeg, injurySeg, distanceMasks, params)
%REFINEVESSELS Advanced vessel post-processing and refinement
%
%   [vesselsRefined, largeVessels] = refineVessels(vesselsSeg, volumes, lungsProcessed, airwaysSeg, injurySeg, distanceMasks, params)
%
%   Performs vessel refinement complementary to segmentVessels:
%   1. Adaptive thresholding to recover vessels in healthy tissue
%   2. Final false positive removal by anisotropy and diameter
%   3. Large vessel reconstruction for structures missed by Jerman filter
%
%   Note: Vessel recovery near airways, fibrotic band removal, and
%   reticulation removal are now handled in segmentVessels.
%
%   Inputs:
%       vesselsSeg     - Structure from segmentVessels with fields:
%                        .final - Final vessel segmentation
%                        .vesselness - Vesselness response map
%       volumes        - Structure with CT volumes (.filtered, .ct)
%       lungsProcessed - Structure with processed lung masks (.binary, .lungs)
%       airwaysSeg     - Structure with airway segmentation (.airways)
%       injurySeg      - Structure with injury segmentation (.dense, .consolidation)
%       distanceMasks  - Distance-based lung masks (.close, .midfar, .far)
%       params         - Processing parameters
%
%   Outputs:
%       vesselsRefined - Final refined vessel segmentation (logical)
%       largeVessels   - Large vessels recovered from intensity (logical)
%
%   See also: segmentVessels, classifyVesselsByDiameter

fprintf('  Post-processing vessel segmentation...\n');

%% Initialize working copy
vessels = vesselsSeg.final;

%% Adaptive thresholding for complete vessels in healthy tissue
fprintf('    Applying adaptive thresholding in healthy tissue...\n');

% Build masks structure for adaptive thresholding
masks.lungs = lungsProcessed.binary;
masks.airways = airwaysSeg.airways;
masks.far = distanceMasks.far;

completeVessels = adaptiveVesselThreshold(volumes.filtered, masks, params);

% Vesselness filtering on adaptive threshold result
fprintf('    Applying vesselness filter to complete vessels...\n');
completeVesselsFrangi = vesselness3D(completeVessels, ...
    [0.5, 1], ...
    [params.voxel.px; params.voxel.py; params.voxel.vz], ...
    0.5, true);

completeVessels = completeVesselsFrangi > params.vessels.threshold.complete;

% Check for fibrotic bands (near consolidation)
if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
    fibroticCheck = imreconstruct(uint8(injurySeg.dense), ...
        uint8(completeVessels), 26);

    % Apply stricter vesselness threshold near pathology
    fibroticCheck(vesselsSeg.vesselness > 0.8) = 0;
    completeVessels = completeVessels & ~fibroticCheck;
end

% Combine with vessels from segmentVessels
vessels = vessels | completeVessels;

%% Final false positive removal
fprintf('    Final false positive removal...\n');

% Segregated objects with large diameter and low anisotropy
vesselsCheck = vessels & ~bwareaopen(vessels, 500);
[diameterImage, ~] = computeVesselDiameter(vesselsCheck, 0);

diameterImageLarge = diameterImage > 3;
diameterImageLarge = bwareaopen(diameterImageLarge, 4);

vesselsCheck = imreconstruct(diameterImageLarge, vessels, 26);

% Filter by anisotropy
[~, vesselsCheckDiscard] = filterByAnisotropy(vesselsCheck, ...
    params.vessels.anisotropy.final);

finalDiscard0 = vesselsCheckDiscard;

% Segregated objects with large diameter and low volume
vesselsCheck = vessels & ~bwareaopen(vessels, 70);
[diameterImage, ~] = computeVesselDiameter(vesselsCheck, 0);

diameterImageLarge = diameterImage > 3;
finalDiscard1 = imreconstruct(diameterImageLarge, vessels, 6);

finalDiscard = finalDiscard0 | finalDiscard1;
vessels = vessels & ~finalDiscard;

%% Distally segregated objects with large diameter
lastSegregated = imreconstruct(distanceMasks.close, vessels, 26);
lastSegregated = vessels & ~lastSegregated;

[diameterImage, ~] = computeVesselDiameter(lastSegregated, 0);

diameterImageMain = diameterImage >= 3;
diameterImageMain = bwareaopen(diameterImageMain, 3, 26);
lastSegregated = imreconstruct(diameterImageMain, lastSegregated, 26);

% Remove if in distal regions (lung apices)
lungsBsApx = lungsProcessed.binary & ~bwareaopen(lungsProcessed.binary, 7000, 8);
lastSegregatedRemove = imreconstruct(lungsBsApx, lastSegregated, 26);
lastSegregated = lastSegregated & ~lastSegregatedRemove;

if isfield(injurySeg, 'dense') && any(injurySeg.dense(:))
    lastSegregated = imreconstruct(injurySeg.dense, lastSegregated, 26);
end

vessels = vessels & ~lastSegregated;

%% Add large vessels that may have been missed by the Jerman filter
fprintf('    Segmenting large vessels from intensity...\n');
largeVessels = addLargeVessels(volumes, lungsProcessed, vessels, airwaysSeg, distanceMasks, params);

% Combine large vessels with refined vessels
vessels = vessels | largeVessels;

%% Final size filtering
vessels = bwareaopen(vessels, 30, 6);

%% Output
vesselsRefined = vessels;

fprintf('    Vessel post-processing complete\n');

end


%% Helper: Adaptive vessel thresholding in healthy tissue
function completeVessels = adaptiveVesselThreshold(volumeFiltered, masks, params)
    % Apply conditional adaptive thresholding to segment vessels in healthy tissue

    C = 180;  % Threshold constant
    ksize = 15;  % Gaussian kernel size
    sigma = 0.3 * ((ksize - 1) * 0.5 - 1) + 0.8;

    ROI = masks.lungs;

    % Exclude airways from processing
    [airwaysDist, ~] = bwdist(single(masks.airways));
    outROI = airwaysDist < 3 & airwaysDist > 0;

    completeVessels = uint8(zeros(size(masks.lungs)));

    % Process slice by slice for memory efficiency
    numSlices = size(volumeFiltered, 3);

    for k = 1:numSlices
        ROIslice = double(ROI(:, :, k));
        src = double(imcomplement(volumeFiltered(:, :, k))) .* (ROIslice > 0);
        roi = double(ones(size(ROIslice)));
        roi = roi .* (ROIslice > 0);

        % Adaptive thresholding
        blur = imgaussfilt(src, sigma, 'FilterSize', ksize);
        roi = imgaussfilt(roi, sigma, 'FilterSize', ksize);
        smooth = blur ./ roi;
        output = 255 .* (src <= (smooth - C));
        healthyVess = uint8(output .* (ROIslice > 0));

        completeVessels(:, :, k) = uint8(healthyVess);
    end

    % Clean up borders
    if isfield(masks, 'far')
        lungs_R = masks.lungs == 1;
        lungs_L = masks.lungs == 2;

        lungs_R_bsapx = lungs_R & ~bwareaopen(lungs_R, 7000, 8);
        lungs_L_bsapx = lungs_L & ~bwareaopen(lungs_L, 7000, 8);
        lungs_bsapx = lungs_R_bsapx | lungs_L_bsapx;

        cleanBorders0 = imdilate(masks.far, strel('disk', 5));
        cleanBorders1 = imdilate(masks.far, strel('disk', 2));
        cleanBorders0(lungs_bsapx ~= 0) = 0;
        cleanBorders1(lungs_bsapx == 0) = 0;
        cleanBorders = cleanBorders0 | cleanBorders1;
        cleanBorders(masks.lungs == 0) = 0;

        completeVessels(cleanBorders ~= 0) = 0;
    end

    completeVessels = logical(completeVessels);
end

%% Helper function: Add large vessels
function largeVessels = addLargeVessels(volumes, lungsProcessed, vesselsUpdate, airwaysSeg, distanceMasks, params)
    %ADDLARGEVESSELS Segment large vessels using intensity thresholding
    %   Large vessels may not be well-captured by the Jerman filter due to
    %   their size. This function uses intensity-based segmentation to
    %   recover them.

    % Threshold based on HU values
    largeVessels0 = volumes.filtered > params.vessels.largeVessels.MinHU & ...
                    volumes.filtered < params.vessels.largeVessels.MaxHU;
    largeVessels0(lungsProcessed.binary == 0) = 0;
    largeVessels0 = largeVessels0 & ~airwaysSeg.airways;

    % Compute diameter
    [diameterImage, ~] = computeVesselDiameter(largeVessels0, 0);

    % Keep only large diameter structures
    diameterImageLarge = diameterImage > params.vessels.largeVessels.diameter;
    diameterImageLarge = bwareaopen(diameterImageLarge, params.vessels.largeVessels.minLength, 26);

    % Expand around large diameter skeleton
    [largeDist, ~] = bwdist(single(diameterImageLarge));
    largeDistMask = largeDist < 3 & largeDist > 0;
    largeVessels = largeVessels0 & (largeDistMask | diameterImageLarge);
    largeVessels = imreconstruct(diameterImageLarge, largeVessels, 8);
    largeVessels = imreconstruct(distanceMasks.close, largeVessels, 26);

    % Clean up
    largeVessels = largeVessels & ~airwaysSeg.airways;
    largeVessels = bwareaopen(largeVessels, params.vessels.largeVessels.minVol, 6);

    % Remove already segmented vessels
    largeVessels = largeVessels & ~vesselsUpdate;
end