function injurySeg = segmentInjuries(volumes, lungsProcessed, airwaysSeg, params)
%SEGMENTINJURIES Segment pathological lung tissue (low/high attenuation injuries)
%
%   injurySeg = segmentInjuries(volumes, lungsProcessed, airwaysSeg, params)
%
%   Segments pathological lung tissue including:
%   1. Air-filled pathological regions (emphysema, honeycombing, air trapping)
%   2. High attenuation abnormalities (GGO, reticulation/consolidations)
%
%   The algorithm adapts based on injury volume:
%     brightness mapping and dual-path analysis
%
%   Inputs:
%       volumes - Structure with fields:
%           .ct          - Original CT volume (int16, HU values)
%           .filtered    - Anisotropic diffusion filtered volume
%           .injury   - Initial injury mask from neural network
%       lungsProcessed - Structure with fields:
%           .binary       - Eroded lung mask
%       airwaysSeg - Structure with fields:
%           .airways     - Clean airway segmentation
%           .spurious    - Spurious airways
%       params - Processing parameters
%
%   Outputs:
%       injurySeg - Structure with fields:
%           .lowAttenuation        - Air-filled pathological regions (logical)
%           .highAttenuation       - Refined highAttenuation abnormalities segmentation (logical)
%           .dense           - Dense opacities (reticulation/consolidation) areas (logical)
%           .ggo                  - Ground glass opacity areas (logical)
%           .brightstructures     - Adaptive threshold result (for vessel seg)
%
%   Algorithm:
%       1. Segment air-filled pathological regions (HU <= -980)
%       2. Complete segmentation of High Attenuation abnormalities:
%          a. Create brightness map using adaptive Gaussian filtering
%          b. Identify pathological areas via histogram peak detection
%          c. Dual-path analysis (with/without air combing)
%          d. Combine and refine injury mask
%       3. Adaptive thresholding for high attenuation abnormalities
%
%   See also: preprocessLungs, refineAirways, segmentVessels

fprintf('  Segmenting pathological lung tissue...\n');

%% Initialize output structures
injurySeg = struct();

%% Get volume dimensions and references
ct = volumes.ct;
ct_filtered = volumes.filtered;
lungsbin = lungsProcessed.binary;
airways = airwaysSeg.airways;
spurious = airwaysSeg.spurious;

% Get initial injury mask
if isfield(volumes, 'injury')
    injury = volumes.injury;
else
    injury = false(size(lungsbin));
end

%% Step 1: Low-Attenuation (Air-filled) Pathological Tissue Segmentation
fprintf('    Segmenting air-filled pathological regions...\n');

% Segment air regions
AirSeg = zeros(size(ct));
AirSeg(ct_filtered <= params.thresholds.air) = 1;
AirSeg(ct_filtered > params.thresholds.air) = 0;
AirSeg(lungsbin == 0) = 0;

% Remove small isolated air regions
AirSeg = bwareaopen(AirSeg, params.air.minSize26, 26);

% Add spurious airways (honeycombing, cysts)
AirSeg = AirSeg | spurious;
% Exclude true airways
AirSeg = AirSeg & ~airways;

% Store air-filled pathology segmentation
injurySeg.lowAttenuation = logical(AirSeg);

% Create a new version of the air-filled path. regions for the honeycombing path
injuryLowAttenuation = bwareaopen(AirSeg, params.air.minSizeInjury26, 26);
injuryLowAttenuation = bwareaopen(injuryLowAttenuation, params.air.minSizeInjury8, 8);

%% Step 2: High-Attenuation Pathological Tissue Segmentation
fprintf('    Segmenting high attenuation abnormalities...\n');

%% Conditional adaptive thresholding to segment brightest structures
sigma = 0.3 * ((params.injury.kernelSize - 1) * 0.5 - 1) + 0.8;
ROI = lungsbin;
roi = double(ROI > 0);
ImSrc = double(imcomplement(ct)) .* (ROI > 0);

lungblur = imgaussfilt(ImSrc, sigma, 'FilterSize', params.injury.kernelSize);
roilungs = imgaussfilt(roi, sigma, 'FilterSize', params.injury.kernelSize);
smooth = lungblur ./ roilungs;
brightstructures = 255 .* (ImSrc <= (smooth - params.injury.adaptiveThresholdSep));
brightstructures = uint8(brightstructures .* (ROI > 0));

%% Create opacity maps for brightness mapping
% Standard opacity map
opacities = ct;
opacities(brightstructures > 0) = -1025;
opacities(airways > 0) = -1025;
opacities(opacities > params.thresholds.consolidationMin) = -1025;

% Opacity map with air combing (accounts for honeycombing)
opacities_combing = opacities;
opacities_combing(injuryLowAttenuation > 0) = -600;

%% brightness mapping with conditional gaussian filtering
% Standard brightness map
ImSrc = double(imcomplement(opacities)) .* (ROI > 0);
lungblur = double(imgaussfilt(ImSrc, sigma, 'FilterSize', params.injury.kernelSize));
roilungs = imgaussfilt(roi, sigma, 'FilterSize', params.injury.kernelSize);
smooth = (lungblur) ./ (roilungs);
smooth = imcomplement(smooth);
smooth(roi == 0) = -1025;

% Brightness map with air combing
ImSrc_combing = double(imcomplement(opacities_combing)) .* (ROI > 0);
lungblur_combing = double(imgaussfilt(ImSrc_combing, sigma, 'FilterSize', params.injury.kernelSize));
roilungs_combing = imgaussfilt(roi, sigma, 'FilterSize', params.injury.kernelSize);
smooth_combing = (lungblur_combing) ./ (roilungs_combing);
smooth_combing = imcomplement(smooth_combing);
smooth_combing(roi == 0) = -1025;

%% Calculate thresholds from histogram peaks
% Threshold for combing path
[thr_combing, ~] = calculateHistogramThreshold(smooth_combing, ...
    params.injury.peakOffset);
injury_new_combing = smooth_combing > thr_combing;

% Threshold for standard path
[thr_standard, ~] = calculateHistogramThreshold(smooth, ...
    params.injury.peakOffset);
injury_new = smooth > thr_standard;

%% Combine injury paths
injury_new = injury_new | injury_new_combing;

% Apply ROI mask
injury_new(roi == 0) = 0;

%% Find pathological voxels within the injury_new area
histo = ct_filtered;
histo(ct < params.injury.histoMinHU) = -1025;
histo(ct > params.injury.histoMaxHU) = -1025;
histo(brightstructures > 0) = -1025;
histo(airways > 0) = -1025;

[thr_o, ~] = calculateHistogramThreshold(histo, ...
    params.injury.vesselPeakOffset);

injury_new_o = ct_filtered > thr_o;

%% Final injury mask
injury = (injury | injury_new) & injury_new_o;
injury(airways > 0) = 0;

% Remove small isolated regions
injury = bwareaopen(injury, params.injury.minSize6, 6);
injury = bwareaopen(injury, params.injury.minSize4, 4);
  
%% Store outputs
injurySeg.highAttenuation = logical(injury);
injurySeg.brightstructures = logical(brightstructures);

%% Distinction between GGO and consolidation based on global thresholding
% Separate consolidation from GGO based on density
% Consolidation: denser areas (HU > -300)
% GGO: less dense areas (HU between -700 and -300)
%consolidation = injury & (ct > params.thresholds.consolidationMin);
%ggo = injury & (ct <= params.thresholds.consolidationMin) & ...
%    (ct > params.thresholds.ggoMin);
%
%injurySeg.dense = consolidation;
%injurySeg.ggo = ggo;

%% Distinction between GGO and consolidation based on adaptive thresholding
sigma = 0.3 * ((params.injury.kernelSize - 1) * 0.5 - 1) + 0.8;
ROI = lungsbin;
roi = double(ROI > 0);
ImSrc = double(imcomplement(ct)) .* (ROI > 0);

lungblur = imgaussfilt(ImSrc,sigma,'FilterSize',params.injury.kernelSize);
roilungs = imgaussfilt(roi,sigma,'FilterSize',params.injury.kernelSize);
smooth = lungblur./roilungs;

dense = 255.*(ImSrc <=(smooth-params.injury.adaptiveThresholdSep));
dense = uint8(dense.*(ROI>0));
dense(injury==0)=0;
dense = bwareaopen(dense,6,8); % REVISE
dense = bwareaopen(dense,35,6); % REVISE
ggo = injury & ~dense;

injurySeg.dense = dense;
injurySeg.ggo = ggo;
    
fprintf('    Injury segmentation complete\n');

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [threshold, peakLocation] = calculateHistogramThreshold(image, offset)
%CALCULATEHISTOGRAMTHRESHOLD Calculate threshold from histogram peak
%
%   Finds the second-highest peak in the histogram and adds an offset.
%   This identifies the typical intensity of the tissue type.

    % Compute histogram
    [counts, binLocations] = imhist(int16(image), 65535);

    % Find highest peak (usually background)
    [~, idx] = max(counts);
    counts(idx) = 0;

    % Find second highest peak (tissue of interest)
    [~, idx] = max(counts);
    peakLocation = binLocations(idx);

    % Calculate threshold with offset
    threshold = peakLocation + offset;
end
