function lowAttenuation = refineAirSeg(volumes, lungsProcessed, lowAttenuation, injurySeg, airways, vesselsRefined, params)
%REFINEAIRSEG Refine low-attenuation (air trapping) segmentation
%
%   lowAttenuation = refineAirSeg(volumes, lungsProcessed, lowAttenuation,
%       injurySeg, airways, vesselsRefined, params)
%
%   Performs two refinement passes on the low-attenuation mask:
%     1. Intensity cleanup + expansion via injury combing mask
%     2. Healthy tissue recovery: small unclassified fragments near AirSeg
%
%   Inputs:
%       volumes        - Structure with .filtered field (anisotropic diffusion CT)
%       lungsProcessed - Structure with .binary field (eroded lung mask)
%       lowAttenuation - Current low-attenuation segmentation mask
%       injurySeg      - Structure with .ggo, .dense, .injuryNewCombing
%       airways        - Airway segmentation mask (full tree)
%       vesselsRefined - Refined vessel segmentation mask
%       params         - Processing parameters (requires .airSeg)
%
%   Output:
%       lowAttenuation - Refined low-attenuation mask
%
%   See also: segmentTrachea, segmentInjuries

fprintf('  Refining low-attenuation segmentation...\n');

ap = params.airSeg;
ctFiltered = volumes.filtered;
lungsbin = lungsProcessed.binary;

%% Pass 1: Intensity cleanup + expansion via injury combing
% Remove non-air voxels from lowAttenuation
lowAttenuation(ctFiltered > ap.maxHU) = 0;

% Expand lowAttenuation into low-HU regions within injury combing mask
airSegExpand = ctFiltered < ap.expandHU;
airSegExpand(injurySeg.injuryNewCombing == 0) = 0;
airSegExpand = imreconstruct(lowAttenuation, airSegExpand, 8);
lowAttenuation = lowAttenuation | airSegExpand;

%% Pass 2: Healthy tissue -> AirSeg recovery
% Identify small unclassified lung fragments and absorb those near AirSeg
healthySeg = lungsbin;
healthySeg(injurySeg.ggo ~= 0) = 0;
healthySeg(injurySeg.dense ~= 0) = 0;
healthySeg(lowAttenuation ~= 0) = 0;
healthySeg(vesselsRefined ~= 0) = 0;
healthySeg(airways ~= 0) = 0;

healthy2Air = healthySeg & ~bwareaopen(healthySeg, ap.healthyMinVol, 8);
airSegDilated = imdilate(lowAttenuation, strel('disk', ap.dilationSize));
airSegFilled = imfill(lowAttenuation, 8, 'holes');
airSegDilated = airSegDilated & healthy2Air;
airSegFilled = airSegFilled & healthy2Air;
lowAttenuation = lowAttenuation | airSegDilated | airSegFilled;

% Ensure lowAttenuation does not overlap airways
lowAttenuation(airways ~= 0) = 0;

fprintf('  AirSeg refinement complete\n\n');

end
