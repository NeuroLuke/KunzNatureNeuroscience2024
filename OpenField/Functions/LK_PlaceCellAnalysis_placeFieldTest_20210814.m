function out = LK_PlaceCellAnalysis_placeFieldTest_20210814(cfg)
%
% LK_PlaceCellAnalysis_placeFieldTest_20210814 identifies a place field
% based on the behavioral data and the firing rates and estimates the
% strength of the place field.
%
% Steps:
% - Estimate the firing-rate map.
% - Smooth the firing-rate map (using a Gaussian kernel).
% - Identify candidate place fields.
% - Estimate the strength of each candidate place field and choose the best
%   candidate place field as the place field.
% - Test whether the firing rates within the place field are higher than
%   the firing rates outside the place field across time.
%
% Input:
%   FR                  --> firing rate per time point
%   xyBins              --> xy bin per time point
%   place               --> place settings including
%       facThreshPF     --> quantile for identifying candidate place fields
%                           (e.g., 0.75)
%       smoothSettings  --> smoothing settings including kernel type,
%                           kernel size, sigma
%
% Output:
%   FRMap       --> original firing-rate map
%   smFRMap     --> smoothed firing-rate map
%   threshPF    --> threshold for candidate place fields (Hz)
%   PF          --> place field
%   strengthPF  --> strength of the place field, estimated as the sum
%                   firing rate within the place field (Hz)
%   pPlace      --> p-value to evaluate whether the firing rates within the
%                   place field are consistently higher than outside the
%                   place field
%   tPlace      --> t-statistic to evaluate whether the firing rates within
%                   the place field are consistently higher than outside
%                   the place field
%
% Lukas Kunz, 2021

% initialize output
out = [];

%% sanity check

% ensure that there are no nans in the firing rates or the behavioral data
if sum(isnan(cfg.FR)) > 0 || sum(isnan(cfg.xyBins)) > 0
    error('There is a nan in the firing rates or in the behavior.');
end

%% firing rate per place bin

% estimate the firing-rate map
FRMap = nan(size(cfg.place.idxTemplate));
for iXY = min(cfg.xyBins):max(cfg.xyBins)
    FRMap(cfg.place.idxTemplate == iXY) = mean(cfg.FR(cfg.xyBins == iXY)); % average firing rate per xy bin
end

%% smoothing
% (following Burgess et al., Hippocampus, 2005)

% smoothing kernel
if size(cfg.place.smoothSettings, 2) == 3
    K   = fspecial(cfg.place.smoothSettings{1}, cfg.place.smoothSettings{2}, cfg.place.smoothSettings{3}); % kernel type, kernel size, sigma: 1 | 2 (default: 0.5)
elseif size(cfg.place.smoothSettings, 2) == 2
    K   = fspecial(cfg.place.smoothSettings{1}, cfg.place.smoothSettings{2});
end

% amount of smoothing per bin
denom                   = ones(size(FRMap));
denom(isnan(FRMap))     = 0;
denom                   = conv2(denom, K, 'same');

% smoothing of firing-rate map
smFRMap                 = FRMap;
smFRMap(isnan(FRMap))   = 0; % need to get rid of nans before smoothing
smFRMap                 = conv2(smFRMap, K, 'same');

% correction for amount of smoothing
smFRMap                 = smFRMap ./ denom;

% re-introduce nans
smFRMap(isnan(FRMap))   = nan;

% collect for output
out.FRMap               = FRMap; % firing-rate map
out.smFRMap             = smFRMap; % smoothed firing-rate map

%% place-field candidates

% firing-rate threshold for place fields
out.threshPF    = quantile(smFRMap(~isnan(smFRMap)), cfg.place.facThreshPF);

% candidate place fields
candidatePF     = smFRMap > out.threshPF;
[L, num]        = bwlabel(candidatePF, 4); % the "4" indicates the type of connectivity (here: pixels are connected if their edges touch)

% strength of each candidate place field
strengthCandidatePF = nan(num, 1);
for iL = 1:num
    strengthCandidatePF(iL) = sum(smFRMap(L == iL));
end

% if there are no candidate place fields, return
if isempty(strengthCandidatePF)
    
    % add essential output
    out.PF          = false(size(smFRMap));
    out.strengthPF  = nan;
    out.pPlace      = nan;
    out.tPlace      = nan;
    
    % skip the rest of this function
    return;
end

% place field = best candidate place field
[strengthPF, idxPF] = max(strengthCandidatePF);
PF                  = L == idxPF;

% xyBins associated with the place field
xyBinsPF            = cfg.place.idxTemplate(PF);

% collect for output
out.PF              = PF; % place field
out.strengthPF      = strengthPF; % strength of the maximum place field

%% test whether firing rates within the place field are higher than the other firing rates

% t-test between firing rates within versus outside the place field
bWithinPF           = ismember(cfg.xyBins, xyBinsPF);
[~, p, ~, stats]    = ttest2(cfg.FR(bWithinPF), cfg.FR(~bWithinPF));

% collect for output
out.pPlace          = p; % p-value for the comparison of place-field firing rates versus other firing rates
out.tPlace          = stats.tstat; % t-value for the comparison of place-field firing rates versus other firing rates
