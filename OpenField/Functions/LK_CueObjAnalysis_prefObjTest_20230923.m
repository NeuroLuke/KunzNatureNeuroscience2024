function out = LK_CueObjAnalysis_prefObjTest_20230923(cfg)
%
% LK_CueObjAnalysis_prefObjTest_20230923 computes a cell's strength of
% object tuning.
%
% Steps:
% - Estimate the firing rate per object;
% - Identify the preferred object;
% - Test whether the firing rates associated with the preferred object are
%   higher than the firing rates not associated with the preferred object.
%
% Input is a structure with fields
%   FR          --> average firing rate per object
%   Object      --> object names (range: 0-7)
%
% Output is a structure with fields
%   prefObjIdx  --> index of the preferred object (range: 1-8)
%   prefObjName --> name of the preferred object (range: 0-7)
%   objT        --> t-statistic when comparing the firing rates during the
%                   preferred vs. other objects
%
% Lukas Kunz, 2021

%% average firing rate per object

% object-specific firing rate and other statistics
objFR       = nan(numel(unique(cfg.Object)), 1);
objFRSEM    = nan(numel(unique(cfg.Object)), 1);
objFRMed    = nan(numel(unique(cfg.Object)), 1);
objFRMin    = nan(numel(unique(cfg.Object)), 1);
objFRMax    = nan(numel(unique(cfg.Object)), 1);
objFRp25    = nan(numel(unique(cfg.Object)), 1);
objFRp75    = nan(numel(unique(cfg.Object)), 1);
for iObj = min(cfg.Object):max(cfg.Object)
    % mean and SEM
    objFR(iObj + 1, 1)      = mean(cfg.FR(cfg.Object == iObj)); % mean
    objFRSEM(iObj + 1, 1)   = std(cfg.FR(cfg.Object == iObj)) / sqrt(sum(cfg.Object == iObj)); % SEM
    % other statistics
    objFRMed(iObj + 1, 1)   = median(cfg.FR(cfg.Object == iObj)); % median
    objFRMin(iObj + 1, 1)   = min(cfg.FR(cfg.Object == iObj)); % minimum
    objFRMax(iObj + 1, 1)   = max(cfg.FR(cfg.Object == iObj)); % maximum
    objFRp25(iObj + 1, 1)   = prctile(cfg.FR(cfg.Object == iObj), 25); % 25th percentile
    objFRp75(iObj + 1, 1)   = prctile(cfg.FR(cfg.Object == iObj), 75); % 75th percentile
end

% preferred object (i.e., the object with the highest firing rate)
[~, prefObjIdx] = max(objFR); % range: 1-8
prefObjName     = prefObjIdx - 1; % range: 0-7

% sanity check
if numel(objFR) ~= 8
    error('"objFR" does not have eight entries.');
end

%% test between firing rates of the preferred vs. unpreferred objects

% two-sided, two-sample t-test
[~, ~, ~, stats]    = ttest2(cfg.FR(cfg.Object == prefObjName), cfg.FR(cfg.Object ~= prefObjName));

% t-statistic
objT                = stats.tstat;

%% output

% create output
out             = [];
out.objFR       = objFR;
out.objFRSEM    = objFRSEM;
out.objFRMed    = objFRMed;
out.objFRMin    = objFRMin;
out.objFRMax    = objFRMax;
out.objFRp25    = objFRp25;
out.objFRp75    = objFRp75;
out.prefObjIdx  = prefObjIdx;
out.prefObjName = prefObjName;
out.objT        = objT;
