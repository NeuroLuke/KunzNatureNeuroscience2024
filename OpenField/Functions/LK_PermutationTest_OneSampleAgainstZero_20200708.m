function out = LK_PermutationTest_OneSampleAgainstZero_20200708(cfg)
%
% LK_PermutationTest_OneSampleAgainstZero_20200708 performs cluster-based
% permutation testing against zero for two-dimensional data (e.g., n x p).
%
% Note that the statistic is calculated across the first dimension of the
% input data matrix.
%
% Use as: out = LK_PermutationTest_OneSampleAgainstZero_20200708(cfg);
%
% Input is a structure with fields
%   mat             --> data matrix (n x p)
% 	alpha           --> first-level and second-level alpha (e.g., 0.05)
%   direction       --> hypothesized direction of the effect (pos | neg)
%   numSurrogates   --> number of surrogates to create
%
% Output is a structure with fields
%   clusStat            --> empirical cluster statistics
%   maxClusStatSurro    --> maximum surrogate cluster statistics
%   clusStatRanks       --> ranks of empirical cluster statistics
%   clusStatP           --> p-values of empirical cluster statistics
%   bSigClusStat        --> significant clusters
%   logIdxSigClus       --> logical indexing of significant clusters
%   maxClusStat         --> maximum empirical cluster statistic
%   logIdxMaxClus       --> logical indexing of largest cluster
%
% References
%   Oostenveld et al., Comput Intell Neurosci, 2011
%   Miller et al., Nat Commun, 2018
% 
% Lukas Kunz, 2021

% report
fprintf('... Cluster-based permutation testing of one sample against zero ...\n');

%% empirical statistic

% one-sample t-test against zero (across the first dimension)
[~, pEmp, ~, statsEmp]  = ttest(cfg.mat, 0);

% identify empirical clusters
if strcmp(cfg.direction, 'pos')
    [LEmp, NUMEmp]  = bwlabel(pEmp < cfg.alpha & statsEmp.tstat > 0); % significantly above zero
elseif strcmp(cfg.direction, 'neg')
    [LEmp, NUMEmp]  = bwlabel(pEmp < cfg.alpha & statsEmp.tstat < 0); % significantly below zero
end
clusStat    = nan(NUMEmp, 1);
for iNUM = 1:NUMEmp
    clusStat(iNUM, 1)   = sum(statsEmp.tstat(LEmp == iNUM));
end

%% surrogate statistics

% create n maximum surrogate cluster statistics
maxClusStatSurro    = zeros(cfg.numSurrogates, 1);
for iSurro = 1:cfg.numSurrogates
    
    % randomly flip the sign of the data to create surrogate data
    randSign                    = datasample([1; -1], size(cfg.mat, 1), 'replace', true);
    [~, pSurro, ~, statsSurro]  = ttest(cfg.mat .* repmat(randSign, 1, size(cfg.mat, 2)), 0);
    
    % identify surrogate clusters
    if strcmp(cfg.direction, 'pos')
        [LSurro, NUMSurro]      = bwlabel(pSurro < cfg.alpha & statsSurro.tstat > 0); % significantly above zero
    elseif strcmp(cfg.direction, 'neg')
        [LSurro, NUMSurro]      = bwlabel(pSurro < cfg.alpha & statsSurro.tstat < 0); % significantly below zero
    end
    clusStatSurro               = nan(NUMSurro, 1);
    for iNUM = 1:NUMSurro
        clusStatSurro(iNUM, 1)  = sum(statsSurro.tstat(LSurro == iNUM));
    end
    
    % maximum cluster statistic
    if ~isempty(clusStatSurro)
        if strcmp(cfg.direction, 'pos')
            maxClusStatSurro(iSurro, 1) = max(clusStatSurro); % maximum positive surrogate cluster statistic
        elseif strcmp(cfg.direction, 'neg')
            maxClusStatSurro(iSurro, 1) = min(clusStatSurro); % maximum negative surrogate cluster statistic
        end
    end
end

%% significance

% determine significance of each empirical cluster statistic
clusStatRanks   = nan(numel(clusStat), 1); % rank of each empirical cluster
clusStatP       = nan(numel(clusStat), 1); % p-value of each empirical cluster
bSigClusStat    = false(numel(clusStat), 1); % boolean indicating the significant clusters
logIdxSigClus   = false(1, numel(LEmp)); % boolean indicating the significant bins
for iClusStat = 1:numel(clusStat)
    
    % rank of empirical statistic within surrogate statistics
    if strcmp(cfg.direction, 'pos')
        clusStatRanks(iClusStat, 1) = sum(clusStat(iClusStat, 1) > maxClusStatSurro) / numel(maxClusStatSurro); % significantly above the surrogates?
    elseif strcmp(cfg.direction, 'neg')
        clusStatRanks(iClusStat, 1) = sum(clusStat(iClusStat, 1) < maxClusStatSurro) / numel(maxClusStatSurro); % significantly below the surrogates?
    end
    clusStatP(iClusStat, 1)     = 1 - clusStatRanks(iClusStat, 1); % p-value = 1 - rank
    bSigClusStat(iClusStat, 1)  = clusStatP(iClusStat, 1) < cfg.alpha; % whether the p-value is below the alpha level
    
    % columns of significant clusters
    if bSigClusStat(iClusStat, 1) == true
        logIdxSigClus(LEmp == iClusStat)    = true;
    end
end

% maximum positive/negative cluster statistic
if strcmp(cfg.direction, 'pos')
    [maxClusStat, maxClusIdx]   = max(clusStat);
elseif strcmp(cfg.direction, 'neg')
    [maxClusStat, maxClusIdx]   = min(clusStat);
end
if isempty(maxClusStat)
    maxClusStat                 = nan;
    maxClusIdx                  = nan;
end

% identify location of maximum cluster
logIdxMaxClus   = LEmp == maxClusIdx;
if sum(logIdxMaxClus) > 0
    fprintf('\tThe maximum cluster statistic is: cluster-t = %.3f, cluster-p = %.3f.\n', maxClusStat, clusStatP(maxClusIdx));
    fprintf('\tThe maximum cluster lasts from index %d to %d.\n', find(logIdxMaxClus, 1, 'first'), find(logIdxMaxClus, 1, 'last'));
end

% report minimum and maximum degrees of freedom of the empirical data
fprintf('\tThe minimum and maximum degrees of freedom of the empirical data is: %d and %d, respectively.\n', min(statsEmp.df), max(statsEmp.df));

%% output

% create output
out                     = [];
out.clusStat            = clusStat; % empirical cluster statistics
out.maxClusStatSurro    = maxClusStatSurro; % maximum surrogate cluster statistics
out.clusStatRanks       = clusStatRanks; % ranks of empirical cluster statistics
out.clusStatP           = clusStatP; % p-values of empirical cluster statistics
out.bSigClusStat        = bSigClusStat; % significant clusters
out.logIdxSigClus       = logIdxSigClus; % logical indexing of significant clusters
out.maxClusStat         = maxClusStat; % maximum empirical cluster statistic
out.logIdxMaxClus       = logIdxMaxClus; % logical indexing of largest cluster
