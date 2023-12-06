function out = LK_3DPermutationTest_OneSampleAgainstZero_20220302(dat)
%
% LK_3DPermutationTest_OneSampleAgainstZero_20220302 performs cluster-based
% permutation testing against zero for three-dimensional data (e.g., time x
% frequency x channels).
%
% Note that the statistics are calculated across the third dimension of the
% input data matrix.
%
% Use as: out = LK_3DPermutationTest_OneSampleAgainstZero_20220302(dat);
%
% Input is a structure with fields
%   mat             --> data matrix (e.g., t x f x c)
%   firstalpha      --> first-level alpha (e.g., 0.05 if it is a one-sided test)
%   secondalpha     --> second-level alpha (e.g., 0.05 if it is a one-sided test)
%   direction       --> expected direction of the effect (pos | neg)
%   conn            --> connectivity of significant bins (4 | 8)
%   numSurrogates   --> number of surrogates to create
%
% Output is a structure with fields
%   clusStat            --> sum statistic for each empirical cluster (1D vector)
%   maxClusStatSurro    --> maximum sum statistic for each surrogate round
%   clusStatRank        --> rank of each empirical cluster within the surrogates (1D vector)
%   clusStatP           --> p-value of each empirical cluster (1D vector)
%   bSigClusStat        --> boolean indicating the significant clusters (1D vector)
%   logIdxSigClus       --> logical indexing of all significant clusters (2D matrix)
%   maxClusStat         --> maximum empirical sum statistic (1 value)
%   logIdxMaxClus       --> logical indexing of the maximum cluster (2D matrix)
%
% References
%   Oostenveld et al., Comput Intell Neurosci, 2011
%   Miller et al., Nat Commun, 2018
%
% Lukas Kunz, 2022

%% setup

% report
fprintf('\nCluster-based permutation testing of one sample against zero across the third dimension ...\n');

% specify the connectivity for identifying clusters
if isfield(dat, 'conn')
    myConn      = dat.conn; % custom connectivity
else
    myConn      = 8; % default connectivity
end
fprintf('Connectivity is: %d.\n', myConn);

% check whether there is enough data for the statistics, otherwise return
try
    ttest(dat.mat, 0, 'dim', 3);
catch
    % create output
    out                     = [];
    out.clusStat            = [];
    out.maxClusStatSurro    = zeros(dat.numSurrogates, 1);
    out.clusStatRank        = [];
    out.clusStatP           = [];
    out.bSigClusStat        = [];
    out.logIdxSigClus       = false(size(dat.mat, 1), size(dat.mat, 2));
    out.maxClusStat         = nan;
    out.logIdxMaxClus       = false(size(dat.mat, 1), size(dat.mat, 2));
    
    % exit the function
    fprintf(2, 'Data sample is not large enough.\n');
    return;
end

%% empirical statistic

% t-test against 0 across the third dimension
[~, pEmp, ~, statsEmp]  = ttest(dat.mat, 0, 'dim', 3);

% identify empirical clusters
if strcmp(dat.direction, 'pos')
    [LEmp, NUMEmp]  = bwlabel(pEmp < dat.firstalpha & statsEmp.tstat > 0, myConn); % significantly above zero
elseif strcmp(dat.direction, 'neg')
    [LEmp, NUMEmp]  = bwlabel(pEmp < dat.firstalpha & statsEmp.tstat < 0, myConn); % significantly below zero
end

% sum up t-values within clusters
clusStat    = nan(NUMEmp, 1);
for iNUM = 1:NUMEmp
    clusStat(iNUM, 1)   = sum(statsEmp.tstat(LEmp == iNUM));
end

%% surrogate statistics

% create random signs (outside the parfor loop)
randSign    = nan(size(dat.mat, 3), dat.numSurrogates); % samples x surrogate rounds
for iSurro = 1:dat.numSurrogates
    randSign(:, iSurro)     = datasample([1; -1], size(dat.mat, 3), 'replace', true);
end

% create n maximum surrogate cluster statistics
maxClusStatSurro    = zeros(dat.numSurrogates, 1);
parfor iSurro = 1:dat.numSurrogates
    
    % randomly flip the sign of the data to create surrogate data
    thisRandSign    = randSign(:, iSurro);
    thisRandSign    = permute(thisRandSign, [3, 2, 1]);
    thisRandSign    = repmat(thisRandSign, size(dat.mat, 1), size(dat.mat, 2), 1);
    
    % t-test against 0 across the third dimension using randomly flipped
    % data
    [~, pSurro, ~, statsSurro]  = ttest(dat.mat .* thisRandSign, 0, 'dim', 3);
    
    % identify surrogate clusters
    LSurro      = [];
    NUMSurro    = [];
    if strcmp(dat.direction, 'pos')
        [LSurro, NUMSurro]      = bwlabel(pSurro < dat.firstalpha & statsSurro.tstat > 0, myConn); % significantly above zero
    elseif strcmp(dat.direction, 'neg')
        [LSurro, NUMSurro]      = bwlabel(pSurro < dat.firstalpha & statsSurro.tstat < 0, myConn); % significantly below zero
    end
    
    % sum up t-values within clusters
    clusStatSurro               = nan(NUMSurro, 1);
    for iNUM = 1:NUMSurro
        clusStatSurro(iNUM, 1)  = sum(statsSurro.tstat(LSurro == iNUM));
    end
    
    % maximum cluster statistic
    if ~isempty(clusStatSurro)
        if strcmp(dat.direction, 'pos')
            maxClusStatSurro(iSurro, 1) = max(clusStatSurro);
        elseif strcmp(dat.direction, 'neg')
            maxClusStatSurro(iSurro, 1) = min(clusStatSurro);
        end
    end
end

%% significance

% determine significance of each empirical cluster statistic
clusStatRank    = nan(numel(clusStat), 1); % rank of each empirical cluster
clusStatP       = nan(numel(clusStat), 1); % p-value of each empirical cluster
bSigClusStat    = false(numel(clusStat), 1); % boolean indicating the significant clusters
logIdxSigClus   = false(size(LEmp, 1), size(LEmp, 2)); % logical area indicating the significant bins
for iClusStat = 1:numel(clusStat)
    
    % rank of empirical statistic within surrogate statistics
    if strcmp(dat.direction, 'pos')
        clusStatRank(iClusStat, 1)  = sum(clusStat(iClusStat, 1) > maxClusStatSurro) / numel(maxClusStatSurro); % significantly above the surrogates?
    elseif strcmp(dat.direction, 'neg')
        clusStatRank(iClusStat, 1)  = sum(clusStat(iClusStat, 1) < maxClusStatSurro) / numel(maxClusStatSurro); % significantly below the surrogates?
    end

    % p-value and significance
    clusStatP(iClusStat, 1)     = 1 - clusStatRank(iClusStat, 1);
    bSigClusStat(iClusStat, 1)  = clusStatP(iClusStat, 1) < dat.secondalpha;
    
    % bins of significant clusters
    if bSigClusStat(iClusStat, 1) == true
        logIdxSigClus(LEmp == iClusStat)    = true;
    end
end

% maximum cluster statistic
if strcmp(dat.direction, 'pos')
    [maxClusStat, maxClusIdx]   = max(clusStat);
elseif strcmp(dat.direction, 'neg')
    [maxClusStat, maxClusIdx]   = min(clusStat);
end
if isempty(maxClusStat)
    maxClusStat             = nan;
    maxClusIdx              = nan;
end

% identify location of maximum cluster
logIdxMaxClus   = LEmp == maxClusIdx;
if sum(logIdxMaxClus(:)) > 0
    fprintf('\tThe maximum cluster statistic is: cluster-t = %.3f, cluster-p = %.3f.\n', maxClusStat, clusStatP(maxClusIdx));
end

% create output
out                     = [];
out.clusStat            = clusStat; % empirical cluster statistics
out.maxClusStatSurro    = maxClusStatSurro; % maximum surrogate cluster statistics
out.clusStatRank        = clusStatRank; % ranks of empirical cluster statistics
out.clusStatP           = clusStatP; % p-values of empirical cluster statistics
out.bSigClusStat        = bSigClusStat; % significant clusters
out.logIdxSigClus       = logIdxSigClus; % logical indexing of significant clusters
out.maxClusStat         = maxClusStat; % maximum empirical cluster statistic
out.logIdxMaxClus       = logIdxMaxClus; % logical indexing of largest cluster
