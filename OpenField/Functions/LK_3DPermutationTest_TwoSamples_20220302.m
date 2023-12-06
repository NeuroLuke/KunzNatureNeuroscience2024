function out = LK_3DPermutationTest_TwoSamples_20220302(dat)
%
% LK_3DPermutationTest_TwoSamples_20220302 performs cluster-based
% permutation testing between two three-dimensional datasets (e.g., time x
% frequency x channels).
%
% Note that the statistics are calculated across the third dimension of the
% input data matrix.
%
% Use as: out = LK_3DPermutationTest_TwoSamples_20220302(dat);
%
% Input is a structure with fields
%   mat1           	--> first data matrix
%   mat2            --> second data matrix
%   firstalpha      --> first-level alpha (e.g., 0.05)
%   secondalpha     --> second-level alpha (e..g, 0.05)
%   direction       --> expected direction of the effect (pos | neg)
%   numSurrogates   --> number of surrogates to create
%
% Output is a structure with fields
%   clusStat            --> sum statistic for each empirical cluster
%   maxClusStatSurro    --> maximum sum statistic for each surrogate round
%   clusStatRank        --> rank of each empirical cluster
%   clusStatP           --> p-value of each empirical cluster
%   bSigClusStat        --> boolean indicating the significant clusters
%   logIdxSigClus       --> logical indexing of all significant clusters
%   maxClusStat         --> maximum empirical sum statistic
%   logIdxMaxClus       --> logical indexing of the maximum cluster
%
% References
%   Oostenveld et al., Comput Intell Neurosci, 2011
%   Miller et al., Nat Commun, 2018
%
% Lukas Kunz, 2022

%% setup

% report
fprintf('\nCluster-based permutation testing of two samples against each other across the third dimension ...\n');

% specify the connectivity for identifying clusters
if isfield(dat, 'conn')
    myConn      = dat.conn;
else
    myConn      = 8;
end
fprintf('Connectivity is: %d.\n', myConn);

%% empirical statistic

% two-sample t-test between the two matrices
[~, pEmp, ~, statsEmp]  = ttest2(dat.mat1, dat.mat2, 'dim', 3);

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

% combine the data from the two different groups
mat12       = cat(3, dat.mat1, dat.mat2);
n1          = size(dat.mat1, 3); % size of the first sample
idx12       = transpose(1:size(mat12, 3)); % indices of the combined sample

% randomly assign indices to group 1
randIdx1    = nan(n1, dat.numSurrogates); % group-1-samples x surrogates
for iSurro = 1:dat.numSurrogates
    randIdx1(:, iSurro)     = datasample(idx12, n1, 'replace', false);
end

% create n maximum surrogate cluster statistics
maxClusStatSurro    = zeros(dat.numSurrogates, 1);
parfor iSurro = 1:dat.numSurrogates
    
    % randomly assign samples to the different groups
    thisRandIdx1    = randIdx1(:, iSurro);
    thisRandIdx2    = idx12(~ismember(idx12, thisRandIdx1));
        
    % two-sample t-test between the two matrices
    [~, pSurro, ~, statsSurro]  = ttest2(mat12(:, :, thisRandIdx1), mat12(:, :, thisRandIdx2), 'dim', 3);
    
    % identify surrogate clusters
    LSurro = [];
    NUMSurro = [];
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

% determine significance of each empirical cluster statistics
clusStatRank    = nan(numel(clusStat), 1); % rank of each empirical cluster
clusStatP       = nan(numel(clusStat), 1); % p-value of each empirical cluster
bSigClusStat    = false(numel(clusStat), 1); % boolean indicating the significant clusters
logIdxSigClus   = false(size(LEmp, 1), size(LEmp, 2)); % logical area indicating the significant time-frequency bins
for iClusStat = 1:numel(clusStat)
    
    % rank of empirical statistic within surrogate statistics
    if strcmp(dat.direction, 'pos')
        clusStatRank(iClusStat, 1)  = sum(clusStat(iClusStat, 1) > maxClusStatSurro) / numel(maxClusStatSurro); % significantly above the surrogates?
    elseif strcmp(dat.direction, 'neg')
        clusStatRank(iClusStat, 1)  = sum(clusStat(iClusStat, 1) < maxClusStatSurro) / numel(maxClusStatSurro); % significantly below the surrogates?
    end

    % p-value and significance
    clusStatP(iClusStat, 1)         = 1 - clusStatRank(iClusStat, 1);
    bSigClusStat(iClusStat, 1)      = clusStatP(iClusStat, 1) < dat.secondalpha;
    
    % logical indices of significant clusters
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
    maxClusStat = nan;
    maxClusIdx  = nan;
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
