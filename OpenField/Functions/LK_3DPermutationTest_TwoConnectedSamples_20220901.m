function out = LK_3DPermutationTest_TwoConnectedSamples_20220901(dat)
%
% LK_3DPermutationTest_TwoConnectedSamples_20220901 performs cluster-based
% permutation testing between two connected 3D datasets (e.g., time x
% frequency x channels).
%
% Note that the statistics are calculated across the third dimension of the
% input data matrix.
%
% Note that this test works on "connected" datasets meaning that the first
% observation in dataset 1 is somehow related to the first observation in
% dataset 2. These two observations should thus not appear in the same
% surrogate dataset.
%
% Use as: out = LK_3DPermutationTest_TwoConnectedSamples_20220901(dat);
%
% Input is a structure with fields
%   mat1                --> first data matrix
%   mat2                --> second data matrix
%   firstalpha          --> first-level alpha (e.g., 0.05)
%   secondalpha         --> second-level alpha (e..g, 0.05)
%   direction           --> expected direction of the effect (pos | neg)
%   numSurrogates       --> number of surrogates to create
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
%   Maris & Oostenveld, J Neurosci Methods, 2007
%
% Lukas Kunz, 2022

%% setup

% report
fprintf('\nCluster-based permutation testing of two connected samples against each other across the third dimension ...\n');

% specify the connectivity for identifying clusters
if isfield(dat, 'conn')
    myConn      = dat.conn;
else
    myConn      = 8;
end
fprintf('Connectivity is: %d.\n', myConn);

% check whether there is enough data for the statistics, otherwise return
try
    ttest2(dat.mat1, dat.mat2, 'dim', 3);
catch
    % create output
    out                     = [];
    out.clusStat            = [];
    out.maxClusStatSurro    = zeros(dat.numSurrogates, 1);
    out.clusStatRank        = [];
    out.clusStatP           = [];
    out.bSigClusStat        = [];
    out.logIdxSigClus       = false(size(dat.mat1, 1), size(dat.mat1, 2));
    out.maxClusStat         = nan;
    out.logIdxMaxClus       = false(size(dat.mat1, 1), size(dat.mat1, 2));
    
    % exit the function
    fprintf(2, 'Data sample is not large enough.\n');
    return;
end

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

% randomly decide which observations to swap
n1          = size(dat.mat1, 3); % size of the first dataset
n2          = size(dat.mat2, 3); % size of the second dataset
bRandSwap   = rand(dat.numSurrogates, n1) > 0.5; % surrogates x samples (0 = keep; 1 = swap)

% sanity check to see whether both datasets have the same number of samples
if n1 ~= n2
    error('Number of samples per dataset is not identical.');
end

% create n maximum surrogate cluster statistics
maxClusStatSurro    = zeros(dat.numSurrogates, 1);
parfor iSurro = 1:dat.numSurrogates
    
    % surrogate datasets
    surroMat1   = dat.mat1;
    surroMat2   = dat.mat2;
    for iSample = 1:size(dat.mat1, 3)
        if bRandSwap(iSurro, iSample) == true % swap data
            surroMat1(:, :, iSample)    = dat.mat2(:, :, iSample);
            surroMat2(:, :, iSample)    = dat.mat1(:, :, iSample);
        end
    end
    
    % two-sample t-test between the two surrogate matrices
    [~, pSurro, ~, statsSurro]  = ttest2(surroMat1, surroMat2, 'dim', 3);
    
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
logIdxSigClus   = false(size(LEmp, 1), size(LEmp, 2)); % logical area indicating the significant bins
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
