%==========================================================================
% This script identifies object cells during the cue period.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                	= [];
param.myRNG            	= 111;
param.numSurrogates    	= 1001; % number of surrogates
param.subjectFileName   = 'subjectdata_20210616.mat';
param.clusterName     	= 'LK_Cluster4Analysis_20210521';
param.timeRes           = 0.01; % time resolution (s)
param.trialTimeEdges    = -1:param.timeRes:3; % time bins for time-resolved firing rates
param.trialTimeCenters  = movmean(param.trialTimeEdges, 2, 'endpoints', 'discard');
param.baselineTimeEdges = [-1, 0]; % time edges for the baseline correction
param.CBPTTimeEdges     = [0, 2]; % time edges for the cluster-based permutation test (on time-resolved firing rates)
param.smoothType        = 'gaussian'; % kernel for smoothing the time-resolved firing rates
param.smoothDur         = 0.5; % smoothing duration (s)
param.smoothFac         = param.smoothDur / param.timeRes + 1; % factor for smoothing the time-resolved firing rates

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.spikes        = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.fieldtrip     = 'E:\fieldtrip\fieldtrip-20210614';
paths.placeCells    = 'E:\OpenField\PlaceCellAnalysis_20210511\20211116_Loc25x25_facThreshPF0.75\';
paths.save          = strcat('E:\OpenField\CueObjAnalysis_20210812\20211116\');
mkdir(paths.save);
addpath(paths.fieldtrip);
ft_defaults;
ft_warning off;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings
save(strcat(paths.save, 'settings'));

%% preallocations

% main results for each cell
allRes      = [];

% bookkeeping
exByWC      = []; % excluded by wave-clus
exByVisInsp = []; % excluded by visual inspection

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % skip if no microwire recordings
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        fprintf('\nSkipping subject %s because no microwires.\n', subjects{iSub});
        continue;
    end
    
    % reset rng and report
    rng(param.myRNG);
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    
    % subject-specific save path
    subjSavePath    = strcat(paths.save, subjects{iSub}, '\');
    mkdir(subjSavePath);
    
    %% behavior
    
    % report
    fprintf('\nProcessing of behavioral data (using microtime).\n');
    
    % load trial information in MICROTIME
    trialInfo   = load(strcat(paths.beh, subjects{iSub}, '\trialInfo.mat'), 'trialInfoMicrotime');
    trialInfo 	= trialInfo.trialInfoMicrotime;
    
    %% channels with microwire data
    
    % available channels
    chans   = LK_GetChans(fullfile(paths.spikes, subjects{iSub}, 'chan*'));
    
    % report
    fprintf('\nProcessing of neural data (%d channels from "%s" to "%s").\n', size(chans, 1), chans(1).name, chans(end).name);
    
    %% loop through wires
    for iChan = 1:size(chans, 1)
        
        %% brain region of this microwire channel
        
        % report
        fprintf('\tChannel: %s.\n', chans(iChan).name);
        
        % original channel index
        chanIdx     = str2double(replace(chans(iChan).name, lettersPattern, ''));
        
        % subject information
        s           = load(strcat(paths.subjects, subjects{iSub}, filesep, param.subjectFileName));
        logIdx  	= any(cell2mat(s.subjectdata.micro2macro(:, 1)) == chanIdx, 2);
        chanRegion  = s.subjectdata.micro2macro{logIdx, 3}; % brain region
        
        %% wave-clus output
        
        % load wave-clus output
        wcFile  = fullfile(chans(iChan).folder, chans(iChan).name, 'times_datacut.mat');
        if exist(wcFile, 'file') > 0
            times_datacut   = load(wcFile);
        else
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC          = cat(1, exByWC, [iSub, chanIdx]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on inspection)
        c4a = load(strcat(chans(iChan).folder, filesep, chans(iChan).name, filesep, param.clusterName, '.mat'));
        
        %% loop through clusters
        for iClus = 1:max(times_datacut.cluster_class(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
            
            % continue if you decided that this cluster is insufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'no')
                fprintf('\t\t- You decided not to analyze this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, chanIdx, iClus]); % bookkeeping
                continue;
            end
            
            % data from this cluster
            thisCluster = times_datacut.cluster_class(times_datacut.cluster_class(:, 1) == iClus, :); % cluster-class, MICROTIME (MSEC)
            thisSpike   = times_datacut.spikes(times_datacut.cluster_class(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            % convert spike times to seconds
            thisCluster(:, 2)   = thisCluster(:, 2) ./ 1000;
            
            %% firing rates during the cue period
            
            % onsets and ends of the cue period
            onsEnds     = [trialInfo.Cue, trialInfo.Retrieval];
            durations   = diff(onsEnds, [], 2);
            
            % sanity check
            if any(isnan(durations))
                error('A duration is a nan.');
            end
            
            % time-resolved trial times
            trOnsEnds   = trialInfo.Cue + param.trialTimeEdges; % relative to cue onset
            trDurations = round(diff(trOnsEnds, [], 2), 6);
            
            % number of spikes
            numSpikes   = nan(size(onsEnds, 1), 1);
            trNumSpikes = nan(size(trOnsEnds, 1), size(trOnsEnds, 2) - 1);
            for iTrial = 1:size(onsEnds, 1)
                numSpikes(iTrial, 1)    = histcounts(thisCluster(:, 2), onsEnds(iTrial, :)); % total number of spikes
                trNumSpikes(iTrial, :)  = histcounts(thisCluster(:, 2), trOnsEnds(iTrial, :)); % time-resolved number of spikes
            end
            
            % firing rates
            FR      = numSpikes ./ durations; % overall firing rate
            trFR    = trNumSpikes ./ trDurations; % time-resolved firing rates
            
            %% test whether a specific object elicits significantly higher firing rates than the other objects
            
            % empirical test
            cfg             = [];
            cfg.FR          = FR;
            cfg.Object      = trialInfo.Object;
            prefObjTest     = LK_CueObjAnalysis_prefObjTest_20230923(cfg);
            
            % indices to randomly shuffle the firing rates
            [~, randIdx]    = sort(rand(size(FR, 1), param.numSurrogates));
            
            % surrogates
            objTSurro       = nan(param.numSurrogates, 1); % surrogate t-values
            for iSurro = 1:param.numSurrogates
                
                % shuffled firing rates relative to the presented objects
                FRSurro                 = FR(randIdx(:, iSurro));
                
                % surrogate test
                cfg                     = [];
                cfg.FR                  = FRSurro;
                cfg.Object              = trialInfo.Object;
                prefObjTestSurro        = LK_CueObjAnalysis_prefObjTest_20230923(cfg);
                
                % surrogate t-statistic
                objTSurro(iSurro, 1)    = prefObjTestSurro.objT;
            end
            
            % empirical p-value
            objTRank	= sum(prefObjTest.objT > objTSurro) / sum(~isnan(objTSurro));
            objTP       = 1 - objTRank;
            
            %% cluster-based permutation test on the time-resolved firing rates
            
            % smooth the time-resolved firing rates over time
            trFR        = smoothdata(trFR, 2, param.smoothType, param.smoothFac);
            
            % baseline correction after smoothing
            logIdx      = param.trialTimeCenters >= min(param.baselineTimeEdges) & param.trialTimeCenters <= max(param.baselineTimeEdges);
            trFR        = trFR - mean(trFR(:, logIdx), 2);
            
            % split into preferred and unpreferred trials (based on the
            % previously identified preferred object)
            trFRPref    = trFR(trialInfo.Object == prefObjTest.prefObjName, :);
            trFRUnpref  = trFR(trialInfo.Object ~= prefObjTest.prefObjName, :);
            
            % re-organize the data for fieldtrip
            timelock1       = [];
            timelock1.label = {'neuron'};
            for iTrial = 1:size(trFRPref, 1)
                timelock1.time{1, iTrial}    = param.trialTimeCenters;
                timelock1.trial{1, iTrial}   = trFRPref(iTrial, :);
            end
            timelock2       = [];
            timelock2.label = {'neuron'};
            for iTrial = 1:size(trFRUnpref, 1)
                timelock2.time{1, iTrial}    = param.trialTimeCenters;
                timelock2.trial{1, iTrial}   = trFRUnpref(iTrial, :);
            end
            
            % fieldtrip configuration
            cfg                     = [];
            cfg.channel             = 'all';
            cfg.latency             = [0, 2]; % time period for the cluster-based permutation test
            cfg.method              = 'montecarlo';
            cfg.statistic           = 'indepsamplesT'; % first-level test
            cfg.correctm            = 'cluster';
            cfg.clusteralpha        = 0.05; % first-level alpha
            cfg.clusterstatistic    = 'maxsum'; % test statistic
            cfg.neighbours          = [];
            cfg.tail                = 1; % H1: stronger activation during preferred than during unpreferred trials
            cfg.alpha               = 0.05; % second-level alpha
            cfg.correcttail         = 'alpha';
            cfg.numrandomization    = param.numSurrogates;
            cfg.design              = [ones(size(timelock1.trial)), 2 .* ones(size(timelock2.trial))]; % design matrix
            cfg.ivar                = 1; % row of the design matrix with the independent variable
            
            % fieldtrip estimation
            objCBPT                 = ft_timelockstatistics(cfg, timelock1, timelock2);
            
            % identify the largest significant cluster
            if isfield(objCBPT, 'posclusters')
                objCBPTSigP         = min(cell2mat({objCBPT.posclusters.prob})); % minimum p-value
                objCBPTSigClus      = objCBPT.mask & (objCBPT.posclusterslabelmat == 1); % cluster "1" is the largest cluster by fieldtrip definition
                objCBPTSigClusTime  = objCBPT.time(objCBPTSigClus);
            else
                objCBPTSigP         = nan;
                objCBPTSigClus      = false(size(objCBPT.time));
                objCBPTSigClusTime  = objCBPT.time(objCBPTSigClus);
            end
                        
            %% figure: object tuning
            
            % create figure
            dt                      = [];
            dt.visible              = 'off';
            dt.param                = param;
            dt.trialInfo            = trialInfo(:, {'Object', 'xCorrect', 'yCorrect'}); % relevant trial information
            dt.prefObjTest          = prefObjTest; % results from testing for a preferred object
            dt.objTP                = objTP; % p-value
            dt.trFR                 = trFR; % time-resolved firing rates
            dt.trNumSpikes          = trNumSpikes; % time-resolved number of spikes
            dt.objCBPTSigClusTime   = objCBPTSigClusTime; % time period of the significant cluster
            dt.objCBPTSigP          = objCBPTSigP; % p-value from the cluster-based permutation test
            dt.thisSpike            = thisSpike; % waveforms
            dt.sr                   = times_datacut.par.sr; % sampling rate
            dt.nspk                 = sum(numSpikes); % total number of spikes during the cue period
            dt.chanRegion           = chanRegion; % channel region
            dt.idx                  = [iSub, chanIdx, iClus]; % indexing
            outFigure               = LK_PlotCueObjTuning_20230923(dt);
            
            % save figure
            LK_print(outFigure.f, strcat(paths.save, subjects{iSub}, '_', chans(iChan).name, '_clus', num2str(iClus), '_obj_20230923'), '-dtiff', '-r300');
            
            % save figure data
            if all(dt.idx == [24, 5, 1], 2)
                save(strcat(paths.save, 'Fig_5a_upperLeft'), 'dt');
            elseif all(dt.idx == [32, 34, 3], 2)
                save(strcat(paths.save, 'Fig_5a_upperRight'), 'dt');
            elseif all(dt.idx == [34, 2, 1], 2)
                save(strcat(paths.save, 'Fig_5a_lowerLeft'), 'dt');
            elseif all(dt.idx == [40, 3, 1], 2)
                save(strcat(paths.save, 'Fig_5a_lowerRight'), 'dt');
            end

            %% collect information for this unit
            
            % basics
            unitRes                     = [];
            unitRes.idx                 = [iSub, chanIdx, iClus];
            unitRes.subjectName         = subjects{iSub, 1};
            unitRes.chanRegion          = chanRegion;
            
            % detailed info
            unitRes.numSpikes           = numSpikes;
            unitRes.durations           = durations;
            unitRes.FR                  = FR;
            unitRes.trNumSpikes         = trNumSpikes;
            unitRes.trDurations         = trDurations;
            unitRes.trFR                = trFR;
            
            % evaluation via average firing rates
            unitRes.objFR               = prefObjTest.objFR;
            unitRes.prefObjIdx          = prefObjTest.prefObjIdx;
            unitRes.prefObjName         = prefObjTest.prefObjName;
            unitRes.objT                = prefObjTest.objT;
            unitRes.objTSurro           = objTSurro;
            unitRes.objTRank            = objTRank;
            unitRes.objTP               = objTP;
            unitRes.mUnprefAndPref      = outFigure.m; % mean firing rate during unpreferred and preferred trials
            unitRes.steUnprefAndPref    = outFigure.ste; % standard error during unpreferred and preferred trials
            
            % evaluation via time-resolved firing rates
            unitRes.objCBPTSigP         = objCBPTSigP;
            unitRes.objCBPTSigClusTime  = objCBPTSigClusTime;
            
            % collapse across units
            allRes                      = cat(1, allRes, unitRes);
            
            %% close all open figures
            close all;
            toc
        end
    end
end

%% results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r = load(strcat(paths.save, 'results.mat'));

%% general

% remove fieldtrip from the path
rmpath(genpath(paths.fieldtrip));
clc;

% cellular indices
allIdx          = cell2mat({r.allRes.idx}');
allSubjectName  = {r.allRes.subjectName}';
fprintf('Number of cells: %d.\n', size(allIdx, 1));

% brain regions
allRegions      = {r.allRes.chanRegion}';

%% object cells

% identify cells that fulfill the criterion regarding average firing rates
allObjT         = cell2mat({r.allRes.objT}');
allObjTSurro    = cell2mat({r.allRes.objTSurro})';
allObjTRank     = sum(allObjT > allObjTSurro, 2) ./ sum(~isnan(allObjTSurro), 2);
bAvObjCell   	= allObjTRank > 0.95;

% identify cells that fulfill the criterion regarding time-resolved firing
% rates (evaluated using fieldtrip's cluster-based permutation test)
allObjCBPTSigP  = cell2mat({r.allRes.objCBPTSigP}');
bTrObjCell      = allObjCBPTSigP < 0.05;

% identify object cells based on both criteria
bObjCell        = bAvObjCell & bTrObjCell;
fprintf('Number of object cells fulfilling both criteria: %d (%.3f%%, P = %.3f).\n', sum(bObjCell), ...
    100 * sum(bObjCell) / numel(bObjCell), myBinomTest(sum(bObjCell), numel(bObjCell), 0.05));

% save cell classification
cellClassification              = cellstr(repmat('NonObjectCell', size(allIdx, 1), 1));
cellClassification(bObjCell)    = {'ObjectCell'};
cellClassification              = table(cellClassification, allSubjectName, allIdx(:, 1), allIdx(:, 2), allIdx(:, 3), ...
    'VariableNames', {'cellType', 'subjectName', 'session', 'channel', 'cluster'});
save(strcat(paths.save, 'cellClassification'), 'cellClassification');

%% object cells per session and subject

% number of object cells per session
numObjCellsPerSess  = nan(size(r.subjects, 1), 1);
numCellsPerSess     = nan(size(r.subjects, 1), 1);
for iSess = 1:size(r.subjects, 1)
    numObjCellsPerSess(iSess, 1)    = sum(bObjCell & allIdx(:, 1) == iSess);
    numCellsPerSess(iSess, 1)       = sum(allIdx(:, 1) == iSess);
end

% sessions with microwires
bMicrowire  = strcmp(r.subjects(:, 2), 'Microwire');
fprintf('Number of microwire sessions with at least one object cell: %d (out of %d).\n', sum(bMicrowire & numObjCellsPerSess > 0), sum(bMicrowire));
LK_ReportMeanAndSEM_20220322('Number of object cells per session', numObjCellsPerSess(bMicrowire));

%% object cells per region

% report
fprintf('\nAnalysis of object cells in different brain regions.\n');

% all regions and unique brain regions
uniqueRegions   = unique(allRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% percentage of object cells per unique region
for iR = 1:size(uniqueRegions, 1)
    
    % cells from this regions
    bThisRegion             = strcmp(allRegions, uniqueRegions{iR});
    numCellsPerReg(iR, 1)   = sum(bThisRegion);
    fprintf('Region: %s. Number of cells: %d.\n', uniqueRegions{iR}, sum(bThisRegion));
    
    % binomial test
    pBinom  = myBinomTest(sum(bObjCell & bThisRegion), sum(bThisRegion), 0.05);
    
    % report
    fprintf('\tNumber of object cells: %d (%.3f%%, P = %.3f).\n', sum(bObjCell & bThisRegion), ...
        100 * sum(bObjCell & bThisRegion) / sum(bThisRegion), pBinom);
end

% restrict unique regions to those with enough units
uniqueRegions   = uniqueRegions(numCellsPerReg > 30);

% figure: percentage of object cells per unique region
dt                      = [];
dt.uniqueRegions        = uniqueRegions;
dt.allRegions           = allRegions;
dt.bCell                = bObjCell;
dt.ylabel               = 'Object cells (%)';
dt.BonferroniCorrection = true;
f                       = LK_PlotCellsPerRegion(dt);

% save figure
LK_print(f, strcat(paths.save, 'objCellsPerRegion_20220329'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_5b'), 'dt');

%% figure: significant time periods from the cluster-based permutation tests

% all significant time periods
allSigTime  = nan(size(r.allRes, 1), size(param.trialTimeCenters, 2));
for iClus = 1:size(r.allRes, 1)
    if ~isempty(r.allRes(iClus).objCBPTSigClusTime)
        sigTime                 = round(r.param.trialTimeCenters, 4) >= round(min(r.allRes(iClus).objCBPTSigClusTime), 4) & ...
            round(r.param.trialTimeCenters, 4) <= round(max(r.allRes(iClus).objCBPTSigClusTime), 4);
        allSigTime(iClus, :)    = sigTime;
    end
end

% report peak time
[maxVal, maxIdx]    = max(mean(allSigTime(bObjCell, :)));
fprintf('\nDistribution of significant time periods of object cells. Max = %.3f at %.3f s.\n', maxVal, r.param.trialTimeCenters(maxIdx));

% create figure showing all significant time periods
dt          = [];
dt.x        = r.param.trialTimeCenters;
dt.m        = mean(allSigTime(bObjCell, :)); % average across object cells
dt.ste      = LK_ste(allSigTime(bObjCell, :));
dt.param    = param;
f = LK_PlotSigTimeForPrefVSUnprefObjects_20231008(dt);

% save figure
LK_print(f, strcat(paths.save, 'objCell_allSigTime_20220329'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_5c'), 'dt');

%% figure: average firing rates for preferred and unpreferred objects

% firing rates during preferred and unpreferred objects for all cells
allUnpref   = cell2mat(cellfun(@(x) x{1}, {r.allRes.mUnprefAndPref}', 'UniformOutput', false));
allPref     = cell2mat(cellfun(@(x) x{2}, {r.allRes.mUnprefAndPref}', 'UniformOutput', false));

% means and standard errors
groups  = {'Unpreferred', [0.5, 0.5, 0.5]; 'Preferred', rgb('orange')};
x       = r.param.trialTimeCenters;
m       = {mean(allUnpref(bObjCell, :), 1), mean(allPref(bObjCell, :), 1)}; % unpreferred, preferred
ste     = {LK_ste(allUnpref(bObjCell, :)), LK_ste(allPref(bObjCell, :))}; % unpreferred, preferred

% report information about the maximum
[maxVal, maxIdx]    = max(m{2});
fprintf('\nAverage firing rate of object cells in response to preferred objects. Max = %.3f Hz at %.3f s.\n', maxVal, r.param.trialTimeCenters(maxIdx));

% create figure
dt          = [];
dt.groups   = groups;
dt.x        = x;
dt.m        = m;
dt.ste      = ste;
dt.param    = param;
f = LK_PlotAvFRForPrefAndUnprefObjects_20231008(dt);

% save figure
LK_print(f, strcat(paths.save, 'objCells_prefVSUnpref_20220329'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_5d'), 'dt');

%% property: temporal stability of the tuning curves

% identify time points during the cue period
bCuePeriod  = r.param.trialTimeCenters >= min(r.param.CBPTTimeEdges) & r.param.trialTimeCenters <= max(r.param.CBPTTimeEdges);

% loop through clusters and estimate temporal stability
temporalStabPref    = nan(size(r.allRes, 1), 1);
temporalStabUnpref  = nan(size(r.allRes, 1), 1);
for iClus = 1:size(r.allRes, 1)
    
    % load trial information
    trialInfo       = load(strcat(paths.beh, subjects{r.allRes(iClus).idx(1)}, '\trialInfo.mat'), 'trialInfoMicrotime');
    trialInfo       = trialInfo.trialInfoMicrotime;
    
    % trials with the preferred object
    bPrefObj        = trialInfo.Object == r.allRes(iClus).prefObjName;
    
    % time-resolved firing rates during trials with the preferred object or
    % the unpreferred object
    trFRPrefObj     = r.allRes(iClus).trFR(bPrefObj, :);
    trFRUnprefObj   = r.allRes(iClus).trFR(~bPrefObj, :);
    
    % correlation between the firing rates from the 1st vs. 2nd data half,
    % preferred objects
    b1stHalf        = 1:size(trFRPrefObj, 1) < median(1:size(trFRPrefObj, 1));
    m1              = mean(trFRPrefObj(b1stHalf, :), 1);
    m2              = mean(trFRPrefObj(~b1stHalf, :), 1);
    temporalStabPref(iClus, 1)      = corr(m1(bCuePeriod)', m2(bCuePeriod)');
    
    % correlation between the firing rates from the 1st vs. 2nd data half,
    % unpreferred objects
    b1stHalf        = 1:size(trFRUnprefObj, 1) < median(1:size(trFRUnprefObj, 1));
    m1              = mean(trFRUnprefObj(b1stHalf, :), 1);
    m2              = mean(trFRUnprefObj(~b1stHalf, :), 1);
    temporalStabUnpref(iClus, 1)    = corr(m1(bCuePeriod)', m2(bCuePeriod)');
end

%% create figure

% t-test of temporal stabilities against zero
[~, p, ~, stats]    = ttest(temporalStabPref(bObjCell));
fprintf('Temporal stability of object cells against zero: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
LK_ReportMeanAndSEM_20220322('Mean correlation between tuning curves', temporalStabPref(bObjCell));

% create figure: histogram of temporal stabilities
dt              = [];
dt.temporalStab = temporalStabPref(bObjCell);
f = LK_PlotTemporalStabilityOfCueObjCells_20231008(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'objCells_temporalStab_20220329'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_5e'), 'dt');
