%==========================================================================
% This script examines single-neuron activity during hippocampal ripples.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;

% settings
param                           = [];
param.myRNG                     = 555;
%- subjects
param.subject.fileName          = 'subjectdata_20210616.mat';
%- microwires
param.micro.c4aFileName         = 'LK_Cluster4Analysis_20210521.mat';
param.micro.ccmFileName         = 'cluster_class_macrotime.mat';
%- preprocessing
param.preproc.type              = 'FBS'; % filtered, bipolar, select
%- ripples
param.ripple.type               = 'ripples';
param.ripple.fileName           = strcat(param.ripple.type, '.mat');
param.ripple.bpfreq             = [80, 140];
param.ripple.threshMinFac       = 2;
param.ripple.timepointType      = 'peak';
param.ripple.twoi4BaselineERP   = [-3, 3]; % time window for baseline correction of the ERP
%- analysis
param.ana.timeRes               = 1/500; % time resolution (in seconds)
param.ana.timeEdges             = -3:param.ana.timeRes:3; % time-bin edges around ripples
param.ana.timeCenters           = movmean(param.ana.timeEdges, 2, 'endpoints', 'discard');
param.ana.smoothType            = 'Gaussian';
param.ana.smoothTime            = 0.2; % smoothing window (s)
param.ana.smoothFac             = param.ana.smoothTime / param.ana.timeRes + 1; % smoothing factor
param.ana.numSurrogates         = 1001; % number of surrogates

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.spikes        = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.ripple        = strcat('E:\OpenField\RipplesMacro_20210614\20220320_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '\');
paths.save          = strcat('E:\OpenField\NeuronalActivityDuringRipples_20210721\20220612_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '_', ...
    param.ripple.type, '_', param.ripple.timepointType, 'Time', filesep);
paths.functions     = 'E:\OpenField\Functions\';
paths.fieldtrip     = 'E:\fieldtrip\fieldtrip-20210614\';
mkdir(paths.save);

% add functions and fieldtrip
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects   	= subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings
save(strcat(paths.save, 'settings'));

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByInspection  = [];
exByNoRipples   = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % skip if it is not a microwire subject
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        fprintf('\nSkipping subject %s because no microwires.\n', subjects{iSub});
        continue;
    end
    
    % report
    fprintf('\n\n==========================================================\n');
    fprintf('SUBJECT: %s.\n', subjects{iSub});
    
    %% behavioral data (to be able to restrict the ripples to specific conditions)
    
    % load trial information in MACROTIME
    trialInfo   = load(strcat(paths.beh, subjects{iSub}, '\trialInfo.mat'));
    trialInfo 	= trialInfo.trialInfoMacrotime;
    fprintf('\tNumber of cue periods: %d.\n', sum(~isnan(trialInfo.Cue)));
    
    % start and end time of the relevant data segment
    expStartEndTime    = [min(trialInfo.Cue), max(trialInfo.Cue)];
    
    %% ripple data
    
    % channels from which ripples were extracted
    rippleChans = LK_dir(strcat(paths.ripple, subjects{iSub}, '\*'));
    
    % assemble the ripple data from all ripple channels into one variable
    % for later use
    rippleData  = repmat(struct(), size(rippleChans, 1), 1);
    for iRippleChan = 1:size(rippleChans, 1)
        
        % load ripples and extract their information
        r           = load(fullfile(rippleChans(iRippleChan).folder, rippleChans(iRippleChan).name, param.ripple.fileName));
        ripples     = r.ripples.ripples;
        eegRipples  = r.ripples.eegRipples;
        
        % determine which ripples you can keep (only keep ripples that are
        % fully between the first and the last cue)
        bRipples2Keep	= cell2mat({ripples.startTime}') >= min(expStartEndTime) & cell2mat({ripples.endTime}') <= max(expStartEndTime);
        
        % remove ripples starting before the first cue or ending after the
        % last cue
        ripples       	= ripples(bRipples2Keep);
        cfg             = [];
        cfg.trials      = bRipples2Keep;
        eegRipples      = ft_selectdata(cfg, eegRipples);
        
        % collect the information from all ripple channels
        rippleData(iRippleChan).ripples     = ripples;
        rippleData(iRippleChan).eegRipples  = eegRipples;
        
        % report
        fprintf('\tRipple channel: %s. Keeping === %d === of %d ripples.\n', rippleChans(iRippleChan).name, sum(bRipples2Keep), size(bRipples2Keep, 1));
    end
    
    %% spike data
    
    % microwire channels
    spikeChans  = LK_GetChans(strcat(paths.spikes, subjects{iSub}, '\chan*'));
    
    %% spike data

    % loop through microwire channels
    for iSpikeChan = 1:size(spikeChans, 1)
        
        % report
        fprintf('\tSpike channel: %s.\n', spikeChans(iSpikeChan).name);
        
        % channel index
        spikeChanIdx    = str2double(replace(spikeChans(iSpikeChan).name, lettersPattern, ''));
        
        % load wave-clus output with spike times in MACROTIME
        c4aFile = fullfile(spikeChans(iSpikeChan).folder, spikeChans(iSpikeChan).name, param.micro.c4aFileName);
        ccmFile = fullfile(spikeChans(iSpikeChan).folder, spikeChans(iSpikeChan).name, param.micro.ccmFileName); % cluster-class, microtime (msec), macrotime (sec)
        
        % load decision which clusters to analyze
        c4a = load(c4aFile);
        
        % if the decision file is not empty, load the spike-time
        % information
        if ~isempty(c4a.Cluster4Analysis)
            ccm = load(ccmFile);
        else
            fprintf('\t... skipping this microwire channel, because the decision-file is empty.\n');
            continue;
        end
        
        % brain region and hemisphere of this wire
        s                   = load(fullfile(paths.subjects, subjects{iSub}, param.subject.fileName));
        logIdx              = any(cell2mat(s.subjectdata.micro2macro(:, 1)) == spikeChanIdx, 2);
        spikeChanRegion     = s.subjectdata.micro2macro{logIdx, 3}; % brain region
        spikeChanHemisphere = s.subjectdata.micro2macro{logIdx, 4}; % hemisphere
        
        %% loop through clusters
        for iClus = 1:max(ccm.cluster_class_macrotime(:, 1))
            
            % skip if visual inspection decided against this cluster
            if ~strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'yes')
                exByInspection  = cat(1, exByInspection, [iSub, spikeChanIdx, iClus]);
                continue;
            end
            
            % data from this cluster
            thisCCM         = ccm.cluster_class_macrotime(ccm.cluster_class_macrotime(:, 1) == iClus, :); % cluster class, microtime (MSEC), macrotime (SEC)
            
            % raw and smoothed firing rate per time point
            macroTimeEdges  = min(expStartEndTime):param.ana.timeRes:max(expStartEndTime);
            macroTime       = round(movmean(macroTimeEdges, 2, 'endpoints', 'discard'), 4);
            rawFR         	= histcounts(thisCCM(:, 3), macroTimeEdges) ./ param.ana.timeRes;
            smFR         	= smoothdata(rawFR, 2, param.ana.smoothType, param.ana.smoothFac);
            
            % compute the mean and the standard deviation of the raw firing
            % rates for later z-scoring
            meanRawFR   	= mean(rawFR, 2);
            stdRawFR     	= std(rawFR, [], 2);
            
            % convert the smoothed firing rates into z-scores
            smFRZ        	= zscore(smFR);
            
            % sanity check regarding spike timing
            if size(thisCCM, 1) > 0 && (min(thisCCM(:, 3)) < min(expStartEndTime) || max(thisCCM(:, 3)) > max(expStartEndTime))
                error('Problem with the spike times in "thisCCM".');
            end
            
            %% loop through ripple channels
            % you obtain results for each combination of a neuron and its
            % simultaneously recorded ripple channels
            
            % number of spikes and duration per time bin
            for iRippleChan = 1:size(rippleChans, 1)
                
                % macro-timepoints of all ripples of this channel
                rippleTimepoints   	= cell2mat({rippleData(iRippleChan).ripples.peakTime}');
                
                % skip if there are no ripples
                if size(rippleTimepoints, 1) < 1
                    fprintf('\tSkipping %s (%s) because there are no ripples.\n', subjects{iSub}, rippleChans(iRippleChan).name);
                    exByNoRipples   = cat(1, exByNoRipples, [iSub, spikeChanIdx, iClus, iRippleChan]);
                    continue;
                end
                
                %% firing rates during ripples
                
                % firing rates around the ripple peak
                rippleFR    = nan(size(rippleData(iRippleChan).ripples, 1), numel(param.ana.timeCenters)); % ripples x time-centers
                rippleFRZ   = nan(size(rippleData(iRippleChan).ripples, 1), numel(param.ana.timeCenters));
                rippleSmFRZ = nan(size(rippleData(iRippleChan).ripples, 1), numel(param.ana.timeCenters));
                for iRipple = 1:size(rippleData(iRippleChan).ripples, 1)
                    
                    % time bins around this ripple
                    thisRippleTimeEdges     = rippleTimepoints(iRipple, 1) + param.ana.timeEdges;
                    thisRippleTime          = round(movmean(thisRippleTimeEdges, 2, 'endpoints', 'discard'), 4);
                    
                    % raw firing rates during ripples
                    rippleFR(iRipple, :)    = histcounts(thisCCM(:, 3), thisRippleTimeEdges) ./ param.ana.timeRes;
                    
                    % z-scored raw firing rates
                    rippleFRZ(iRipple, :)   = (rippleFR(iRipple, :) - meanRawFR) ./ stdRawFR;
                    
                    % z-scored smoothed firing rates (needs a different
                    % indexing method)
                    [~, tmpIdx]             = min(abs(macroTime - rippleTimepoints(iRipple, 1)));
                    macroTimeIdx            = (tmpIdx - numel(param.ana.timeCenters) / 2):(tmpIdx + numel(param.ana.timeCenters) / 2 - 1);
                    if min(macroTimeIdx) < 1 || max(macroTimeIdx) > numel(smFRZ)
                        warning('Skipping ripple %d because indices go beyond the data.', iRipple);
                    else
                        rippleSmFRZ(iRipple, :) = smFRZ(macroTimeIdx);
                    end
                end
                
                % average the ripple-related firing rates across ripples
                meanRippleFR    = mean(rippleFR, 1, 'omitnan');
                meanRippleFRZ   = mean(rippleFRZ, 1, 'omitnan');
                meanRippleSmFRZ = mean(rippleSmFRZ, 1, 'omitnan');
                
                %% ripple data
                
                % baseline correction
                cfg             = [];
                cfg.baseline    = param.ripple.twoi4BaselineERP;
                rippleERP       = ft_timelockbaseline(cfg, rippleData(iRippleChan).eegRipples);
                
                % averaging
                cfg             = [];
                rippleERP       = ft_timelockanalysis(cfg, rippleERP);
                
                %% results
                
                % results from this unit
                thisRes                         = [];
                thisRes.idx                     = [iSub, spikeChanIdx, iClus];
                thisRes.subjectName             = subjects{iSub};
                thisRes.spikeChanRegion         = spikeChanRegion;
                thisRes.spikeChanHemisphere     = spikeChanHemisphere;
                
                % ripple information
                thisRes.rippleChanName          = rippleChans(iRippleChan).name;
                thisRes.rippleTimepoints        = rippleTimepoints;
                thisRes.rippleERPTime           = rippleERP.time;
                thisRes.rippleERPAvg            = rippleERP.avg; % average ripple
                
                % firing rates during ripples
                thisRes.rippleFR                = rippleFR; % raw firing rates during ripples
                thisRes.meanRippleFR            = meanRippleFR; % average raw firing rate during ripples
                thisRes.rippleFRZ               = rippleFRZ; % z-scored raw firing rates during ripples
                thisRes.meanRippleFRZ           = meanRippleFRZ; % average z-scored raw firing rate during ripples
                thisRes.rippleSmFRZ             = rippleSmFRZ; % z-scored smoothed firing rates during all ripples
                thisRes.meanRippleSmFRZ         = meanRippleSmFRZ; % average z-scored smoothed firing rate during ripples
                
                % average firing rate in this unit
                thisRes.overallFR               = size(thisCCM, 1) / range(expStartEndTime); % overall firing rate
                
                % behavior
                thisRes.trialInfo               = trialInfo;
                
                % assemble across different units
                allRes                          = cat(1, allRes, thisRes);
                
                %% close open figures
                close all;
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r   = load(strcat(paths.save, 'results.mat'));

%% general information

% all and unique cellular indices
allIdx      = cell2mat({r.allRes.idx}');
uniqueIdx   = unique(allIdx, 'rows');

% labels of all hippocampal ripple channels
allRippleChanName       = {r.allRes.rippleChanName}';
uniqueRippleChanName    = unique(allRippleChanName);

% all and unique brain regions of all wires
allSpikeChanRegion      = {r.allRes.spikeChanRegion}';
uniqueSpikeChanRegion   = unique(allSpikeChanRegion);

% hemispheres of all wires
allSpikeChanHemisphere  = {r.allRes.spikeChanHemisphere}';

% information whether neurons and ripples were recorded from the same
% hemisphere
bSameHemisphere = false(size(r.allRes, 1), 1);
for iRes = 1:size(r.allRes, 1)
    if (strcmp(allSpikeChanHemisphere{iRes}, 'left') && contains(allRippleChanName{iRes}, 'L')) || ...
            (strcmp(allSpikeChanHemisphere{iRes}, 'right') && contains(allRippleChanName{iRes}, 'R'))
        bSameHemisphere(iRes, 1)  = true;
    end
end

% all ripple averages
allRippleERPAvg     = cell2mat({r.allRes.rippleERPAvg}');
allRippleERPTime    = cell2mat({r.allRes.rippleERPTime}');

% report
fprintf('\n==== Results.\n');
fprintf('Number of unit-ripple combinations: %d.\n', size(allIdx, 1));
fprintf('Number of unique units: %d.\n', size(uniqueIdx, 1));
fprintf('Number of units recorded from the same hemisphere as their hippocampal ripples: %d.\n', sum(bSameHemisphere));
fprintf('Names of the ripple channels: %s.\n', uniqueRippleChanName{:});

%% normalized firing rates during ripples, separately for different brain regions

% reset rng
rng(r.param.myRNG);

% update the parameters for this analysis
thisParam               = [];
thisParam.regionGroups  = {'All regions'; 'AMY'; 'EC'; 'HC'; 'PHC'; 'TP'}; % brain regions to examine
thisParam.xLim          = [-1, 1];
thisParam.yLim          = [-0.02, 0.07];
thisParam.yLabel        = {'Firing rate', '(z-scored)'};

% all normalized smoothed firing rates during ripples
allMeanRippleSmFRZ      = cell2mat({r.allRes.meanRippleSmFRZ}');

% loop through brain regions
for iRG = 1:size(thisParam.regionGroups, 1)
    
    % select neurons from this brain region
    if strcmp(thisParam.regionGroups{iRG}, 'All regions')
        bMask 	= true(size(allSpikeChanRegion, 1), 1);
    else
        bMask   = strcmp(allSpikeChanRegion, thisParam.regionGroups{iRG});
    end
    fprintf('\nRegion: %s. Number of combinations: %d.\n', thisParam.regionGroups{iRG}, sum(bMask));
    
    %% statistics: firing rates against 0
    
    % test z-scored firing rates against zero
    cfg               	= [];
    cfg.mat             = allMeanRippleSmFRZ(bMask, :);
    cfg.alpha        	= 0.05;
    cfg.direction      	= 'pos';
    cfg.numSurrogates 	= param.ana.numSurrogates;
    permTest          	= LK_PermutationTest_OneSampleAgainstZero_20200708(cfg);
    fprintf('\tTime of maximum cluster: %.3f to %.3f s.\n', min(r.param.ana.timeCenters(permTest.logIdxSigClus & permTest.logIdxMaxClus)), max(r.param.ana.timeCenters(permTest.logIdxSigClus & permTest.logIdxMaxClus)));
    
    %% figure
    
    % create figure
    dt                  = [];
    dt.tRipple          = mean(allRippleERPTime(bMask, :), 1);
    dt.mRipple          = mean(allRippleERPAvg(bMask, :), 1, 'omitnan');
    dt.regionGroup      = thisParam.regionGroups{iRG};
    dt.nSame            = sum(bMask & bSameHemisphere);
    dt.nDiff            = sum(bMask & ~bSameHemisphere);
    dt.t                = r.param.ana.timeCenters;
    dt.mSame            = mean(allMeanRippleSmFRZ(bMask & bSameHemisphere, :), 1, 'omitnan');
    dt.steSame          = LK_ste(allMeanRippleSmFRZ(bMask & bSameHemisphere, :));
    dt.mDiff            = mean(allMeanRippleSmFRZ(bMask & ~bSameHemisphere, :), 1, 'omitnan');
    dt.steDiff          = LK_ste(allMeanRippleSmFRZ(bMask & ~bSameHemisphere, :));
    dt.logIdxSigClus    = permTest.logIdxSigClus;
    dt.xLim             = thisParam.xLim;
    dt.yLim             = thisParam.yLim;
    dt.yLabel           = thisParam.yLabel;
    f = LK_PlotNeuronalActivityDuringRipples_20231005(dt);
    
    % save figure
    fileName = strcat(r.paths.save, 'zScoredFiringRatesDuringRipples_', thisParam.regionGroups{iRG, 1});
    LK_print(f, fileName, '-dtiff', '-r300');

    % save figure data
    if strcmp(thisParam.regionGroups{iRG}, 'All regions')
        save(strcat(paths.save, 'Fig_4d'), 'dt');
    end
end

%% additional trial characteristics, such as memory performance

% report
fprintf('\nEstimating additional trial characteristics including memory performance:\n');

% random locations in the environment for normalizing drop errors
cfg         = [];
cfg.maxR    = 5000;
cfg.minR    = 0;
cfg.N       = 10000001;
cfg.centerX = 0;
cfg.centerY = 0;
randLocs    = LK_RandomPointsInCircle(cfg);

% loop through neuron-ripple combinations
for iRes = 1:size(r.allRes, 1)
    
    % report
    if mod(iRes, 100) == 1
        fprintf('... for neuron-ripple combination #%d.\n', iRes);
    end
    
    % estimate memory performance
    if iRes == 1 || r.allRes(iRes).idx(1) ~= r.allRes(iRes - 1).idx(1)
        
        % estimate memory performance
        dropError       = sqrt((r.allRes(iRes).trialInfo.xCorrect - r.allRes(iRes).trialInfo.xResponse) .^ 2 + (r.allRes(iRes).trialInfo.yCorrect - r.allRes(iRes).trialInfo.yResponse) .^ 2);
        dropErrorSurro  = pdist2([r.allRes(iRes).trialInfo.xCorrect, r.allRes(iRes).trialInfo.yCorrect], randLocs);
        memPerf         = sum(dropError < dropErrorSurro, 2) ./ sum(~isnan(dropErrorSurro), 2);
    else
        
        % use memory performance from previous cell
        memPerf         = r.allRes(iRes - 1).trialInfo.memPerf;
    end
    
    % assign memory performance to trials
    r.allRes(iRes).trialInfo.memPerf  = memPerf;
    r.allRes(iRes).trialInfo.bGoodMem = memPerf >= median(memPerf, 'omitnan'); % median split across trials
end

%% firing rates during ripples, separately for ripples from different trial phases and different memory performance

% report
fprintf('\nExamining firing rates during ripples from different trial phases.\n');

% parameters for this analysis
thisParam               = [];
thisParam.numSurrogates = param.ana.numSurrogates; % number of surrogates
thisParam.xLim          = [-1, 1]; % time window for statistics and for visualization
thisParam.poi           = {'AllPhases', 'All phases'; 'ITI', 'ITI'; 'Cue', 'Cue'; 'Retrieval', 'Retrieval'; 'Feedback', 'Feedback'; 'Reencoding', 'Re-encoding'}; % trial phases of interest
thisParam.colors        = {rgb('gray'), [0, 0, 0.8]; rgb('gray'), rgb('purple'); rgb('gray'), rgb('red'); rgb('gray'), rgb('orange'); rgb('gray'), rgb('limegreen'); rgb('gray'), rgb('darkgreen')}; % column 1, colors for bad trials
thisParam.FRType        = 'smFRZ'; % type of firing rates
thisParam.yLim          = [-0.05, 0.07];
thisParam.yLabel        = {'Firing rate', '(z-scored)'};
thisParam.base          = 'no'; % whether baseline correction shall be applied

% loop through trial phases
for iPOI = 1:size(thisParam.poi, 1)
    
    % reset rng for reproducibility
    rng(r.param.myRNG);
    
    % report
    fprintf('\n============================================================ Trial phase: %s.\n', thisParam.poi{iPOI, 1});
    
    % preallocate average firing rates of interest
    allMeanFROI         = nan(size(r.allRes, 1), size(r.param.ana.timeCenters, 2)); % average across all available ripples (cells x timepoints)
    allMeanFROIGood     = nan(size(allMeanFROI)); % average across ripples from good trials
    allMeanFROIBad      = nan(size(allMeanFROI)); % average across ripples from bad trials
    
    % loop through neuron-ripple combinations
    for iRes = 1:size(r.allRes, 1)
        
        % trial information for this neuron-ripple combination
        thisTrialInfo   = r.allRes(iRes).trialInfo;
        
        % temporal borders of specific trial phase
        if strcmp(thisParam.poi{iPOI, 1}, 'AllPhases')
            thisTimeEdges   = [thisTrialInfo.ITI, thisTrialInfo.Grab];
        elseif strcmp(thisParam.poi{iPOI, 1}, 'ITI')
            thisTimeEdges   = [thisTrialInfo.ITI, thisTrialInfo.Cue];
        elseif strcmp(thisParam.poi{iPOI, 1}, 'Cue')
            thisTimeEdges   = [thisTrialInfo.Cue, thisTrialInfo.Retrieval];
        elseif strcmp(thisParam.poi{iPOI, 1}, 'Retrieval')
            thisTimeEdges   = [thisTrialInfo.Retrieval, thisTrialInfo.Feedback];
        elseif strcmp(thisParam.poi{iPOI, 1}, 'Feedback')
            thisTimeEdges   = [thisTrialInfo.Feedback, thisTrialInfo.Reencoding];
        elseif strcmp(thisParam.poi{iPOI, 1}, 'Reencoding')
            thisTimeEdges   = [thisTrialInfo.Reencoding, thisTrialInfo.Grab];
        end
        
        % identify ripples belonging to this trial phase
        thisRippleTimepoints    = r.allRes(iRes).rippleTimepoints;
        bROI                    = any(thisRippleTimepoints >= thisTimeEdges(:, 1)' & thisRippleTimepoints <= thisTimeEdges(:, 2)', 2);
        bROIGoodMem             = any(thisRippleTimepoints >= thisTimeEdges(thisTrialInfo.bGoodMem, 1)' & thisRippleTimepoints <= thisTimeEdges(thisTrialInfo.bGoodMem, 2)', 2);
        
        % select specific type of firing rates
        if strcmp(thisParam.FRType, 'smFRZ')
            thisRippleFR            = r.allRes(iRes).rippleSmFRZ; % smoothing was performed before z-scoring
        elseif strcmp(thisParam.FRType, 'FRZ')
            thisRippleFR            = r.allRes(iRes).rippleFRZ;
            thisRippleFR            = smoothdata(thisRippleFR, 2, r.param.ana.smoothType, r.param.ana.smoothFac); % smooth across time
        elseif strcmp(thisParam.FRType, 'FR')
            thisRippleFR            = r.allRes(iRes).rippleFR;
            thisRippleFR            = smoothdata(thisRippleFR, 2, r.param.ana.smoothType, r.param.ana.smoothFac); % smooth across time
        end
        
        % baseline correction if necessary
        if ~strcmp(thisParam.base, 'no')
            bBaseTime               = r.param.ana.timeCenters >= min(thisParam.base) & r.param.ana.timeCenters <= max(thisParam.base);
            baseFR                  = mean(thisRippleFR(:, bBaseTime), 2, 'omitnan');
            thisRippleFR            = thisRippleFR - repmat(baseFR, 1, size(thisRippleFR, 2));
        end
        
        % sanity check
        if size(thisRippleFR, 1) ~= size(bROI, 1)
            error('Problem with data dimensions.');
        end
        
        % average firing rates of interest
        allMeanFROI(iRes, :)        = mean(thisRippleFR(bROI, :), 1, 'omitnan');
        allMeanFROIGood(iRes, :)    = mean(thisRippleFR(bROI & bROIGoodMem, :), 1, 'omitnan');
        allMeanFROIBad(iRes, :)     = mean(thisRippleFR(bROI & ~bROIGoodMem, :), 1, 'omitnan');
    end
    
    %% statistics: firing rates against 0
    
    % test firing rates against zero
    cfg               	= [];
    cfg.mat             = allMeanFROI;
    cfg.alpha        	= 0.05;
    cfg.direction      	= 'pos';
    cfg.numSurrogates 	= thisParam.numSurrogates;
    permTest          	= LK_PermutationTest_OneSampleAgainstZero_20200708(cfg);
    fprintf('\tTime of maximum cluster: %.3f to %.3f s.\n', min(r.param.ana.timeCenters(permTest.logIdxSigClus & permTest.logIdxMaxClus)), max(r.param.ana.timeCenters(permTest.logIdxSigClus & permTest.logIdxMaxClus)));
    
    %% figure: comparison of the firing rates from this condition against 0
    
    % number of combinations from same vs. different hemispheres
    nSame   = sum(bSameHemisphere);
    nDiff   = sum(~bSameHemisphere);
    
    % create figure
    f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1]);
    
    % ripples from hippocampal channels contributing to these combinations
    x   = mean(allRippleERPTime, 1);
    m   = mean(allRippleERPAvg, 1, 'omitnan');
    % subpanel
    ax1 = axes('units', 'centimeters', 'position', [1.75, 4, 3.75, 1]);
    hold on;
    xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
    plot(x, m, 'Color', [0, 0, 0]);
    tx = text(0.575, 0.7, 'HC ripple', 'units', 'normalized');
    tl = title({thisParam.poi{iPOI, 2}, ['(\color[rgb]{0, 0, 0.8}', num2str(nSame), '\color[rgb]{0, 0, 0}/\color[rgb]{0.5, 0.5, 0.5}', num2str(nDiff), '\color[rgb]{0, 0, 0} combinations)']}, ...
        'units', 'normalized', 'position', [0.5, 1.05]);
    set(gca, 'xlim', thisParam.xLim, 'xticklabel', [], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
    set([gca, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    axis off;
    
    % single-neuron activity
    x       = r.param.ana.timeCenters;
    mSame   = mean(allMeanFROI(bSameHemisphere, :), 1, 'omitnan');
    steSame = LK_ste(allMeanFROI(bSameHemisphere, :));
    mDiff   = mean(allMeanFROI(~bSameHemisphere, :), 1, 'omitnan');
    steDiff = LK_ste(allMeanFROI(~bSameHemisphere, :));
    % subpanel
    ax2 = axes('units', 'centimeters', 'position', [1.75, 0.9, 3.75, 3]);
    hold on;
    xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
    yline(0, '--', 'Color', [0.3, 0.3, 0.3]);
    % different hemispheres
    patch([x, fliplr(x)], [mDiff + steDiff, fliplr(mDiff - steDiff)], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, mDiff, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
    % same hemispheres
    patch([x, fliplr(x)], [mSame + steSame, fliplr(mSame - steSame)], [0, 0, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, mSame, 'Color', [0, 0, 0.8], 'LineWidth', 1);
    % significance
    LK_SigLine(x, [max(thisParam.yLim); max(thisParam.yLim) - range(thisParam.yLim) * 0.025], permTest.logIdxSigClus);
    % enhance axes
    xl = xlabel('Time (s)', 'units', 'normalized', 'position', [0.5, -0.09]);
    yl = ylabel(thisParam.yLabel, 'units', 'normalized', 'position', [-0.05, 0.5]);
    set(gca, ...
        'xlim', thisParam.xLim, 'xtick', [min(thisParam.xLim), 0, max(thisParam.xLim)], 'xticklabel', {num2str(min(thisParam.xLim)), '', num2str(max(thisParam.xLim))}, ...
        'ylim', thisParam.yLim, 'ytick', thisParam.yLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
    set([gca, xl, yl], 'FontUnits', 'centimeters', 'fontsize', 0.4);
    drawnow;
    
    % link single-neuron activity and ripples
    linkaxes([ax1, ax2], 'x');
    
    % save figure
    fileName = strcat(r.paths.save, 'allMeanFROI_', thisParam.FRType, '_SameAndDiffHemi_', thisParam.poi{iPOI, 1});
    LK_print(f, fileName, '-dtiff', '-r300');
    
    %% statistics: comparison of firing rates during good-performance ripples vs. bad-performance ripples
    
    % process data for statistics and figure
    time        = r.param.ana.timeCenters;
    n           = {sum(~isnan(allMeanFROIBad), 1), sum(~isnan(allMeanFROIGood), 1)}; % sample sizes: bad, good
    m           = {mean(allMeanFROIBad, 1, 'omitnan'), mean(allMeanFROIGood, 1, 'omitnan')}; % means: bad, good
    ste         = {LK_ste(allMeanFROIBad), LK_ste(allMeanFROIGood)}; % SEMs: bad, good
    
    % re-organize data for fieldtrip
    timelock1   = cell(1, size(allMeanFROIGood, 1));
    timelock2  	= cell(1, size(allMeanFROIBad, 1));
    for iCell = 1:size(allMeanFROIGood, 1)
        
        % good-performance trials
        timelock1{iCell}.avg       = allMeanFROIGood(iCell, :);
        timelock1{iCell}.label     = {'neuron'};
        timelock1{iCell}.time      = time;
        
        % bad-performance trials
        timelock2{iCell}.avg       = allMeanFROIBad(iCell, :);
        timelock2{iCell}.label     = {'neuron'};
        timelock2{iCell}.time      = time;
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = thisParam.xLim;
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'depsamplesT';
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum';
    cfg.neighbours          = [];
    cfg.tail                = 0; % test for both firing-rate increases and decreases
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = thisParam.numSurrogates;
    numCells                = numel(timelock1);
    design                  = zeros(2, numCells * 2);
    design(1, :)            = [1:numCells, 1:numCells]; % cell indices
    design(2, :)            = [ones(1, numCells), ones(1, numCells) * 2]; % 1 = good; 2 = bad
    cfg.design              = design; % design matrix
    cfg.uvar                = 1; % row of the design matrix that contains the units of observation
    cfg.ivar                = 2; % row of the design matrix that contains the independent variable
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    
    % report significant clusters
    LK_report_ft_timelockstatistics(outFT);
    
    %% figure: comparison of firing rates during good-performance ripples vs. bad-performance ripples
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 5.5, 6], 'visible', 'off');
    axes('units', 'centimeters', 'position', [1.5, 1.5, 3.5, 3.5]);
    hold on;
    xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
    yline(0, '--', 'Color', [0.3, 0.3, 0.3]);
    for iM = 1:numel(m) % loop through bad and good trials
        patch([time, fliplr(time)], [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], thisParam.colors{iPOI, iM}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(time, m{iM}, 'Color', thisParam.colors{iPOI, iM});
    end
    % enhance axes
    set(gca, 'xlim', thisParam.xLim, 'ylim', thisParam.yLim, 'ytick', thisParam.yLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    % significance
    tmpAx = get(gca);
    LK_SigLine(outFT.time, [max(tmpAx.YLim); max(tmpAx.YLim) - 0.025 * range(tmpAx.YLim)], outFT.mask, [0, 0, 0]);
    % axis labels and font size
    xl = xlabel('Time (s)');
    yl = ylabel(thisParam.yLabel, 'units', 'normalized', 'position', [-0.05, 0.5]);
    tl = title(['(\color[rgb]{', num2str(thisParam.colors{iPOI, 2}), '}', num2str(min(n{2})), '\color[rgb]{0, 0, 0}, \color[rgb]{', num2str(thisParam.colors{iPOI, 1}), '}', num2str(min(n{1})), '\color[rgb]{0, 0, 0})'], ...
        'units', 'normalized', 'position', [0.5, 1.05]); % report sample size of good data first
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    % save figure
    LK_print(f, strcat(r.paths.save, 'allMeanFROI_', thisParam.FRType, '_GoodVSBad_', thisParam.poi{iPOI, 1}), '-dpng', '-r300');
    close(f);
end
