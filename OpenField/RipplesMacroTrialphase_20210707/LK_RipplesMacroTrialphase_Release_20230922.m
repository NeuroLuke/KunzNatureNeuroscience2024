%==========================================================================
% This script examines the following topics:
% (1) Rate of ripples during different trial phases.
% (2) Relationship between ripple rates, memory, and experience (separately
%     for the different trial phases).
% (3) Development of ripples over time during specific trial phases.
% (4) Influence of IEDs on memory performance and their interaction with
%     ripple rates.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                       = [];
%- subject information
param.info.name             = 'subjectdata_20210616';
%- behavior
param.beh.fileName          = 'trialInfo.mat';
param.beh.trialPhases       = {'ITI'; 'Cue'; 'Retrieval'; 'Feedback'; 'Reencoding'};
param.beh.trialPhasesColor  = [rgb('purple'); rgb('red'); rgb('orange'); rgb('limegreen'); rgb('darkgreen')]; % colors for the different phases
param.beh.PSTHTimeRes       = 0.001; % time resolution for ripple PSTH
param.beh.PSTHStartEdges    = -2:param.beh.PSTHTimeRes:4; % timeline for start-locked ripple PSTHs
param.beh.PSTHStartTime     = movmean(param.beh.PSTHStartEdges, 2, 'endpoints', 'discard');
param.beh.PSTHEndEdges      = -4:param.beh.PSTHTimeRes:2; % timeline for end-locked ripple PSTHs
param.beh.PSTHEndTime       = movmean(param.beh.PSTHEndEdges, 2, 'endpoints', 'discard');
param.beh.maxNumTrials      = 160; % maximum number of trials
%- channels
param.channel.choi          = {'HCleftBipolar'; 'HCrightBipolar'};
%- ripples
param.ripple.folder         = 'RipplesMacro_20210614';
param.ripple.specs          = '20220320_FBS_80to140Hz_tMF2';
param.ripple.type           = 'ripples';
param.ripple.fileName       = 'ripples.mat';
param.ripple.dataFileName   = 'data4Ripples.mat';
%- IEDs
param.IED.type              = 'artifacts';
param.IED.fileName          = strcat(param.IED.type, '.mat');
%- analysis
param.ana.bRankMemPerf      = true; % whether to rank the memory-performance values
param.ana.bCenterPredictors = true; % whether to center the predictors
param.ana.corrType          = 'Pearson'; % type of correlation
param.ana.numSurrogates     = 1001; % number of surrogates

% set rng for reproducibility
param.myRNG                 = 7777;
rng(param.myRNG);

% paths
paths       	= [];
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.beh    	= 'E:\OpenField\BehPreprocessing_20210707\';
paths.ripple 	= strcat('E:\OpenField\', param.ripple.folder, '\', param.ripple.specs, '\');
paths.save    	= strcat('E:\OpenField\RipplesMacroTrialphase_20210707\20220625_', param.ripple.folder, '_', param.ripple.specs, '_', param.ripple.type, filesep);
mkdir(paths.save);

% add fieldtrip
addpath(paths.fieldtrip);
ft_defaults;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% save settings
save(strcat(paths.save, 'settings'));

%% surrogate locations to convert drop errors into memory performances

% random locations in the environment
cfg         = [];
cfg.maxR    = 5000;
cfg.minR    = 0;
cfg.N       = 10000001;
cfg.centerX = 0;
cfg.centerY = 0;
randLocs    = LK_RandomPointsInCircle(cfg);

%% preallocations

% main results for each channel
allRes      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)

    %% subject information
    
    % report
    fprintf('\n============================================================ Subject: %s.\n', subjects{iSub});
    
    % load subject information
    subjectdata = load(fullfile(paths.subjects, subjects{iSub}, param.info.name));
    subjectdata = subjectdata.subjectdata;
    
    %% trial information
    
    % trial information in MACROTIME
    trialInfo       = load(fullfile(paths.beh, subjects{iSub}, param.beh.fileName));
    trialInfo       = trialInfo.trialInfoMacrotime;
    fprintf('Original number of trials: %d.\n', size(trialInfo, 1));
    
    % restrict to 160 trials
    if size(trialInfo, 1) > param.beh.maxNumTrials
        trialInfo   = trialInfo(1:param.beh.maxNumTrials, :);
    end
    fprintf('Corrected number of trials: %d.\n', size(trialInfo, 1));
    
    %% channels of interest
    
    % channels of interest from this subject
    choi        = [];
    if any(strcmp(param.channel.choi, 'HCleftBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCleftBipolar);
    end
    if any(strcmp(param.channel.choi, 'HCrightBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCrightBipolar);
    end
    
    %% loop through channels of interest
    for iChan = 1:numel(choi)
        
        %% ripple and IED data
        
        % load ripple data
        ripples = load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, param.ripple.fileName));
        ripples = ripples.ripples;
        
        % load IED data
        IEDs    = load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, param.IED.fileName));
        IEDs    = IEDs.artifacts;
        
        % report
        fprintf('\nChannel: %s. Total number of ripples: %d.\n', choi{iChan}, size(ripples.ripples, 1));
        
        %% ripple, IED, and artifact data
        
        % load ripple data
        data4Ripples   	= load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, param.ripple.dataFileName));
        data4Ripples 	= data4Ripples.data4Ripples;
                
        % time periods with artifacts
        bArtifact       = isnan(data4Ripples.dataEB{1});
        
        % identify time periods with IEDs
        bIED            = IEDs.bArtifact;
        
        %% ripple rate and fraction of artifacts overall
        
        % time of ripple peaks
        rippleTimepoints  	= transpose([ripples.ripples.peakTime]);
        
        % total number of ripples
        totNumRipples       = size(ripples.ripples, 1);
        totExpDur           = numel(ripples.bRipple) / ripples.eegRipples.fsample;
        meanRippleRate      = totNumRipples / totExpDur;
        
        % duration of all single ripples
        durationPerRipple   = transpose([ripples.ripples.endTime] - [ripples.ripples.startTime]);
               
        % time fraction of artifacts
        totFracArtifacts    = sum(bArtifact) / numel(bArtifact);
        
        % time fraction of IEDs
        totFracIEDs         = sum(bIED) / numel(bIED);
        
        % report
        fprintf('Channel: %s. Mean ripple rate: %.3f (Hz).\n', choi{iChan}, meanRippleRate);
        fprintf('Channel: %s. Total fraction of artifacts: %.3f (based on "dataEB").\n', choi{iChan}, totFracArtifacts);
        
        %% ripple rate in the first vs. second half of all trials
        
        % temporal borders of first and second half
        numTrials   = size(trialInfo, 1);
        half1Time   = [min(trialInfo{1, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}}), ...
            max(trialInfo{floor(numTrials / 2), {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}})];
        half2Time   = [min(trialInfo{floor(numTrials / 2 + 1), {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}}), ...
            max(trialInfo{end, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}})];
        
        % ripple rate in first and second trial half
        meanRippleRateHalves    = [sum(rippleTimepoints >= min(half1Time) & rippleTimepoints <= max(half1Time)) / range(half1Time), ...
            sum(rippleTimepoints >= min(half2Time) & rippleTimepoints <= max(half2Time)) / range(half2Time)];
        
        %% ripple rate and fraction of artifacts in different trial phases
        
        % characteristics of ripples per trial and trial phase
        numRipplesPerPhase      = nan(size(trialInfo, 1), numel(param.beh.trialPhases)); % ripple number per trial phase (e.g., size = 100 trials x 5 phases)
        durationPerPhase        = nan(size(trialInfo, 1), numel(param.beh.trialPhases)); % duration per trial phase
        rippleDurPerPhase       = nan(size(trialInfo, 1), numel(param.beh.trialPhases)); % ripple duration per trial phase
        rippleFreqPerPhase      = nan(size(trialInfo, 1), numel(param.beh.trialPhases)); % ripple frequency per trial phase
        
        % peristimulus time histogram of ripples
        ripplePSTHAbsStart      = nan(size(trialInfo, 1), numel(param.beh.trialPhases), numel(param.beh.PSTHStartTime)); % ripple PSTH, start-locked
        ripplePSTHAbsEnd        = nan(size(trialInfo, 1), numel(param.beh.trialPhases), numel(param.beh.PSTHEndTime)); % ripple PSTH, end-locked
        
        % fraction of IEDs and artifacts per trial phase
        fracIEDsPerPhase        = nan(size(trialInfo, 1), numel(param.beh.trialPhases)); % e.g., size = 100 trials x 5 phases
        fracArtifactsPerPhase   = nan(size(trialInfo, 1), numel(param.beh.trialPhases));
        
        % loop through trials and trial phases
        for iTrial = 1:size(trialInfo, 1)
            for iPhase = 1:numel(param.beh.trialPhases)
                
                % temporal borders of this trial phase
                thisPhaseName   = param.beh.trialPhases{iPhase};
                thisBorders     = [nan, nan]; % reset borders
                switch thisPhaseName
                    case 'ITI'
                        thisBorders = [trialInfo.ITI(iTrial), trialInfo.Cue(iTrial)];
                    case 'Cue'
                        thisBorders = [trialInfo.Cue(iTrial), trialInfo.Retrieval(iTrial)];
                    case 'Retrieval'
                        thisBorders = [trialInfo.Retrieval(iTrial), trialInfo.Feedback(iTrial)];                        
                    case 'Feedback'
                        thisBorders = [trialInfo.Feedback(iTrial), trialInfo.Reencoding(iTrial)];
                    case 'Reencoding'
                        thisBorders = [trialInfo.Reencoding(iTrial), trialInfo.Grab(iTrial)];
                end
                
                % continue if there is a nan in the borders (can happen for
                % last trials, for example)
                if any(isnan(thisBorders))
                    fprintf('\tSkipping trial #%d, phase "%s".\n', iTrial, thisPhaseName);
                    continue;
                end
                
                %% IEDs and artifacts in this trial phase
                
                % time within this trial phase
                bTimeWithinBorders                      = data4Ripples.time{1} >= min(thisBorders) & data4Ripples.time{1} <= max(thisBorders);
                
                % percentage of IEDs in this trial phase
                fracIEDsPerPhase(iTrial, iPhase)        = sum(bIED(bTimeWithinBorders)) / sum(bTimeWithinBorders);
                
                % percentage of artifacts in this trial phase
                fracArtifactsPerPhase(iTrial, iPhase)   = sum(bArtifact(bTimeWithinBorders)) / sum(bTimeWithinBorders);                
                
                %% ripple characteristics in this trial phase
                
                % duration of this trial phase
                durationPerPhase(iTrial, iPhase)        = range(thisBorders);
                
                % identify the ripples occurring in this trial phase
                bRippleWithinBorders                    = rippleTimepoints >= min(thisBorders) & rippleTimepoints < max(thisBorders);
                
                % number of ripples in this trial phase
                numRipplesPerPhase(iTrial, iPhase)      = sum(bRippleWithinBorders);
                
                % average ripple duration in this trial phase
                rippleDurPerPhase(iTrial, iPhase)       = mean(cell2mat({ripples.ripples(bRippleWithinBorders).duration}'));
                
                % average ripple frequency in this trial phase
                rippleFreqPerPhase(iTrial, iPhase)      = mean(cell2mat({ripples.ripples(bRippleWithinBorders).freq}'));
                
                %% ripple PSTH, relative to start
                
                % time bins relative to the start of the phase
                thisBorders4PSTH                        = min(thisBorders) + param.beh.PSTHStartEdges;
                
                % ripple PSTH relative to start
                ripplePSTHAbsStart(iTrial, iPhase, :)   = histcounts(rippleTimepoints, thisBorders4PSTH);
                
                %% ripple PSTH, relative to end
                
                % time bins relative to the end of the phase
                thisBorders4PSTH                        = max(thisBorders) + param.beh.PSTHEndEdges;
                
                % ripple PSTH relative to end
                ripplePSTHAbsEnd(iTrial, iPhase, :)     = histcounts(rippleTimepoints, thisBorders4PSTH);
                
            end
        end
        
        % ripple rate
        rippleRatePerPhase         	= numRipplesPerPhase ./ durationPerPhase; % (Hz)
        
        % ripple rate during artifact-free periods (note: contains nans if
        % a period is fully covered by artifacts)
        rippleRatePerPhaseArtFree   = numRipplesPerPhase ./ (durationPerPhase .* (1 - fracArtifactsPerPhase)); % (Hz)
        
        %% fraction of artifacts and IEDs per trial
        
        % preallocate
        fracArtifactsPerTrial   = nan(size(trialInfo, 1), 1); % e.g., size = 100 trials x 1
        fracIEDsPerTrial        = nan(size(trialInfo, 1), 1);
        
        % loop through trials
        for iTrial = 1:size(trialInfo, 1)
            
            % temporal borders of this trial
            thisBorders     = [trialInfo.ITI(iTrial), trialInfo.Grab(iTrial)];
            
            % time within this trial phase
            bTimeWithinBorders                  = data4Ripples.time{1} >= min(thisBorders) & data4Ripples.time{1} <= max(thisBorders);
            
            % fraction of artifacts in this trial
            fracArtifactsPerTrial(iTrial, 1)    = sum(bArtifact(bTimeWithinBorders)) / sum(bTimeWithinBorders);
            
            % fraction of IEDs in this trial
            fracIEDsPerTrial(iTrial, 1)         = sum(bIED(bTimeWithinBorders)) / sum(bTimeWithinBorders);
        end
        
        %% results
        
        % results for this channel
        thisRes                             = [];
        thisRes.subject                     = subjects{iSub};
        thisRes.chan                        = choi{iChan};
        thisRes.idx                         = [iSub, iChan];
        
        % behavior
        thisRes.trialInfo                   = trialInfo;
        thisRes.durationPerPhase            = durationPerPhase; % trials x phases
        thisRes.totExpDur                   = totExpDur;
        
        % overall ripple characteristics
        thisRes.rippleTimepoints            = rippleTimepoints; % timepoints of the ripple peaks
        thisRes.totNumRipples               = totNumRipples;
        thisRes.meanRippleRate              = meanRippleRate;
        thisRes.durationPerRipple           = durationPerRipple;
        
        % ripple rate per experiment half
        thisRes.meanRippleRateHalves        = meanRippleRateHalves;
        
        % ripple characteristics (per trial and per trial phase)
        thisRes.numRipplesPerPhase          = numRipplesPerPhase; % trials x phases
        thisRes.rippleRatePerPhase          = rippleRatePerPhase;
        thisRes.rippleRatePerPhaseArtFree   = rippleRatePerPhaseArtFree;
        thisRes.rippleDurPerPhase           = rippleDurPerPhase;
        thisRes.rippleFreqPerPhase          = rippleFreqPerPhase;
        
        % ripple PSTHs
        thisRes.ripplePSTHAbsStart          = ripplePSTHAbsStart;
        thisRes.ripplePSTHAbsEnd            = ripplePSTHAbsEnd;
        
        % artifact characteristics (overall and per trial/phase)
        thisRes.totFracArtifacts            = totFracArtifacts; % total fraction of artifacts across the entire experiment
        thisRes.fracArtifactsPerTrial       = fracArtifactsPerTrial;
        thisRes.fracArtifactsPerPhase       = fracArtifactsPerPhase;
        
        % IED characteristics (overall and per trial/phase)
        thisRes.totFracIEDs                 = totFracIEDs;
        thisRes.fracIEDsPerTrial            = fracIEDsPerTrial;
        thisRes.fracIEDsPerPhase            = fracIEDsPerPhase;
        
        % collect results across all sessions
        allRes                              = cat(1, allRes, thisRes);
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r   = load(strcat(paths.save, 'results.mat'));

% general information
chanIdx                 = cell2mat({r.allRes.idx}');
chanSubjects            = {r.allRes.subject}';
chanBMicroSubject       = ismember(chanSubjects, r.subjects(strcmp(r.subjects(:, 2), 'Microwire'), 1));
chanBRightHemisphere    = contains({r.allRes.chan}', 'R');
chanBFirstSession       = cellfun(@(x) ~contains(x(end), {'b', 'c'}), chanSubjects); % first sessions
chanBSecondSession      = cellfun(@(x) contains(x(end), 'b'), chanSubjects);
chanBThirdSession       = cellfun(@(x) contains(x(end), 'c'), chanSubjects);

% more specific information
chanNumTrials           = cellfun(@(x) size(x, 1), {r.allRes.ripplePSTHAbsStart}');
totalNumTrials          = sum(chanNumTrials);

% report
fprintf('\n\n============================================================== Results.\n');
fprintf('Number of channels: %d.\n', size(r.allRes, 1));
fprintf('Number of channels with microwires: %d.\n', sum(chanBMicroSubject));
fprintf('Number of channels from the right hemisphere: %d.\n', sum(chanBRightHemisphere));
fprintf('Number of sessions: %d.\n', size(unique(chanSubjects), 1));
fprintf('Number of second sessions: %d.\n', sum(chanBSecondSession));
fprintf('Number of third sessions: %d.\n', sum(chanBThirdSession));
fprintf('Total number of trials: %d.\n', totalNumTrials);

%% mean ripple rate per channel

% overall average ripple rate per channel
chanMeanRippleRate          = cell2mat({r.allRes.meanRippleRate}');
LK_ReportMeanAndSEM_20220322('Ripple rate', chanMeanRippleRate);

% overall average ripple rate per channel, only considering artifact-free
% data
chanMeanRippleRateArtFree   = [r.allRes.totNumRipples]' ./ ([r.allRes.totExpDur]' .* (1 - [r.allRes.totFracArtifacts]'));
LK_ReportMeanAndSEM_20220322('Ripple rate, artifact-free periods only', chanMeanRippleRateArtFree);

%% effect of trial phase on ripple rate, ripple duration, and ripple frequency (repeated measures ANOVA)

% data groups
groups  = {'chanRippleRatePerPhase'; 'chanRippleRatePerPhaseArtFree'; ...
    'chanRippleDurPerPhase'; 'chanRippleFreqPerPhase'; ...
    'chanFracArtifactsPerPhase'; 'chanFracIEDsPerPhase'};

% average ripple rate, duration, and frequency per phase and per channel
chanRippleRatePerPhase          = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.rippleRatePerPhase}', 'UniformOutput', false));
chanRippleRatePerPhaseArtFree   = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.rippleRatePerPhaseArtFree}', 'UniformOutput', false));
chanRippleDurPerPhase           = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.rippleDurPerPhase}', 'UniformOutput', false));
chanRippleFreqPerPhase          = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.rippleFreqPerPhase}', 'UniformOutput', false));

% average percentage-of-artifacts per phase and per channel
chanFracArtifactsPerPhase       = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.fracArtifactsPerPhase}', 'UniformOutput', false));
chanFracIEDsPerPhase            = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), {r.allRes.fracIEDsPerPhase}', 'UniformOutput', false));

% loop through the different data groups
for iGroup = 1:numel(groups)
    
    %% data
    
    % select data
    if strcmp(groups{iGroup}, 'chanRippleRatePerPhase')
        myData      = chanRippleRatePerPhase;
        myXLabel    = 'Rate (Hz)';
        myTitle     = {'Ripple rate'};
        myPrec      = 3;
    elseif strcmp(groups{iGroup}, 'chanRippleRatePerPhaseArtFree')
        myData      = chanRippleRatePerPhaseArtFree;
        myXLabel    = 'Rate (Hz)';
        myTitle     = {'Ripple rate (adj)'};
        myPrec      = 3;
    elseif strcmp(groups{iGroup}, 'chanRippleDurPerPhase')
        myData      = chanRippleDurPerPhase;
        myXLabel    = 'Duration (s)';
        myTitle     = {'Ripple duration'};
        myPrec      = 3;
    elseif strcmp(groups{iGroup}, 'chanRippleFreqPerPhase')
        myData      = chanRippleFreqPerPhase;
        myXLabel    = 'Frequency (Hz)';
        myTitle     = {'Ripple frequency'};
        myPrec      = 1;
    elseif strcmp(groups{iGroup}, 'chanFracArtifactsPerPhase')
        myData    = chanFracArtifactsPerPhase;
        myXLabel    = 'Fraction of time';
        myTitle     = {'Artifact prevalence'};
        myPrec      = 2;
    elseif strcmp(groups{iGroup}, 'chanFracIEDsPerPhase')
        myData    = chanFracIEDsPerPhase;
        myXLabel    = 'Fraction of time';
        myTitle     = {'IED prevalence'};
        myPrec      = 2;
    end
    
    %% stats
    
    % convert data into a table and rename trial phases
    withinNames     = {'p1', 'p2', 'p3', 'p4', 'p5'}; % phases 1-5
    data4ANOVA      = array2table(myData, 'VariableNames', withinNames);
    
    % create within design
    withinDesign    = table(categorical(transpose(1:numel(withinNames))), 'VariableNames', {'Phase'});
    
    % fit repeated measures model and calculate repeated measures ANOVA
    rm              = fitrm(data4ANOVA, 'p1-p5~1', 'WithinDesign', withinDesign);
    ranovaTbl       = ranova(rm);
    
    % post-hoc tests controlling for multiple comparisons
    ranovaMC        = multcompare(rm, 'Phase', 'ComparisonType', 'tukey-kramer'); % returns multiple comparisons of the estimated marginal means based on the repeated measures model rm
    
    % report
    fprintf('\nRepeated measures ANOVA testing for the effect of trial phase on "%s":\n', groups{iGroup});
    disp(sum(~isnan(myData))); % sample size
    disp(ranovaTbl);
    disp(ranovaMC);
    
    %% figure
    
    % create figure
    dt          = [];
    dt.beh      = param.beh;
    dt.myData   = myData;
    dt.myPrec   = myPrec;
    dt.group    = groups{iGroup};
    dt.myXLabel = myXLabel;
    dt.myTitle  = myTitle;
    f = LK_PlotRipplePropertyPerPhase_20231030(dt);

    % save figure
    LK_print(f, strcat(paths.save, groups{iGroup}, '_20230922'), '-dtiff', '-r300');

    % save figure data
    if strcmp(groups{iGroup}, 'chanRippleRatePerPhase')
        save(strcat(r.paths.save, 'Fig_3a_left'), 'dt');
    elseif strcmp(groups{iGroup}, 'chanRippleDurPerPhase')
        save(strcat(r.paths.save, 'Fig_3a_middle'), 'dt');
    elseif strcmp(groups{iGroup}, 'chanRippleFreqPerPhase')
        save(strcat(r.paths.save, 'Fig_3a_right'), 'dt');
    end
end

%% get information about drop errors and memory performance

% report
fprintf('\nExtracting information about drop errors and memory performance.\n');

% trial-wise memory performance for later use
chanDropError       = cell(size(r.allRes, 1), 1); % drop errors per channel
chanMemPerf         = cell(size(r.allRes, 1), 1);
chanBGoodMem        = cell(size(r.allRes, 1), 1);
chanTrialIdx        = cell(size(r.allRes, 1), 1);
chanTrialIdxMemPerf = cell(size(r.allRes, 1), 1);

% loop through channels
for iChan = 1:size(r.allRes, 1)
    
    % report
    fprintf('- %s.\n', chanSubjects{iChan});
        
    % trial info for this channel
    trialInfo   = r.allRes(iChan).trialInfo;
    
    % trial index
    trialIdx    = trialInfo.TrialIdx;
    
    % trial-wise drop error and conversion into memory performance
    dropError                   = sqrt((trialInfo.xCorrect - trialInfo.xResponse) .^ 2 + (trialInfo.yCorrect - trialInfo.yResponse) .^ 2);
    dropErrorSurro              = pdist2([trialInfo.xCorrect, trialInfo.yCorrect], randLocs); % surrogate drop errors
    memPerf                     = sum(dropError < dropErrorSurro, 2) ./ sum(~isnan(dropErrorSurro), 2); % convert drop error into memory performance
    memPerf(isnan(dropError))   = nan;
    
    % rank memory performance
    if param.ana.bRankMemPerf == true
        [~, ~, memPerf]   	    = unique(memPerf);
    end
    
    % center the continuous predictor variables
    if param.ana.bCenterPredictors == true
        memPerf                 = memPerf - mean(memPerf, 'omitnan');
        trialIdx          	    = trialIdx - mean(trialIdx, 'omitnan');
    end
    
    % interaction term
    trialIdxMemPerf      	    = memPerf .* trialIdx;
    
    % collect trial-wise information across all channels
    chanDropError{iChan}        = dropError;
    chanMemPerf{iChan}          = memPerf;
    chanBGoodMem{iChan}         = memPerf > median(memPerf, 'omitnan'); % median split
    chanTrialIdx{iChan}         = trialIdx;
    chanTrialIdxMemPerf{iChan}  = trialIdxMemPerf;
end

%% correlation between ripple rates, memory performance, and trial index using partial correlations

% data groups
dataGroups  = {'rippleRatePerPhase', [-0.5, 0.5], {'Ripple rate', ''}; ...
    'rippleRatePerPhaseArtFree', [-0.5, 0.5], {'Ripple rate', '(adj)'}; ...
    'rippleDurPerPhase', [-1, 1], {'Ripple duration', ''}; ...
    'rippleFreqPerPhase', [-1, 1], {'Ripple frequency', ''}; ...
    'fracArtifactsPerPhase', [-0.6, 0.6], {'Artifact prevalence', ''}; ...
    'fracIEDsPerPhase', [-0.6, 0.6], {'IED prevalence', ''}};

% correlation groups
corrGroups  = {'chanCorr2Memory'; 'chanCorr2Trial'; 'chanCorr2TrialMemory'};

% loop through the different data groups
for iDG = 1:size(dataGroups, 1)
    
    % report
    fprintf('\n==== Analysis of the correlation between neural data, memory performance, and trial index. Neural data: %s\n', dataGroups{iDG, 1});
    
    % preallocate correlations to ripple rates (channels x trial phases)
    chanCorr2Trial          = nan(size(r.allRes, 1), numel(param.beh.trialPhases)); % correlation with trial index
    chanCorr2Memory         = nan(size(r.allRes, 1), numel(param.beh.trialPhases)); % correlation with memory performance
    chanCorr2TrialMemory    = nan(size(r.allRes, 1), numel(param.beh.trialPhases)); % correlation with trial index * memory performance
    
    % loop through channels
    for iChan = 1:size(r.allRes, 1)
        
        % select the neural data for computing the correlations
        if strcmp(dataGroups{iDG, 1}, 'rippleRatePerPhase')
            myData  = r.allRes(iChan).rippleRatePerPhase;
        elseif strcmp(dataGroups{iDG, 1}, 'rippleRatePerPhaseArtFree')
            myData  = r.allRes(iChan).rippleRatePerPhaseArtFree;
        elseif strcmp(dataGroups{iDG, 1}, 'rippleDurPerPhase')
            myData  = r.allRes(iChan).rippleDurPerPhase;
        elseif strcmp(dataGroups{iDG, 1}, 'rippleFreqPerPhase')
            myData  = r.allRes(iChan).rippleFreqPerPhase;
        elseif strcmp(dataGroups{iDG, 1}, 'fracArtifactsPerPhase')
            myData  = r.allRes(iChan).fracArtifactsPerPhase;
        elseif strcmp(dataGroups{iDG, 1}, 'fracIEDsPerPhase')
            myData  = r.allRes(iChan).fracIEDsPerPhase;
        end

        % phase-specific correlation of neural data with memory
        % performance, controlling for trial index and the interaction
        mainPred    = chanMemPerf{iChan};
        otherPred 	= [chanTrialIdx{iChan}, chanTrialIdxMemPerf{iChan}];
        for iPhase = 1:numel(param.beh.trialPhases)
            try
                chanCorr2Memory(iChan, iPhase)  = partialcorr(mainPred, myData(:, iPhase), otherPred, 'rows', 'complete', 'type', param.ana.corrType);
            catch
                fprintf('\tSkipping: %s, channel %d, %s.\n', dataGroups{iDG, 1}, iChan, r.param.beh.trialPhases{iPhase});
            end
        end
        
        % phase-specific correlation of neural data with trial index,
        % controlling for memory performance and the interaction
        mainPred    = chanTrialIdx{iChan};
        otherPred   = [chanMemPerf{iChan}, chanTrialIdxMemPerf{iChan}];
        for iPhase = 1:numel(param.beh.trialPhases)
            try
                chanCorr2Trial(iChan, iPhase)   = partialcorr(mainPred, myData(:, iPhase), otherPred, 'rows', 'complete', 'type', param.ana.corrType);
            catch
                fprintf('\tSkipping: %s, channel %d, %s.\n', dataGroups{iDG, 1}, iChan, r.param.beh.trialPhases{iPhase});
            end
        end
        
        % phase-specific correlation of neural data with trial index *
        % memory performance, controlling for trial index and memory
        % performance
        mainPred        = chanTrialIdxMemPerf{iChan};
        otherPred       = [chanTrialIdx{iChan}, chanMemPerf{iChan}];
        for iPhase = 1:numel(param.beh.trialPhases)
            try
                chanCorr2TrialMemory(iChan, iPhase) = partialcorr(mainPred, myData(:, iPhase), otherPred, 'rows', 'complete', 'type', param.ana.corrType);
            catch
                fprintf('\tSkipping: %s, channel %d, %s.\n', dataGroups{iDG, 1}, iChan, r.param.beh.trialPhases{iPhase});
            end
        end
    end
    
    %% evaluate the correlations across channels
    
    % loop through the different correlation groups
    for iCG = 1:size(corrGroups, 1)
        
        %% data
        
        % select specific data
        if strcmp(corrGroups{iCG}, 'chanCorr2Memory')
            myData    = chanCorr2Memory;
        elseif strcmp(corrGroups{iCG}, 'chanCorr2Trial')
            myData    = chanCorr2Trial;
        elseif strcmp(corrGroups{iCG}, 'chanCorr2TrialMemory')
            myData    = chanCorr2TrialMemory;
        end
        
        %% stats: test correlation values against 0
        
        % t-test of correlation values against zero, separately for each
        % trial phase
        [~, p_corr2Beh, ~, stats_corr2Beh]  = ttest(myData, 0);
        
        % Bonferroni correction of the p-values
        p_corr2Beh  = p_corr2Beh .* numel(p_corr2Beh);
        
        % report
        fprintf('\nT-test on the correlation values from "%s" (p-values are Bonferroni corrected):\n', corrGroups{iCG});
        disp(sum(~isnan(myData)));
        disp(stats_corr2Beh);
        disp(p_corr2Beh);
        
        %% figure
        
        % create figure
        dt              = [];
        dt.myData       = myData;
        dt.beh          = param.beh;
        dt.p_corr2Beh   = p_corr2Beh;
        dt.dataGroup    = dataGroups(iDG, :);
        dt.corrGroup    = corrGroups{iCG};
        f = LK_PlotRippleCorrelationsWithBeh_20231030(dt);

        % save figure
        LK_print(f, strcat(paths.save, corrGroups{iCG}, '_', dataGroups{iDG, 1}, '_20230923'), '-dtiff', '-r300');

        % save figure data
        if strcmp(dataGroups{iDG}, 'rippleRatePerPhase') && strcmp(corrGroups{iCG}, 'chanCorr2Memory')
            save(strcat(r.paths.save, 'Fig_3b'), 'dt');
        elseif strcmp(dataGroups{iDG}, 'rippleRatePerPhase') && strcmp(corrGroups{iCG}, 'chanCorr2Trial')
            save(strcat(r.paths.save, 'Fig_3d'), 'dt');
        end
    end
end

%% linear mixed model for the effects of memory performance, trial phase, and trial index on ripple rate

% report
fprintf('\nLinear mixed model to test the influence of various predictors on ripple rates.\n');

% preallocate response and predictor variables
RR4LME          = []; % ripple rate in each channel*trial*phase
RRArtFree4LME   = []; % ripple rate during artifact-free periods in each channel*trial*phase
AF4LME          = []; % artifact fraction in each channel*trial*phase
IEDFrac4LME     = []; % IED fraction in each channel*trial*phase
memPerf4LME     = []; % memory performance in each channel*trial*phase
trialIdx4LME    = []; % trial index in each channel*trial*phase
trialPhase4LME  = []; % trial phase in each channel*trial*phase
chanIdx4LME     = []; % channel index  in each channel*trial*phase
subjIdx4LME     = []; % subject index in each channel*trial*phase
bFirstSess4LME  = []; % whether from a first session
for iChan = 1:size(r.allRes, 1)
    
    %% get data
    
    % each ripple rate, each artifact fraction, each trial phase, and each
    % trial index for this channel
    thisRR          = r.allRes(iChan).rippleRatePerPhase; % trials x phase; e.g., size = 100 x 5
    thisRRArtFree   = r.allRes(iChan).rippleRatePerPhaseArtFree;
    thisAF          = r.allRes(iChan).fracArtifactsPerPhase;
    thisIEDFrac     = r.allRes(iChan).fracIEDsPerPhase;
    thisTrialPhase  = ones(size(thisRR, 1), numel(param.beh.trialPhases)) .* (1:numel(param.beh.trialPhases));
    thisTrialIdx    = repmat(r.allRes(iChan).trialInfo.TrialIdx, 1, size(thisRR, 2));
    
    % each memory performance for this channel
    thisMemPerf   	= repmat(chanMemPerf{iChan}, 1, size(thisRR, 2)); % same for all trial phases
    
    %% reshape data
    
    % reshape data to organize it as phase 1:p, trial 1:t
    thisRR          = reshape(transpose(thisRR), numel(thisRR), 1);
    thisRRArtFree   = reshape(transpose(thisRRArtFree), numel(thisRRArtFree), 1);
    thisAF          = reshape(transpose(thisAF), numel(thisAF), 1);
    thisIEDFrac     = reshape(transpose(thisIEDFrac), numel(thisIEDFrac), 1);
    thisTrialPhase  = reshape(transpose(thisTrialPhase), numel(thisTrialPhase), 1);
    thisTrialIdx    = reshape(transpose(thisTrialIdx), numel(thisTrialIdx), 1);
    thisMemPerf     = reshape(transpose(thisMemPerf), numel(thisMemPerf), 1);
        
    % each channel index and each subject index for all trials and phases
    thisChanIdx     = repmat(iChan, size(thisRR, 1), 1);
    thisSubjIdx     = repmat(r.allRes(iChan).idx(1), size(thisRR, 1), 1);
    
    % whether from a first session
    thisBFirstSess  = repmat(chanBFirstSession(iChan), size(thisRR, 1), 1);
    
    %% center continuous predictor variables
    
    % center continuous predictor variables
    thisTrialIdx    = thisTrialIdx - mean(thisTrialIdx, 'omitnan');
    thisMemPerf     = thisMemPerf - mean(thisMemPerf, 'omitnan');
    
    %% concatenate for LME
    
    % concatenate all data
    RR4LME       	= cat(1, RR4LME, thisRR);
    RRArtFree4LME   = cat(1, RRArtFree4LME, thisRRArtFree);
    AF4LME          = cat(1, AF4LME, thisAF);
    IEDFrac4LME     = cat(1, IEDFrac4LME, thisIEDFrac);
    memPerf4LME  	= cat(1, memPerf4LME, thisMemPerf);
    trialIdx4LME  	= cat(1, trialIdx4LME, thisTrialIdx);
    trialPhase4LME  = cat(1, trialPhase4LME, thisTrialPhase);
    chanIdx4LME   	= cat(1, chanIdx4LME, thisChanIdx);
    subjIdx4LME    	= cat(1, subjIdx4LME, thisSubjIdx);
    bFirstSess4LME  = cat(1, bFirstSess4LME, thisBFirstSess);
end

% data table for LME
data4LME    = table(RR4LME, RRArtFree4LME, AF4LME, IEDFrac4LME, memPerf4LME, trialIdx4LME, categorical(trialPhase4LME), categorical(chanIdx4LME), categorical(subjIdx4LME), categorical(bFirstSess4LME), ...
    'VariableNames', {'RR', 'RRArtFree', 'AF', 'IEDFrac', 'MemPerf', 'TrialIdx', 'TrialPhase', 'ChanIdx', 'SubjIdx', 'bFirstSess'});

%% linear mixed model for the effects of memory performance, trial phase, and trial index on ripple rate

% fit linear mixed model to analyze ripple rates
fprintf('\n==== LME investigating the modulation of ripple rates.\n');
LMERippleRate       = fitlme(data4LME, 'RR ~ 1 + MemPerf * TrialPhase * TrialIdx + (1 | ChanIdx)');
statsLMERippleRate  = anova(LMERippleRate);
disp(LMERippleRate);
disp(statsLMERippleRate);

%% linear mixed model for the effects of memory performance, trial phase, and trial index on ripple rate - controlling for artifact fraction

% fit linear mixed model to analyze ripple rates
fprintf('\n==== LME investigating the modulation of ripple rates, controlling for artifacts.\n');
LMERippleRate       = fitlme(data4LME, 'RR ~ 1 + AF + MemPerf * TrialPhase * TrialIdx + (1 | ChanIdx)');
statsLMERippleRate  = anova(LMERippleRate);
disp(LMERippleRate);
disp(statsLMERippleRate);

%% linear mixed model for the effects of memory performance, trial phase, and trial index on IED fraction

% fit linear mixed model to analyze IED fractions
fprintf('\n==== LME investigating the modulation of IED fractions.\n');
LMEArtifacts        = fitlme(data4LME, 'IEDFrac ~ 1 + MemPerf * TrialPhase * TrialIdx + (1 | ChanIdx)');
statsLMEArtifacts   = anova(LMEArtifacts);
disp(LMEArtifacts);
disp(statsLMEArtifacts);

%% ripple rates during the first vs. second half of all trials, separately for first and later sessions

% mean ripple rate in the first vs. second half of all trials
chanMeanRippleRateHalves    = cell2mat({r.allRes.meanRippleRateHalves}');

% different session groups
sessGroups  = {...
    'allSess', true(size(chanBFirstSession)), false, 1, 'All sessions'; ... % session label, mask, Bonferroni correction, strength of Bonferroni correction, title
    'firstSess', chanBFirstSession, true, 2, 'First sessions'; ...
    'laterSess', ~chanBFirstSession, true, 2, 'Later sessions'};
for iSG = 1:size(sessGroups, 1)
    
    % figure: mean ripple rates in trial halves
    x   = [1, 2];
    d   = {chanMeanRippleRateHalves(sessGroups{iSG, 2} == 1, 1), chanMeanRippleRateHalves(sessGroups{iSG, 2} == 1, 2)};
    m   = {mean(d{1}, 1), mean(d{2}, 1)};
    ste = {LK_ste(d{1}), LK_ste(d{2})};
    
    % stats: comparison of ripple rates between first and second half
    [~, pRRHalves, ~, statsRRHalves]    = ttest(d{1}, d{2});
    if sessGroups{iSG, 3} == true
        pRRHalves   = pRRHalves * sessGroups{iSG, 4}; % Bonferroni correction
    end
    fprintf('Are the ripple rates different between the first and second half of all trials? t(%d) = %.3f, p = %.3f.\n', statsRRHalves.df, statsRRHalves.tstat, pRRHalves);
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 5, 6]);
    axes('units', 'centimeters', 'position', [1.5, 1.4, 3, 3.5]);
    hold on;
    for iM = 1:numel(m)
        bar(x(iM), m{iM}, 'EdgeColor', [0, 0, 0], 'FaceColor', [0.9, 0.9, 0.9]);
        plot([x(iM), x(iM)], [m{iM} - ste{iM}, m{iM} + ste{iM}], 'Color', [0, 0, 0], 'LineWidth', 2);
        plot(x(iM) + rand(size(d{iM}, 1), 1) .* 0.5 - 0.25, d{iM}, '.', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 3); % individual data points
    end
    xl = xlabel('trials');
    yl = ylabel('Ripple rate (Hz)');
    t1 = text(0.05, 0.9, sprintf('\\itP\\rm = %.3f', pRRHalves), 'units', 'normalized');
    tl = title({sessGroups{iSG, 5}, sprintf('(%d channels)', size(d{1}, 1))});
    set(gca, 'xlim', [min(x) - 0.8, max(x) + 0.8], 'xtick', x, 'xticklabel', {'Early', 'Late'}, 'ylim', [0, 0.3], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    LK_print(f, strcat(paths.save, 'chanRippleRate_ComparisonBetweenSessionHalves_', sessGroups{iSG, 1}, '_20220705'), '-dtiff', '-r300');
end

%% absolute timing of ripples, start-locked

% report
fprintf('\nAnalysis of ripple rates during the different trial phases, start-locked.\n');
rng(r.param.myRNG + 1); % reset rng for reproducibility

% settings for this analysis
thisParam           	= [];
thisParam.smoothType    = 'gaussian'; % how to smooth the average ripple rates
thisParam.smoothFac     = 0.5 / r.param.beh.PSTHTimeRes + 1; % smoothing kernel has a duration of 0.5 s
thisParam.xLim          = [-1, 3]; % (s)
thisParam.yLim          = [0, 0.3]; % (Hz)

% loop through trial phases
for iPhase = 1:numel(param.beh.trialPhases)
    
    % report
    fprintf('\n============================================================ Ripple rates, absolute timing, start locked: %s.\n\n', param.beh.trialPhases{iPhase});
    
    % ripple PSTH from this trial phase
    ripplePSTH  = cellfun(@(x) squeeze(x(:, iPhase, :)), {r.allRes.ripplePSTHAbsStart}', 'UniformOutput', false);
    
    % transform the ripple counts into ripple rates
    ripplePSTH  = cellfun(@(x) x ./ r.param.beh.PSTHTimeRes, ripplePSTH, 'UniformOutput', false);
    
    % all ripple PSTHs (for plotting) and associated memory performance
    allPSTH         = cell2mat(ripplePSTH);
    allBGoodMemory  = cell2mat(chanBGoodMem);
    
    % mean PSTH per channel (all, good, bad)
    meanPSTH        = nan(size(ripplePSTH, 1), numel(r.param.beh.PSTHStartTime)); % channels x timepoints
    meanPSTHGood    = nan(size(meanPSTH));
    meanPSTHBad     = nan(size(meanPSTH));
    for iChan = 1:size(ripplePSTH, 1)
        
        % only consider trials phases that are long enough
        bValid  = r.allRes(iChan).durationPerPhase(:, iPhase) > 1; % at least 1 second
        
        % sanity check
        if size(bValid, 1) ~= size(ripplePSTH{iChan, 1}, 1) || size(bValid, 1) ~= size(chanBGoodMem{iChan, 1}, 1)
            error('Problem with the size of "bValid".');
        end
        
        % mean PSTH from all trials (ignoring trials that are too short)
        meanPSTH(iChan, :)      = mean(ripplePSTH{iChan, 1}(bValid, :), 'omitnan');
        
        % mean PSTH from good trials
        meanPSTHGood(iChan, :)  = mean(ripplePSTH{iChan, 1}(chanBGoodMem{iChan, 1} == 1 & bValid, :), 'omitnan');
        
        % mean PSTH from bad trials
        meanPSTHBad(iChan, :)   = mean(ripplePSTH{iChan, 1}(chanBGoodMem{iChan, 1} == 0 & bValid, :), 'omitnan');
    end
    
    %% smoothing
    
    % smooth over time
    meanPSTH        = smoothdata(meanPSTH, 2, thisParam.smoothType, thisParam.smoothFac);
    meanPSTHGood    = smoothdata(meanPSTHGood, 2, thisParam.smoothType, thisParam.smoothFac);
    meanPSTHBad     = smoothdata(meanPSTHBad, 2, thisParam.smoothType, thisParam.smoothFac);
    
    %% significance testing
    
    % re-organize data for fieldtrip
    timelock1   = cell(1, size(meanPSTHGood, 1));
    timelock2  	= cell(1, size(meanPSTHBad, 1));
    for iChan = 1:size(timelock1, 2)
        
        % good trials
        timelock1{iChan}.avg       = meanPSTHGood(iChan, :);
        timelock1{iChan}.label     = {'chan'};
        timelock1{iChan}.time      = r.param.beh.PSTHStartTime;
        
        % bad trials
        timelock2{iChan}.avg       = meanPSTHBad(iChan, :);
        timelock2{iChan}.label     = {'chan'};
        timelock2{iChan}.time      = r.param.beh.PSTHStartTime;
    end
    
    % fieldtrip configuration
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = thisParam.xLim; % time window for significance testing (same window as for plotting)
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'depsamplesT';
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum';
    cfg.neighbours          = [];
    cfg.tail                = 0;
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = r.param.ana.numSurrogates;
    numChans                = size(timelock1, 2);
    design                  = zeros(2, numChans * 2);
    design(1, :)            = [1:numChans, 1:numChans]; % channel indices
    design(2, :)            = [ones(1, numChans), ones(1, numChans) * 2]; % 1 = good; 2 = bad
    cfg.design              = design; % design matrix
    cfg.uvar                = 1; % row of the design matrix that contains the units of observation
    cfg.ivar                = 2; % row of the design matrix that contains the independent variable
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
        
    % report fieldtrip significance
    LK_report_ft_timelockstatistics(outFT);
    
    %% figure
    
    % create figure
    dt                  = [];
    dt.beh              = r.param.beh;
    dt.thisParam        = thisParam;
    dt.iPhase           = iPhase;
    dt.allBGoodMemory   = allBGoodMemory;
    dt.meanPSTHBad      = meanPSTHBad;
    dt.meanPSTHGood     = meanPSTHGood;
    dt.outFT            = outFT;
    dt.allPSTH          = allPSTH;
    f = LK_PlotRipplePSTHPerPhaseStart_20231030(dt);
    
    % save figure
    LK_print(f, strcat(paths.save, 'chanRippleRate_AbsTimeFromStart_', param.beh.trialPhases{iPhase}), '-dtiff', '-r450');

    % save figure data
    if strcmp(r.param.beh.trialPhases{iPhase}, 'Cue')
        save(strcat(r.paths.save, 'Fig_3c_Cue'), 'dt');
    elseif strcmp(r.param.beh.trialPhases{iPhase}, 'Feedback')
        save(strcat(r.paths.save, 'Fig_3c_Feedback'), 'dt');
    end
end

%% absolute timing of ripples, end-locked

% report
fprintf('\nAnalysis of ripple rates during the different trial phases, end-locked.\n');
rng(r.param.myRNG + 2); % reset rng for reproducibility

% settings for this analysis
thisParam           	= [];
thisParam.smoothType    = 'gaussian';
thisParam.smoothFac     = 0.5 / r.param.beh.PSTHTimeRes + 1;
thisParam.xLim          = [-3, 1]; % (s)
thisParam.yLim          = [0, 0.3]; % (Hz)

% loop through trial phases
for iPhase = 1:numel(param.beh.trialPhases)
    
    % report
    fprintf('\n============================================================ Ripple rates, absolute timing, end locked: %s.\n\n', param.beh.trialPhases{iPhase});
    
    % ripple PSTH from this trial phase
    ripplePSTH  = cellfun(@(x) squeeze(x(:, iPhase, :)), {r.allRes.ripplePSTHAbsEnd}', 'UniformOutput', false);
    
    % transform the ripple counts into ripple rates
    ripplePSTH  = cellfun(@(x) x ./ r.param.beh.PSTHTimeRes, ripplePSTH, 'UniformOutput', false);
    
    % all ripple PSTHs (for plotting) and associated memory performance
    allPSTH         = cell2mat(ripplePSTH);
    allBGoodMemory  = cell2mat(chanBGoodMem);
    
    % mean PSTH per channel (all, good, bad)
    meanPSTH        = nan(size(ripplePSTH, 1), numel(r.param.beh.PSTHEndTime)); % channels x timepoints
    meanPSTHGood    = nan(size(meanPSTH));
    meanPSTHBad     = nan(size(meanPSTH));
    for iChan = 1:size(ripplePSTH, 1)
        
        % only consider trial phases that are long enough
        bValid  = r.allRes(iChan).durationPerPhase(:, iPhase) > 1; % at least one second
        
        % sanity check
        if size(bValid, 1) ~= size(ripplePSTH{iChan, 1}, 1) || size(bValid, 1) ~= size(chanBGoodMem{iChan, 1}, 1)
            error('Problem with the size of "bValid".');
        end
        
        % mean PSTH from all trials (ignoring trials that are too short)
        meanPSTH(iChan, :)      = mean(ripplePSTH{iChan, 1}(bValid, :), 'omitnan');
        
        % mean PSTH from good trials
        meanPSTHGood(iChan, :)  = mean(ripplePSTH{iChan, 1}(chanBGoodMem{iChan, 1} == 1 & bValid, :), 'omitnan');
        
        % mean PSTH from bad trials
        meanPSTHBad(iChan, :)   = mean(ripplePSTH{iChan, 1}(chanBGoodMem{iChan, 1} == 0 & bValid, :), 'omitnan');
    end
    
    %% smoothing
    
    % smooth over time
    meanPSTH        = smoothdata(meanPSTH, 2, thisParam.smoothType, thisParam.smoothFac);
    meanPSTHGood    = smoothdata(meanPSTHGood, 2, thisParam.smoothType, thisParam.smoothFac);
    meanPSTHBad     = smoothdata(meanPSTHBad, 2, thisParam.smoothType, thisParam.smoothFac);
    
    %% significance testing
    
    % re-organize data for fieldtrip
    timelock1   = cell(1, size(meanPSTHGood, 1));
    timelock2  	= cell(1, size(meanPSTHBad, 1));
    for iChan = 1:size(timelock1, 2)
        
        % good trials
        timelock1{iChan}.avg       = meanPSTHGood(iChan, :);
        timelock1{iChan}.label     = {'chan'};
        timelock1{iChan}.time      = r.param.beh.PSTHEndTime;
        
        % bad trials
        timelock2{iChan}.avg       = meanPSTHBad(iChan, :);
        timelock2{iChan}.label     = {'chan'};
        timelock2{iChan}.time      = r.param.beh.PSTHEndTime;
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
    cfg.tail                = 0;
    cfg.alpha               = 0.05;
    cfg.correcttail         = 'alpha';
    cfg.numrandomization    = r.param.ana.numSurrogates;
    numChans                = size(timelock1, 2);
    design                  = zeros(2, numChans * 2);
    design(1, :)            = [1:numChans, 1:numChans];
    design(2, :)            = [ones(1, numChans), ones(1, numChans) * 2];
    cfg.design              = design;
    cfg.uvar                = 1;
    cfg.ivar                = 2;
    
    % fieldtrip estimation
    outFT = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
        
    % report fieldtrip significance
    LK_report_ft_timelockstatistics(outFT);
    
    %% figure
    
    % create figure
    dt                  = [];
    dt.beh              = r.param.beh;
    dt.thisParam        = thisParam;
    dt.iPhase           = iPhase;
    dt.allBGoodMemory   = allBGoodMemory;
    dt.meanPSTHBad      = meanPSTHBad;
    dt.meanPSTHGood     = meanPSTHGood;
    dt.outFT            = outFT;
    dt.allPSTH          = allPSTH;
    f = LK_PlotRipplePSTHPerPhaseEnd_20231030(dt);
    
    % save figure
    LK_print(f, strcat(paths.save, 'chanRippleRate_AbsTimeFromEnd_', param.beh.trialPhases{iPhase}), '-dtiff', '-r450');

    % save figure data
    if strcmp(r.param.beh.trialPhases{iPhase}, 'ITI')
        save(strcat(r.paths.save, 'Fig_3c_ITI'), 'dt');
    elseif strcmp(r.param.beh.trialPhases{iPhase}, 'Retrieval')
        save(strcat(r.paths.save, 'Fig_3c_Retrieval'), 'dt');
    elseif strcmp(r.param.beh.trialPhases{iPhase}, 'Reencoding')
        save(strcat(r.paths.save, 'Fig_3c_Reencoding'), 'dt');    
    end
end

%% total fraction of IEDs (or artifacts) per channel

% total fraction of IEDs for each channel
chanTotFracIEDs         = cell2mat({r.allRes.totFracIEDs}');
LK_ReportMeanAndSEM_20220322('\nIED fraction (%%)', chanTotFracIEDs .* 100);

% total fraction of artifacts for each channel
chanTotFracArtifacts    = cell2mat({r.allRes.totFracArtifacts}');
LK_ReportMeanAndSEM_20220322('Artifact fraction (%%)', chanTotFracArtifacts .* 100);

% overall relationship between artifact fractions and ripple rates
[rho, pval] = corr(chanTotFracIEDs, chanMeanRippleRate, 'type', 'spearman');
fprintf('Relationship between "chanTotFracArtifacts" and "chanMeanRippleRate": rho = %.3f, p = %.3f.\n', rho, pval);
[rho, pval] = corr(chanTotFracIEDs, chanMeanRippleRateArtFree, 'type', 'spearman');
fprintf('Relationship between "chanTotFracArtifacts" and "chanMeanRippleRateArtFree": rho = %.3f, p = %.3f.\n', rho, pval);

%% effect of trial index on IED fraction

% correlation between artifact fraction and trial index
rhoIED2TrialIdx = nan(size(r.allRes, 1), 1);
for iChan = 1:size(r.allRes, 1)
    thisTrialIdx                = transpose(1:size(r.allRes(iChan).fracIEDsPerTrial, 1));
    rhoIED2TrialIdx(iChan, 1)   = corr(thisTrialIdx, r.allRes(iChan).fracIEDsPerTrial, 'type', 'pearson', 'rows', 'complete');
end

% statistics: test general trend of correlation values
[~, p, ~, stats]    = ttest(rhoIED2TrialIdx);
fprintf('Are the correlation values in "rhoIED2TrialIdx" generelly above zero? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% create histogram figure: correlation coefficients
cfg             = [];
cfg.data        = rhoIED2TrialIdx;
cfg.binEdges    = -1:0.1:1;
cfg.xlabel      = 'Pearson''s \itr';
cfg.ylabel      = 'Count';
cfg.title       = {['(', LK_IndicateThousands(sum(~isnan(rhoIED2TrialIdx))), ' channels)']};
cfg.fontSize    = 0.4;
cfg.figPosition = [2, 2, 6, 6];
cfg.axPosition  = [1.75, 1.4, 3.5, 3.5];
f               = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'rhoIED2TrialIdx_20220321'), '-dtiff', '-r300');
