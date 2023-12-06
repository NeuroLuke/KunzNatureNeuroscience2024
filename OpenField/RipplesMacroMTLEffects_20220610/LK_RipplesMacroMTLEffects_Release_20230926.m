%==========================================================================
% This script investigates the MTL-wide effects of hippocampal ripples
% including cross-correlations and time-frequency power analyses.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all; clc;

% settings
param                           = [];
param.myRNG                     = 555;
%- subjects
param.subject.fileName          = 'subjectdata_20210616.mat';
%- preprocessing
param.preproc.type              = 'FBS'; % filtered, bipolar, select
%- LFPs
param.LFP.choi                  = {'AMYleftBipolar', 'AMYrightBipolar', 'ECleftBipolar', 'ECrightBipolar', 'HCleftBipolar', 'HCrightBipolar', 'PHCleftBipolar', 'PHCrightBipolar', 'TPleftBipolar', 'TPrightBipolar'};
param.LFP.method                = 'wavelet';
param.LFP.foi                   = logspace(log10(1), log10(200), 50);
param.LFP.pad                   = 'nextpow2';
param.LFP.width                 = 7;
param.LFP.output                = 'pow';
param.LFP.toi                   = 'all';
param.LFP.twoi4Data             = [-3, 3]; % time window for data extraction around ripples
param.LFP.fsample               = 2000; % sampling frequency
%- artifacts
param.artifact.time2Exclude     = 1; % additional time to exclude around artifacts (seconds)
param.artifact.amp.idx          = 1; % index for labeling
param.artifact.amp.threshFac    = 4; % factor to calculate amplitude threshold
param.artifact.gra.idx         	= 2; % index for labeling
param.artifact.gra.threshFac    = 4; % factor to calculate gradient threshold
param.artifact.pow.idx          = 3; % index for labeling
param.artifact.pow.threshFac    = 4; % factor to calculate power threshold
param.artifact.pow.lowerFreq    = 1; % lower frequency for calculating the power spectrum
param.artifact.pow.upperFreq    = 60; % upper frequency for calculating the power spectrum
param.artifact.pow.numFreqs     = 30; % number of frequencies for calculating the power spectrum
param.artifact.bInspect         = false; % whether to inspect the artifacts
%- ripples
param.ripple.type               = 'ripples';
param.ripple.fileName           = strcat(param.ripple.type, '.mat');
param.ripple.timepointType      = 'peak';
param.ripple.bpfilter           = 'yes'; % whether to use a bandpass filter
param.ripple.bpfilttype         = 'but'; % filter type
param.ripple.bpfreq             = [80, 140];
param.ripple.envSmoothTime      = 0.02; % time window for smoothing the envelope
param.ripple.threshMinFac       = 2;
param.ripple.threshMaxFac       = Inf;
param.ripple.threshPeakMinFac   = 3;
param.ripple.durationMin        = 0.02; % minimum ripple duration
param.ripple.durationMax        = 0.5; % maximum ripple duration
param.ripple.numCyclesMin       = 3; % minimum number of cycles
param.ripple.twoi4Data          = [-3, 3]; % time window for data extraction around ripples
param.ripple.twoi4BaselineERP   = [-3, 3]; % time window for baseline correction of the ripple ERP
%- for false-positive detection
param.ripple.falsePos.method    = 'wavelet';
param.ripple.falsePos.width     = 7;
param.ripple.falsePos.pad       = 'nextpow2';
param.ripple.falsePos.output    = 'pow';
param.ripple.falsePos.foi     	= 30:2:190; % frequencies for calculating the power spectrum (Hz)
param.ripple.falsePos.timeRes   = 'all';
%- for trial-phase specific analyses
param.trial.phases              = {'AllPhases'; 'ITI'; 'Cue'; 'Retrieval'; 'Feedback'; 'Reencoding'}; % trial phases
param.trial.maxNum              = 160; % maximum number of trials
param.trial.speedCutoff         = 0.001; % speed cutoff for defining movement periods
%- cross-correlation
param.xcorr.maxLagTime          = 5; % maximum time lag for the cross-correlation (seconds)
%- region colors in figures
param.fig.colorAMY              = [255, 0, 0] ./ 255; % red
param.fig.colorEC               = [231, 182, 100] ./ 255; % yellow
param.fig.colorHC               = [0, 119, 0] ./ 255; % dark green
param.fig.colorPHC              = [131, 217, 131] ./ 255; % light green
param.fig.colorTP               = [146, 146, 158] ./ 255; % gray
param.fig.colorAll              = [0, 0, 0] ./ 255; % black

% paths
paths           = [];
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.beh       = 'E:\OpenField\BehPreprocessing_20210707\';
paths.macro     = 'E:\OpenField\MacroPreprocessing_20210627\';
paths.ripple    = strcat('E:\OpenField\RipplesMacro_20210614\20220320_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '\');
paths.save      = strcat('E:\OpenField\RipplesMacroMTLEffects_20220610\20230415_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '_', param.ripple.type, '_', param.ripple.timepointType, 'Time', filesep);
paths.functions = 'E:\OpenField\Functions\';
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';
mkdir(paths.save);

% add functions and fieldtrip
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
ft_warning off;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects   	= subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings
save(strcat(paths.save, 'settings'));

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\n==========================================================\n');
    fprintf('SUBJECT: %s.\n', subjects{iSub});
    tic;
    
    % subject data
    subjectdata = load(strcat(paths.subjects, subjects{iSub}, filesep, param.subject.fileName));
    subjectdata = subjectdata.subjectdata;
    
    %% ripple data from hippocampal channels
    
    % hippocampal channels from which ripples were extracted
    HCrippleChans   = LK_dir(strcat(paths.ripple, subjects{iSub}, '\*'));
    
    % assemble the ripple data from all hippocampal ripple channels for
    % later use
    HCrippleData    = repmat(struct(), size(HCrippleChans, 1), 1);
    for iHCRippleChan = 1:size(HCrippleChans, 1)
        
        % report
        fprintf('\tHippocampal ripple channel: %s.\n', HCrippleChans(iHCRippleChan).name);
        
        % load ripples and extract their information
        r   = load(fullfile(HCrippleChans(iHCRippleChan).folder, HCrippleChans(iHCRippleChan).name, param.ripple.fileName));
        
        % collect the information from all hippocampal ripple channels
        HCrippleData(iHCRippleChan).ripples       = r.ripples.ripples; % structure with various information about each ripple
        HCrippleData(iHCRippleChan).eegRipples    = r.ripples.eegRipples; % fieldtrip structure with EEG data from around each ripple
        HCrippleData(iHCRippleChan).bRipple       = r.ripples.bRipple; % boolean vector indicating the ripple timepoints
    end
    
    %% LFP channels of interest
    
    % LFP channels of interest
    LFPChans    = [];
    for iChoi = 1:numel(param.LFP.choi)
        switch param.LFP.choi{iChoi}
            case 'AMYleftBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.AMYleftBipolar}); % AMY
            case 'AMYrightBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.AMYrightBipolar});
            case 'ECleftBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.ECleftBipolar}); % EC
            case 'ECrightBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.ECrightBipolar});
            case 'HCleftBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.HCleftBipolar}); % HC
            case 'HCrightBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.HCrightBipolar});
            case 'PHCleftBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.PHCleftBipolar}); % PHC
            case 'PHCrightBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.PHCrightBipolar});
            case 'TPleftBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.TPleftBipolar}); % TP
            case 'TPrightBipolar'
                LFPChans    = cat(1, LFPChans, {param.LFP.choi{iChoi}, subjectdata.macro.TPrightBipolar});
        end
    end
    LFPChans        = LFPChans(~cellfun(@isempty, LFPChans(:, 2)), :); % remove channels-of-interest that the subject does not have
    LFPChans(:, 2)  = cellfun(@(x) x{:}, LFPChans(:, 2), 'UniformOutput', false); % simplify the channel labels
    
    % report
    fprintf('\nSubject "%s" has --- %d --- LFP channels to examine during hippocampal ripples:\n', subjects{iSub}, size(LFPChans, 1));
    disp(LFPChans);
    
    %% skip subject if the results are already available
    
    % identify previous result files
    resFiles = dir(strcat(paths.save, subjects{iSub}, '*_res*'));
    if size(resFiles, 1) == (size(HCrippleChans, 1) * size(LFPChans, 1))
        fprintf(2, '\nSkipping subject "%s" because --- %d --- result files have already been produced.\n\n', subjects{iSub}, size(resFiles, 1));
        continue;
    end
    
    %% LFP data
    
    % load entire LFP data from this subject
    eegAll      = load(strcat(paths.macro, subjects{iSub}, '\eeg', param.preproc.type, '.mat'));
    eegAll      = eegAll.eegFBS;
    
    % load the grand-average signal
    eegGA       = load(strcat(paths.macro, subjects{iSub}, '\eegFAS.mat')); % filtered, grand average, select
    eegGA     	= eegGA.eegFAS;
    
    %% behavioral data
    
    % load trial information
    trialInfo   = load(strcat(paths.beh, subjects{iSub}, '\trialInfo.mat'));
    trialInfo 	= trialInfo.trialInfoMacrotime; % in MACROTIME
    
    % restrict number of trials
    if size(trialInfo, 1) > param.trial.maxNum
        trialInfo   = trialInfo(1:param.trial.maxNum, :);
    end
    fprintf('\tNumber of cue periods: %d.\n', sum(~isnan(trialInfo.Cue)));
    
    % start and end time of the experiment
    expStartEndTime     = [min(trialInfo{1, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}}), ...
        max(trialInfo{end, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}})];
    fprintf('\tExperiment start and end time: %.1f s ... %.1f s.\n', min(expStartEndTime), max(expStartEndTime));
    
    % timepoints during the experiment
    bExperiment         = eegAll.time{1} >= min(expStartEndTime) & eegAll.time{1} <= max(expStartEndTime);
    
    % load behavioral information
    behInfo = load(strcat(paths.beh, subjects{iSub}, '\behInfoRes10Hz.mat')); % use resampled behavioral data
    behInfo = behInfo.behInfoMacrotimeRes; % in MACROTIME
    
    %% artifact rejection: ripples in the grand-average signal
    
    % detect artifacts
    cfg                 = param.artifact;
    artifactsGA         = LK_ArtifactDetection_20220319(cfg, eegGA);
    
    % preprocessing of data for ripple detection
    cfg                 = param.ripple;
    cfg.bArtifact       = artifactsGA.bArtifact; % remove artifacts
    data4RipplesGA      = LK_RippleDetection_Preprocessing_20220319(cfg, eegGA);
    
    % detection of ripples, results in "idxRipples" and "idxFalsePositives"
    cfg                 = param.ripple;
    ripplesGA           = LK_RippleDetection_DetectRipples_20220319(cfg, data4RipplesGA);
    
    % characterization of false positives
    cfg                 = [];
    cfg.idxRipples      = ripplesGA.idxFalsePositives; % m x 2 matrix with onsets and ends of false positives
    falsePositivesGA    = LK_RippleDetection_CharacterizeRipples_20220319(cfg, data4RipplesGA);
    
    % data extraction for false positives
    cfg                 = param.ripple;
    cfg.ripples         = falsePositivesGA.ripples; % m x 1 structure with multiple fields for characterizing ripples
    falsePositivesGA    = LK_RippleDetection_ExtractData_20220319(cfg, data4RipplesGA);
    
    % characterization of ripples
    cfg                 = [];
    cfg.idxRipples      = ripplesGA.idxRipples; % m x 2 matrix with onsets and ends of ripples
    ripplesGA           = LK_RippleDetection_CharacterizeRipples_20220319(cfg, data4RipplesGA);
    
    % data extraction for ripples
    cfg                 = param.ripple;
    cfg.ripples         = ripplesGA.ripples; % m x 1 structure with multiple fields for characterizing ripples 
    ripplesGA           = LK_RippleDetection_ExtractData_20220319(cfg, data4RipplesGA);
    
    % save grand-average ripples
 	save(strcat(paths.save, subjects{iSub}, '_GrandAverage_ripples.mat'), 'ripplesGA', 'falsePositivesGA');
        
    %% loop through LFP channels of interest
    for iLFPChan = 1:size(LFPChans, 1)
        
        %% select data from this LFP channel of interest
        
        % macro EEG data from this LFP channel
        cfg         = [];
        cfg.channel = LFPChans(iLFPChan, 2);
        eegLFP      = ft_selectdata(cfg, eegAll);
        
        % report
        fprintf('\n\n------------------------------------------------------\n');
        fprintf('Working on "%s", channel "%s".\n', subjects{iSub}, eegLFP.label{:});
        fprintf('------------------------------------------------------\n\n');
        
        %% artifact rejection: epileptic activity
        
        % detect epileptic artifacts in this channel
        cfg         	= param.artifact;
        LFPartifacts   	= LK_ArtifactDetection_20220319(cfg, eegLFP);
                        
        %% ripple detection on this LFP channel
        
        % preprocessing of data for ripple detection
        cfg           	= param.ripple;
        cfg.bArtifact 	= LFPartifacts.bArtifact | ripplesGA.bRipple | falsePositivesGA.bRipple; % remove epileptic artifacts, grand-average ripples, and grand-average false positives
        data4Ripples   	= LK_RippleDetection_Preprocessing_20220319(cfg, eegLFP);
        
        % detection of ripples, results in "idxRipples"
        cfg         	= param.ripple;
        LFPripples   	= LK_RippleDetection_DetectRipples_20220319(cfg, data4Ripples);
        
        % characterization of ripples
        cfg           	= [];
        cfg.idxRipples 	= LFPripples.idxRipples;
        LFPripples     	= LK_RippleDetection_CharacterizeRipples_20220319(cfg, data4Ripples);
                
        % data extraction for ripples
        cfg           	= param.ripple;
        cfg.ripples    	= LFPripples.ripples;
        LFPripples     	= LK_RippleDetection_ExtractData_20220319(cfg, data4Ripples);
        
        % save LFP ripples
        save(strrep(strcat(paths.save, subjects{iSub}, '_', LFPChans{iLFPChan, 2}, '_ripples.mat'), '-', ''), 'LFPripples');
                
        %% power over time
        
        % time-frequency power across the entire recording
        cfg       	= [];
        cfg.method  = param.LFP.method;
        cfg.width 	= param.LFP.width;
        cfg.foi   	= param.LFP.foi;
        cfg.pad    	= param.LFP.pad;
        cfg.output  = param.LFP.output;
        cfg.toi    	= param.LFP.toi;
        LFPtfData 	= ft_freqanalysis(cfg, eegLFP);
        
        % power at all time points
        LFPpow   	= squeeze(LFPtfData.powspctrm); % freq x time
        
        % sanity check
        if size(LFPpow, 2) ~= size(LFPartifacts.bArtifact, 2)
            error('"LFPpow" does not have the same size as "bArtifact".');
        end
        
        % set artifact timepoints to nan
        LFPpow(:, LFPartifacts.bArtifact)   = nan; % row 1, lowest frequency; row n, highest frequency
        
        % normalize power
        LFPpowNorm  = normalize(LFPpow, 2, 'zscore'); % row 1, lowest frequency; row n, highest frequency
                
        %% MTL effects of hippocampal ripples
        
        % you obtain results for each combination of an LFP channel and its
        % simultaneously recorded hippocampal ripple channels
        for iHCRippleChan = 1:size(HCrippleChans, 1)
            
            % file name for saving
            saveFileName        = strrep(strcat(paths.save, subjects{iSub}, '_lfp', LFPChans{iLFPChan, 2}, '_rip', HCrippleChans(iHCRippleChan).name, '_res.mat'), '-', '');
            
            % timepoints of all ripples of this channel (in MACROTIME)
            HCrippleTimepoints  = cell2mat({HCrippleData(iHCRippleChan).ripples.peakTime}'); % cave: some ripples may be outside the experiment
            
            % skip if there are no ripples
            if size(HCrippleTimepoints, 1) < 1
                fprintf(2, '\tSkipping --- %s (%s) --- because there are no ripples in the hippocampus.\n', subjects{iSub}, HCrippleChans(iHCRippleChan).name);
                continue;
            end
            
            % assign ripples to trial phases
            HCrippleTrial   = nan(size(HCrippleTimepoints, 1), 1);
            for iTP = 1:numel(param.trial.phases)
                
                % starts and ends of this trial phase
                if strcmp(param.trial.phases{iTP}, 'AllPhases')
                    startEnd    = [trialInfo.ITI, trialInfo.Grab];
                elseif strcmp(param.trial.phases{iTP}, 'ITI')
                    startEnd    = [trialInfo.ITI, trialInfo.Cue];
                elseif strcmp(param.trial.phases{iTP}, 'Cue')
                    startEnd    = [trialInfo.Cue, trialInfo.Retrieval];
                elseif strcmp(param.trial.phases{iTP}, 'Retrieval')
                    startEnd    = [trialInfo.Retrieval, trialInfo.Feedback];
                elseif strcmp(param.trial.phases{iTP}, 'Feedback')
                    startEnd    = [trialInfo.Feedback, trialInfo.Reencoding];
                elseif strcmp(param.trial.phases{iTP}, 'Reencoding')
                    startEnd    = [trialInfo.Reencoding, trialInfo.Grab];
                end
                
                % ripples from this trial phase
                bPhase                  = any(HCrippleTimepoints >= startEnd(:, 1)' & HCrippleTimepoints <= startEnd(:, 2)', 2);
                HCrippleTrial(bPhase)   = iTP;
            end
            
            %% cross-correlation between ripple timelines, overall
            % use an "unbiased" cross-correlation that accounts for
            % different numbers of samples at different lags
            
            % compute cross correlations at reasonable lags
            [cc, ccLags]    = xcorr(HCrippleData(iHCRippleChan).bRipple, LFPripples.bRipple, param.xcorr.maxLagTime * LFPripples.eegRipples.fsample, 'unbiased');
            
            % set the cross correlations to NaN if there are no ripples
            if sum(HCrippleData(iHCRippleChan).bRipple) == 0 || sum(LFPripples.bRipple) == 0
                cc  = nan(size(cc));
                warning('Setting the cross correlation "cc" to NaNs because there are no ripples.');
            end
            
            %% cross-correlation between ripple timelines, separately for the different trial phases
            % note: calculate it separately for each trial, then average
            % across trials
            
            % preallocate
            ccTrial     = nan(size(param.trial.phases, 1), param.xcorr.maxLagTime * LFPripples.eegRipples.fsample * 2 + 1, size(trialInfo, 1)); % trial phases x lags x trials
            ccLagsTrial = nan(size(ccTrial));
            
            % loop through the different trial phases
            for iTP = 1:size(param.trial.phases, 1)
                
                % starts and ends of this trial phase
                if strcmp(param.trial.phases{iTP}, 'AllPhases')
                    startEnd    = [trialInfo.ITI, trialInfo.Grab];
                elseif strcmp(param.trial.phases{iTP}, 'ITI')
                    startEnd    = [trialInfo.ITI, trialInfo.Cue];
                elseif strcmp(param.trial.phases{iTP}, 'Cue')
                    startEnd    = [trialInfo.Cue, trialInfo.Retrieval];
                elseif strcmp(param.trial.phases{iTP}, 'Retrieval')
                    startEnd    = [trialInfo.Retrieval, trialInfo.Feedback];
                elseif strcmp(param.trial.phases{iTP}, 'Feedback')
                    startEnd    = [trialInfo.Feedback, trialInfo.Reencoding];
                elseif strcmp(param.trial.phases{iTP}, 'Reencoding')
                    startEnd    = [trialInfo.Reencoding, trialInfo.Grab];
                end
                
                % loop through the different trials
                for iTrial = 1:size(startEnd, 1)
                    
                    % identify all timepoints belonging to this trial phase
                    bPhase  = eegLFP.time{1} >= startEnd(iTrial, 1) & eegLFP.time{1} <= startEnd(iTrial, 2);
                    
                    % compute cross-correlation at reasonable lags for this
                    % trial phase
                    [thisCC, thisCCLags]  = xcorr(HCrippleData(iHCRippleChan).bRipple(bPhase), LFPripples.bRipple(bPhase), param.xcorr.maxLagTime * LFPripples.eegRipples.fsample, 'unbiased');
                    
                    % set the cross correlation to NaN if there are no
                    % ripples
                    if sum(HCrippleData(iHCRippleChan).bRipple(bPhase)) == 0 || sum(LFPripples.bRipple(bPhase)) == 0
                        thisCC              = nan(size(thisCC));
                    end
                    
                    % collect across phases and trials
                    ccTrial(iTP, :, iTrial)     = thisCC;
                    ccLagsTrial(iTP, :, iTrial) = thisCCLags;
                end
            end
            
            % average across trials
            ccTrial     = mean(ccTrial, 3, 'omitnan'); % trial phases x lags
            ccLagsTrial = mean(ccLagsTrial, 3, 'omitnan');
            
            %% power spectrum: normalized LFP power during the hippocampal ripples
            
            % average power in this channel during ripples (artifacts and
            % timepoints outside the experiment are nanned out)
            avLFPpowNormDuringHCRipples = mean(LFPpowNorm(:, HCrippleData(iHCRippleChan).bRipple), 2, 'omitnan'); % during hippocampal ripples
            
            %% power spectrum: normalized LFP power during the hippocampal ripples, separately for the different trial phases
            
            % preallocate power spectra during ripples from different trial
            % phases
            avLFPpowNormDuringHCRipplesTrial    = nan(size(LFPpowNorm, 1), size(param.trial.phases, 1)); % frequencies x trial phases
            
            % loop through the different trial phases
            for iTP = 1:size(param.trial.phases, 1)
                
                % starts and ends of this trial phase
                if strcmp(param.trial.phases{iTP}, 'AllPhases')
                    startEnd    = [trialInfo.ITI, trialInfo.Grab];
                elseif strcmp(param.trial.phases{iTP}, 'ITI')
                    startEnd    = [trialInfo.ITI, trialInfo.Cue];
                elseif strcmp(param.trial.phases{iTP}, 'Cue')
                    startEnd    = [trialInfo.Cue, trialInfo.Retrieval];
                elseif strcmp(param.trial.phases{iTP}, 'Retrieval')
                    startEnd    = [trialInfo.Retrieval, trialInfo.Feedback];
                elseif strcmp(param.trial.phases{iTP}, 'Feedback')
                    startEnd    = [trialInfo.Feedback, trialInfo.Reencoding];
                elseif strcmp(param.trial.phases{iTP}, 'Reencoding')
                    startEnd    = [trialInfo.Reencoding, trialInfo.Grab];
                end
            
                % identify the timepoints that belong to this trial phase
                bPhase  = transpose(any(eegLFP.time{1}(:) >= startEnd(:, 1)' & eegLFP.time{1}(:) <= startEnd(:, 2)', 2));
                
                % average power in this channel during ripples (artifacts
                % are nanned out)
                avLFPpowNormDuringHCRipplesTrial(:, iTP)    = mean(LFPpowNorm(:, bPhase & HCrippleData(iHCRippleChan).bRipple), 2, 'omitnan'); % during hippocampal ripples
            end
            
            %% power spectrogram: normalized LFP power within larger time window around each hippocampal ripple
            
            % loop through ripples and collect power values
            LFPtimeDuringHCRipples      = (eegLFP.fsample * min(param.LFP.twoi4Data):eegLFP.fsample * max(param.LFP.twoi4Data)) ./ eegLFP.fsample;
            LFPpowNormDuringHCRipples   = nan(size(HCrippleTimepoints, 1), size(LFPpowNorm, 1), size(LFPtimeDuringHCRipples, 2)); % ripples x freq x time-during-ripples
            for iRipple = 1:size(HCrippleTimepoints, 1)
                
                % indices for this data segment
                [~, minIdx] = min(abs(LFPtfData.time - HCrippleTimepoints(iRipple, 1)));
                idx         = minIdx + eegLFP.fsample * min(param.LFP.twoi4Data):minIdx + eegLFP.fsample * max(param.LFP.twoi4Data);
                
                % sanity check regarding data length
                if numel(idx) ~= size(LFPtimeDuringHCRipples, 2)
                    error('Incorrect size of "idx".');
                end
                
                % skip this hippocampal ripple, if the indices extend
                % beyond the data during the experiment
                if min(idx) < 1 || max(idx) > size(LFPpowNorm, 2)
                    fprintf('\tSkipping hippocampal ripple #%d, because the data indices go beyond the data.\n', iRipple);
                    continue;
                end
                
                % store the normalized power from this data segment
                LFPpowNormDuringHCRipples(iRipple, :, :)    = LFPpowNorm(:, idx); % ripples x freq x time
            end
            
            % get the power spectrograms into fieldtrip format
            LFPtfDuringHCRipples            = [];
            LFPtfDuringHCRipples.label      = LFPChans(iLFPChan, 2);
            LFPtfDuringHCRipples.dimord     = 'chan_freq_time';
            LFPtfDuringHCRipples.freq       = LFPtfData.freq;
            LFPtfDuringHCRipples.time       = LFPtimeDuringHCRipples;
            LFPtfDuringHCRipples.powspctra  = LFPpowNormDuringHCRipples; % all power spectrograms separately
            LFPtfDuringHCRipples.powspctrm  = mean(LFPpowNormDuringHCRipples, 1, 'omitnan'); % average across ripples
            
            % sanity check
            if size(LFPtfDuringHCRipples.powspctra, 1) ~= size(HCrippleTrial, 1)
                error('Problem with the size of the power spectrograms.');
            end
            
            %% hippocampal ripple ERP data (for figures)
            
            % baseline correction
            cfg             = [];
            cfg.baseline    = param.ripple.twoi4BaselineERP;
            HCrippleERP     = ft_timelockbaseline(cfg, HCrippleData(iHCRippleChan).eegRipples);
            
            % averaging
            cfg             = [];
            HCrippleERP     = ft_timelockanalysis(cfg, HCrippleERP);
            
            %% results
            
            % results from this LFP channel
            res                                     = [];
            res.idx                                 = [iSub, iLFPChan, iHCRippleChan]; % subject index, LFP channel index, HC channel index
            res.subjectName                         = subjects{iSub}; % subject name
            res.LFPChanRegion                       = LFPChans{iLFPChan, 1}; % standardized region of this LFP channel
            res.LFPChanName                         = LFPChans{iLFPChan, 2}; % name of this LFP channel
            res.LFPripples                          = LFPripples.ripples; % information about ripples from this LFP channel
            
            % hippocampal ripple information
            res.HCrippleChanName                    = HCrippleChans(iHCRippleChan).name; % name of this ripple channel
            res.HCrippleTimepoints                  = HCrippleTimepoints; % timepoints when the ripples occurred
            res.HCrippleERPTime                     = HCrippleERP.time; % time of the hippocampal ripple ERP
            res.HCrippleERPAvg                      = HCrippleERP.avg; % hippocampal ripple ERP (averaged across all ripples including ripple-wise baseline correction)
            res.HCrippleTrial                       = HCrippleTrial; % trial phase during which the ripple occurs
            
            % cross correlations (general and trial-phase specific)
            res.cc                                  = cc; % cross-correlation between hippocampal ripples and the ripples from this channel
            res.ccLags                              = ccLags ./ eegLFP.fsample; % lags for the cross-correlation (seconds)
            res.ccTrial                             = ccTrial; % cross-correlation between hippocampal ripples and the ripples from this channel, trial-phase specific
            res.ccLagsTrial                         = ccLagsTrial ./ eegLFP.fsample; % lags for the cross-correlation (seconds), trial-phase specific
            
            % power spectra during HC ripples (general and trial-phase
            % specific)
            res.avLFPpowNormDuringHCRipples         = avLFPpowNormDuringHCRipples; % normalized power (z-scored) from this LFP channel during hippocampal ripples (freq x 1 vector)
            res.avLFPpowNormDuringHCRipplesTrial    = avLFPpowNormDuringHCRipplesTrial;
            
            % save this result including the power spectrograms during each
            % HC ripple
            save(saveFileName, 'res', 'LFPtfDuringHCRipples', '-v7.3');
            
            % clear large variables
            clear res LFPtfDuringHCRipples LFPpowNormDuringHCRipples;
            
            %% close open figures
            close all;
        end
    end
    toc;
end

%% collect results

% report
fprintf('\nCollecting previously estimated results:\n');

% previously saved result files
resFiles    = dir(strcat(paths.save, '*_res.mat'));

% loop through result files and assemble results
allRes      = [];
for iFile = 1:size(resFiles, 1)
    
    % report
    if mod(iFile, 20) == 1
        fprintf('\tLoading data from file %d: "%s"\n', iFile, resFiles(iFile).name);
    end
    
    % load data
    res     = load(fullfile(resFiles(iFile).folder, resFiles(iFile).name), 'res');
    res     = res.res;
            
    % collect data across channels
    allRes  = cat(1, allRes, res);
end

%% general information

% all and unique cellular indices
allIdx                      = cell2mat({allRes.idx}'); % all combinations of subject indices, LFP channels, and HC channels
allIdx                      = allIdx(:, 1:2); % all combinations of subject indices and LFP channels
[uniqueIdx, uniqueIdxIA]    = unique(allIdx, 'rows'); % unique combinations of subject indices and LFP channels
allBUnique                  = ismember(transpose(1:size(allIdx, 1)), uniqueIdxIA);

% labels of all hippocampal ripple channels
allHCRippleChanName         = {allRes.HCrippleChanName}';
uniqueHCRippleChanName      = unique(allHCRippleChanName); % unique labels of the hippocampal channels

% brain region of all LFP channels
allLFPChanRegion            = {allRes.LFPChanRegion}';
uniqueLFPChanRegion         = unique(allLFPChanRegion); % unique brain regions of the LFP channels

% information whether LFP channels and ripple channels were recorded from
% the same hemisphere
bSameHemisphere = false(size(allRes, 1), 1);
for iRes = 1:size(allRes, 1)
    if (contains(allLFPChanRegion{iRes}, 'left') && contains(allHCRippleChanName{iRes}, 'L')) || ...
            (contains(allLFPChanRegion{iRes}, 'right') && contains(allHCRippleChanName{iRes}, 'R'))
        bSameHemisphere(iRes, 1)  = true;
    end
end

% all ripple averages for plotting
allHCRippleERPAvg   = cell2mat({allRes.HCrippleERPAvg}'); % average hippocampal ripple per channel combination
allHCRippleERPTime  = cell2mat({allRes.HCrippleERPTime}');

% report
fprintf('\n==== Results.\n');
fprintf('Number of session-LFP channel-ripple channel combinations: %d.\n', size(allIdx, 1));
fprintf('Number of unique session-LFP channel combinations (irrespective of the number of ripple channels in a given subject): %d.\n', size(uniqueIdx, 1));
fprintf('Number of LFP channels recorded from the same hemisphere as their hippocampal ripple channels: %d.\n', sum(bSameHemisphere));
fprintf('Names of the hippocampal ripple channels:\n');
disp(uniqueHCRippleChanName);

%% cross-correlations between hippocampal and extrahippocampal ripples

% reset rng
rng(param.myRNG);

% report
fprintf('\n=== Cross-correlations between hippocampal ripples and ripples in various regions:\n');

% settings for this analysis
thisParam               = [];
thisParam.LFPROIs       = {...
    'AMY', {'AMYleftBipolar'; 'AMYrightBipolar'}, param.fig.colorAMY; ...
    'EC', {'ECleftBipolar'; 'ECrightBipolar'}, param.fig.colorEC; ...
    'HC', {'HCleftBipolar'; 'HCrightBipolar'}, param.fig.colorHC; ...
    'PHC', {'PHCleftBipolar'; 'PHCrightBipolar'}, param.fig.colorPHC; ...
    'TP', {'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorTP; ...
    'Non-HC regions', {'AMYleftBipolar'; 'AMYrightBipolar'; 'ECleftBipolar'; 'ECrightBipolar'; 'PHCleftBipolar'; 'PHCrightBipolar'; 'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorAll};
thisParam.twoi          = [-0.5, 0.5]; % time lags of interest
thisParam.fsample       = param.LFP.fsample;
thisParam.smoothDur     = 0.2; % smoothing duration (s)
thisParam.smoothFac     = thisParam.fsample * thisParam.smoothDur + 1; % smoothing factor
thisParam.smoothType    = 'Gaussian'; % type of smoothing
thisParam.trialPhases   = {...
    'AllPhases', [-1, 1]; ... % trial-phase name, y-limits
    'ITI', [-0.5, 0.5]; ...
    'Cue', [-0.5, 0.5]; ...
    'Retrieval', [-0.5, 0.5]; ...
    'Feedback', [-0.5, 0.5]; ...
    'Reencoding', [-0.5, 0.5]};

% loop through the different trial phases
for iTP = 1:size(thisParam.trialPhases, 1)

    % report trial phase
    fprintf('\n--- Trial phase: %s.\n', thisParam.trialPhases{iTP});
    
    % select data
    if strcmp(thisParam.trialPhases{iTP, 1}, 'AllPhases')
        thisCCLags  = cell2mat({allRes.ccLags}');
        thisCC      = cell2mat({allRes.cc}');
    else
        thisCCLags  = cell2mat(cellfun(@(x) x(strcmp(param.trial.phases, thisParam.trialPhases{iTP, 1}), :), {allRes.ccLagsTrial}', 'UniformOutput', false));
        thisCC      = cell2mat(cellfun(@(x) x(strcmp(param.trial.phases, thisParam.trialPhases{iTP, 1}), :), {allRes.ccTrial}', 'UniformOutput', false));
    end
    
    % normalize cross correlations across time lags
    uniqueCCLags    = round(mean(thisCCLags, 1), 4);
    thisCCLags     	= thisCCLags(:, uniqueCCLags >= min(thisParam.twoi) & uniqueCCLags <= max(thisParam.twoi)); % cave: restrict to time window for displaying
    thisNormCC      = smoothdata(thisCC, 2, thisParam.smoothType, thisParam.smoothFac); % smoothing over all time lags
    thisNormCC    	= thisNormCC(:, uniqueCCLags >= min(thisParam.twoi) & uniqueCCLags <= max(thisParam.twoi)); % cave: restrict to time window for displaying
    thisNormCC   	= normalize(thisNormCC, 2, 'zscore'); % z-score across the time window of interest
    
    % loop through LFP-channel regions
    for iRegion = 1:size(thisParam.LFPROIs, 1)
        
        %% select data
        
        % report
        fprintf('\nRegion: %s.\n', thisParam.LFPROIs{iRegion, 1});
        
        % masks for region and hemisphere
        bMask           = ismember(allLFPChanRegion, thisParam.LFPROIs{iRegion, 2}); % LFP channels from this region
        bMaskSameHemi   = bMask & bSameHemisphere; % LFP channels from this region and the same hemisphere
        bMaskDiffHemi   = bMask & ~bSameHemisphere; % LFP channels from this region and the other hemisphere
        
        %% stats: z-scored cross-correlations against 0
        
        % perform permutation test
        cfg                 = [];
        cfg.mat             = thisNormCC(bMask, :); % test cross-correlations from both hemispheres
        cfg.direction       = 'pos';
        cfg.alpha           = 0.05;
        cfg.numSurrogates   = 1001;
        permTest            = LK_PermutationTest_OneSampleAgainstZero_20200708(cfg);
        
        %% figure
        
        % create figure
        dt                  = [];
        dt.meanRippleTime   = mean(allHCRippleERPTime(bMask, :), 1, 'omitnan');
        dt.meanRipple       = mean(allHCRippleERPAvg(bMask, :), 1, 'omitnan');
        dt.LFPROI           = thisParam.LFPROIs{iRegion, 1};
        dt.nSame            = sum(bMaskSameHemi);
        dt.nDiff            = sum(bMaskDiffHemi);
        dt.twoi             = thisParam.twoi;
        dt.x                = mean(thisCCLags(bMask, :), 1);
        dt.mSame            = mean(thisNormCC(bMaskSameHemi, :), 1, 'omitnan');
        dt.steSame          = LK_ste(thisNormCC(bMaskSameHemi, :));
        dt.mDiff            = mean(thisNormCC(bMaskDiffHemi, :), 1, 'omitnan');
        dt.steDiff          = LK_ste(thisNormCC(bMaskDiffHemi, :));
        dt.yLim             = thisParam.trialPhases{iTP, 2};
        dt.logIdxSigClus    = permTest.logIdxSigClus;
        f = LK_PlotCrossCorrelations_20231005(dt);
        
        % save figure
        LK_print(f, strcat(paths.save, 'CrossCorrelation_LFPRipples2HCRipples_', thisParam.trialPhases{iTP, 1}, '_', thisParam.LFPROIs{iRegion, 1}, '_20230409'), '-dpng', '-r300');

        % save figure data
        if strcmp(thisParam.trialPhases{iTP}, 'AllPhases') && strcmp(thisParam.LFPROIs{iRegion}, 'Non-HC regions')
            save(strcat(paths.save, 'Fig_4a'), 'dt');
        end
    end
    
    % close open figures
    close all;
end

%% power spectrograms

% report
fprintf(2, '\n\nCollecting the power spectrograms.\n');

% load or collect the power spectrograms
if exist(strcat(paths.save, 'allRes.mat'), 'file') > 0
    
    % use previous results
    tmp                             = load(strcat(paths.save, 'allRes'));
    allRes                          = tmp.allRes;
    allLFPtfDuringHCRipples         = tmp.allLFPtfDuringHCRipples;
    allLFPtfDuringHCRipplesTrial    = tmp.allLFPtfDuringHCRipplesTrial;
else
    
    % preallocate
    allLFPtfDuringHCRipples         = nan(numel(param.LFP.foi), range(param.LFP.twoi4Data) * param.LFP.fsample + 1, size(resFiles, 1)); % freq x time x channels
    allLFPtfDuringHCRipplesTrial    = nan(numel(param.LFP.foi), range(param.LFP.twoi4Data) * param.LFP.fsample + 1, size(resFiles, 1), numel(param.trial.phases)); % freq x time x channels x phases
    
    % assemble the power spectrograms for specific conditions
    for iFile = 1:size(resFiles, 1)
        
        % report
        if mod(iFile, 20) == 1
            fprintf('\tLoading data from file %d: "%s"\n', iFile, resFiles(iFile).name);
        end
        
        % load data
        res             = load(fullfile(resFiles(iFile).folder, resFiles(iFile).name), 'LFPtfDuringHCRipples');
        thisPowspctrm   = res.LFPtfDuringHCRipples.powspctrm;
        thisPowspctra   = res.LFPtfDuringHCRipples.powspctra; % ripples x freq x time
        bValid          = all(all(~isnan(thisPowspctra), 3), 2); % exclude power spectrograms with any NaNs
        
        % power spectrogram across the entire recording
        allLFPtfDuringHCRipples(:, :, iFile)        = squeeze(thisPowspctrm);
        
        % collect average power spectrograms from different trial phases
        for iTP = 1:numel(param.trial.phases)
            
            % select relevant ripples
            if strcmp(param.trial.phases{iTP}, 'AllPhases')
                bPhase  = ~isnan(allRes(iFile).HCrippleTrial);
            else
                bPhase  = allRes(iFile).HCrippleTrial == find(strcmp(param.trial.phases, param.trial.phases{iTP}));
            end
            
            % process trial-phase-specific power spectrograms
            if sum(bValid & bPhase) >= 1
                allLFPtfDuringHCRipplesTrial(:, :, iFile, iTP)  = squeeze(mean(thisPowspctra(bPhase & bValid, :, :), 1)); % average across ripples
            end
        end
        
        % free memory
        clear res;
    end
    
    % save the assembled results
    save(strcat(paths.save, 'allRes'), 'allRes', 'allLFPtfDuringHCRipples', 'allLFPtfDuringHCRipplesTrial', '-v7.3');
end

%% ripple-related power spectrogram for each region

% reset rng
rng(param.myRNG);

% report
fprintf('\n=== Ripple-related power spectrograms for various regions:\n');

% set rng
rng(param.myRNG);

% settings for this analysis
thisParam               = [];
thisParam.LFPROIs       = {...
    'Non-HC regions', {'AMYleftBipolar'; 'AMYrightBipolar'; 'ECleftBipolar'; 'ECrightBipolar'; 'PHCleftBipolar'; 'PHCrightBipolar'; 'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorAll; ...
    'Non-HC regions IH', {'AMYleftBipolar'; 'AMYrightBipolar'; 'ECleftBipolar'; 'ECrightBipolar'; 'PHCleftBipolar'; 'PHCrightBipolar'; 'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorAll; ...
    'Non-HC regions CH', {'AMYleftBipolar'; 'AMYrightBipolar'; 'ECleftBipolar'; 'ECrightBipolar'; 'PHCleftBipolar'; 'PHCrightBipolar'; 'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorAll; ...    
    'AMY', {'AMYleftBipolar'; 'AMYrightBipolar'}, param.fig.colorAMY; ... % region title, channel labels, color
    'AMY IH', {'AMYleftBipolar'; 'AMYrightBipolar'}, param.fig.colorAMY; ...
    'AMY CH', {'AMYleftBipolar'; 'AMYrightBipolar'}, param.fig.colorAMY; ...
    'EC', {'ECleftBipolar'; 'ECrightBipolar'}, param.fig.colorEC; ...
    'EC IH', {'ECleftBipolar'; 'ECrightBipolar'}, param.fig.colorEC; ...
    'EC CH', {'ECleftBipolar'; 'ECrightBipolar'}, param.fig.colorEC; ...
    'HC', {'HCleftBipolar'; 'HCrightBipolar'}, param.fig.colorHC; ...
    'HC IH', {'HCleftBipolar'; 'HCrightBipolar'}, param.fig.colorHC; ...
    'HC CH', {'HCleftBipolar'; 'HCrightBipolar'}, param.fig.colorHC; ...
    'PHC', {'PHCleftBipolar'; 'PHCrightBipolar'}, param.fig.colorPHC; ...
    'PHC IH', {'PHCleftBipolar'; 'PHCrightBipolar'}, param.fig.colorPHC; ...
    'PHC CH', {'PHCleftBipolar'; 'PHCrightBipolar'}, param.fig.colorPHC; ...
    'TP', {'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorTP; ...
    'TP IH', {'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorTP; ...
    'TP CH', {'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorTP; ...
    };
thisParam.time          = (param.LFP.fsample * min(param.LFP.twoi4Data):param.LFP.fsample * max(param.LFP.twoi4Data)) ./ param.LFP.fsample; % original time of the power spectrograms
thisParam.twoi          = [-0.5, 0.5]; % time window of interest
thisParam.bTwoi         = thisParam.time >= min(thisParam.twoi) & thisParam.time <= max(thisParam.twoi);
thisParam.timeOI        = thisParam.time(thisParam.bTwoi); % time of interest
thisParam.twoiBase      = [-2.5, -0.5]; % baseline time window
thisParam.bTwoiBase     = thisParam.time >= min(thisParam.twoiBase) & thisParam.time <= max(thisParam.twoiBase);
thisParam.freq          = param.LFP.foi;
thisParam.smoothDur     = 0.2; % smooth duration
thisParam.fsample       = param.LFP.fsample; % sampling rate
thisParam.smoothFac     = thisParam.fsample * thisParam.smoothDur + 1; % smoothing factor
thisParam.smoothType    = 'Gaussian'; % smoothing kernel type
thisParam.trialPhases   = {...
    'All', [-0.05, 0.05]; ... % trial-phase name, caxis limits
    'AllPhases', [-0.05, 0.05]; ...
    'ITI', [-0.1, 0.05]; ...
    'Cue', [-0.1, 0.05]; ...
    'Retrieval', [-0.05, 0.05]; ...
    'Feedback', [-0.2, 0.1]; ...
    'Reencoding', [-0.05, 0.05]; ...
    };

% loop through the different trial phases
for iTP = 1:size(thisParam.trialPhases, 1)
    
    % report
    fprintf('\n\nPower spectrograms for phase: "%s".\n', thisParam.trialPhases{iTP, 1});
    
    % select data
    if strcmp(thisParam.trialPhases{iTP}, 'All')
        allPowspctrm    = allLFPtfDuringHCRipples;
    else
        allPowspctrm    = allLFPtfDuringHCRipplesTrial(:, :, :, strcmp(param.trial.phases, thisParam.trialPhases{iTP}));
    end
    
    % smooth across time
    allPowspctrm    = smoothdata(allPowspctrm, 2, thisParam.smoothType, thisParam.smoothFac);
    
    % baseline correction
    if ~ismember(thisParam.trialPhases{iTP}, {'All'})
        fprintf(2, 'Applying baseline correction to the power spectrograms.\n');
        baseVal         = repmat(mean(allPowspctrm(:, thisParam.bTwoiBase, :), 2, 'omitnan'), 1, size(allPowspctrm, 2), 1);
        allPowspctrm    = allPowspctrm - baseVal;
    end
    
    % shorten to time window of interest
    allPowspctrm    = allPowspctrm(:, thisParam.bTwoi, :); % freq x time x combination
    
    % loop through LFP-channel regions
    for iRegion = 1:size(thisParam.LFPROIs, 1)
        
        % select LFP channels from this region
        bMask       = ismember(allLFPChanRegion, thisParam.LFPROIs{iRegion, 2}); % channels from this region
        
        % select specific hemisphere
        if contains(thisParam.LFPROIs{iRegion, 1}, 'IH')
            bMask   = bMask & bSameHemisphere; % LFP channel and HC ripple channel are in the same hemisphere
        elseif contains(thisParam.LFPROIs{iRegion, 1}, 'CH')
            bMask   = bMask & ~bSameHemisphere; % LFP channel and HC ripple channel are in different hemispheres
        end
        
        % report
        fprintf('\n- Working on "%s". Number of channels: %d.\n', thisParam.LFPROIs{iRegion, 1}, sum(bMask));
        
        %% statistical evaluation of the power spectrogram across channels against 0
        
        % name of the file containing the stats from the CBPT
        CBPTFile    = strrep(strcat(paths.save, 'CBPT_powSpectrogramDuringHCRipples_', thisParam.trialPhases{iTP, 1}, '_', thisParam.LFPROIs{iRegion, 1}, '.mat'), ' ', '_');
        
        % load or perform CBPT
        if exist(CBPTFile, 'file')
            
            % load previously performed CBPT
            CBPT    = load(CBPTFile);
            posCBPT = CBPT.posCBPT; % cluster-based permutation test against 0 (positive direction)
            fprintf('\tThe maximum cluster statistic is: cluster-t = %.3f, cluster-p = %.3f.\n', max(posCBPT.clusStat), min(posCBPT.clusStatP));
            negCBPT = CBPT.negCBPT; % cluster-based permutation test against 0 (negative direction)
            fprintf('\tThe maximum cluster statistic is: cluster-t = %.3f, cluster-p = %.3f.\n', min(negCBPT.clusStat), min(negCBPT.clusStatP));
        else
            % cluster-based permutation test against 0 (positive direction)
            cfg                 = [];
            cfg.mat             = allPowspctrm(:, :, bMask);
            cfg.firstalpha      = 0.025; % threshold for identifying clusters
            cfg.secondalpha     = 0.025; % threshold for estimating the significance of clusters
            cfg.direction       = 'pos';
            cfg.numSurrogates   = 1001;
            posCBPT             = LK_3DPermutationTest_OneSampleAgainstZero_20220302(cfg);
            
            % cluster-based permutation test against 0 (negative direction)
            cfg.direction       = 'neg';
            negCBPT             = LK_3DPermutationTest_OneSampleAgainstZero_20220302(cfg);
            
            % save CBPT results
            save(CBPTFile, 'posCBPT', 'negCBPT');
        end
        
        %% power spectrogram from this trial phase and this region

        % create figure
        dt              = [];
        dt.t            = mean(allHCRippleERPTime(bMask, :), 1, 'omitnan');
        dt.m            = mean(allHCRippleERPAvg(bMask, :), 1, 'omitnan');
        dt.LFPROI       = thisParam.LFPROIs{iRegion, 1};
        dt.n            = sum(bMask);
        dt.twoi         = thisParam.twoi;
        dt.mPow         = mean(allPowspctrm(:, :, bMask), 3, 'omitnan');
        dt.timeOI       = thisParam.timeOI;
        dt.freq         = thisParam.freq;
        dt.posCBPT      = posCBPT;
        dt.negCBPT      = negCBPT;
        dt.cbLim        = thisParam.trialPhases{iTP, 2};
        f = LK_PlotPowSpectrogramDuringHCRipples_20231005(dt);
        
        % save figure
        LK_print(f, strcat(paths.save, 'powSpectrogramDuringHCRipples_', thisParam.trialPhases{iTP, 1}, '_', thisParam.LFPROIs{iRegion, 1}, '_20230409'), '-dpng', '-r300');

        % save figure data
        if strcmp(thisParam.trialPhases{iTP}, 'All') && strcmp(thisParam.LFPROIs{iRegion}, 'Non-HC regions IH')
            save(strcat(paths.save, 'Fig_4b_left'), 'dt');
        elseif strcmp(thisParam.trialPhases{iTP}, 'All') && strcmp(thisParam.LFPROIs{iRegion}, 'Non-HC regions CH')
            save(strcat(paths.save, 'Fig_4b_right'), 'dt');
        end
    end
    
    % close open figures
    close all;
end

%% normalized ripple power in the different regions of interest, exactly during hippocampal ripples

% define your LFP regions of interest
thisParam           = [];
thisParam.LFPROIs   = {'AMY', {'AMYleftBipolar'; 'AMYrightBipolar'}, param.fig.colorAMY; ...
    'EC', {'ECleftBipolar'; 'ECrightBipolar'}, param.fig.colorEC; ...
    'PHC', {'PHCleftBipolar'; 'PHCrightBipolar'}, param.fig.colorPHC; ...
    'TP', {'TPleftBipolar'; 'TPrightBipolar'}, param.fig.colorTP};

% loop through the different trial phases
for iTP = 1:size(param.trial.phases, 1)
    
    % all power spectra (normalized)
    if strcmp(param.trial.phases{iTP}, 'AllPhases')
        thisSpectra = cell2mat({allRes.avLFPpowNormDuringHCRipples})'; % channel-combinations x frequencies
    else
        thisSpectra = cell2mat(cellfun(@(x) x(:, strcmp(param.trial.phases, param.trial.phases{iTP})), {allRes.avLFPpowNormDuringHCRipplesTrial}, 'UniformOutput', false))';
    end

    % create figure
    dt                  = [];
    dt.foi              = param.LFP.foi;
    dt.LFPROIs          = thisParam.LFPROIs;
    dt.allLFPChanRegion = allLFPChanRegion;
    dt.thisSpectra      = thisSpectra;
    f = LK_PlotNormalizedPowerDuringHCRipples_20231005(dt);
    
    % save figure
    LK_print(f, strcat(paths.save, 'NonHCRegions_NormalizedPowerDuringHCRipples_', param.trial.phases{iTP}, '_20230412'), '-dtiff', '-r300');
    
    % save figure data
    if strcmp(param.trial.phases{iTP}, 'AllPhases')
        save(strcat(paths.save, 'Fig_4c'), 'dt');
    end
end

%% ripple frequency in different MTL regions

% all ripple frequencies
allFreq     = cell(size(allRes, 1), 1);
for iRes = 1:size(allRes, 1)
    allFreq{iRes, 1}    = [allRes(iRes).LFPripples.freq]';
end
allMeanFreq = cellfun(@mean, allFreq);

% remove duplicates
wdMeanFreq   	= allMeanFreq(uniqueIdxIA); % "wd" stands for "without duplicates"
wdLFPChanRegion = allLFPChanRegion(uniqueIdxIA);

% statistics: ripple frequencies in the hippocampus vs. ripple frequencies
% in other regions
bHC                                 = ismember(wdLFPChanRegion, {'HCleftBipolar', 'HCrightBipolar'});
[~, pHCvsNonHC, ~, statsHCvsNonHC]  = ttest2(wdMeanFreq(bHC), wdMeanFreq(~bHC));
fprintf('Are ripple frequencies in the hippocampus different from those in extrahippocampal regions? t(%d) = %.3f, p = %.3f.\n', ...
    statsHCvsNonHC.df, statsHCvsNonHC.tstat, pHCvsNonHC);

% create shorter figure labels for the LFP channel regions
myLFPChanRegion = {...
    'HCleftBipolar', 'H_{left}'; ...
    'HCrightBipolar', 'H_{right}'; ...
    'AMYleftBipolar', 'A_{left}'; ...
    'AMYrightBipolar', 'A_{right}'; ...
    'ECleftBipolar', 'E_{left}'; ...
    'ECrightBipolar', 'E_{right}'; ...
    'PHCleftBipolar', 'P_{left}'; ...
    'PHCrightBipolar', 'P_{right}'; ...
    'TPleftBipolar', 'T_{left}'; ...
    'TPrightBipolar', 'T_{right}'; ...
    };

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 12, 6]);
axes('units', 'centimeters', 'position', [1.7, 1.6, 10, 4]);
hold on;
for iU = 1:size(myLFPChanRegion, 1)
    
    % LFP channels from this region
    bThisRegion = strcmp(wdLFPChanRegion, myLFPChanRegion{iU, 1});
    
    % mean, STE, individual data points
    m   = mean(wdMeanFreq(bThisRegion));
    ste = LK_ste(wdMeanFreq(bThisRegion));
    bar(iU, m, 'FaceColor', [1, 1, 1]);
    plot([iU, iU], [m - ste, m + ste], '-', 'Color', [0, 0, 0], 'LineWidth', 2);
    plot(iU + rand(sum(bThisRegion), 1) .* 0.5 - 0.25, wdMeanFreq(bThisRegion), '.', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 3);
end
% show significance
if pHCvsNonHC < 0.001
    plot([1, 2, 1.5, 1.5, 6.5, 6.5, 3, 10], [117, 117, 117, 119, 119, 117, 117, 117], '-', 'Color', [0, 0, 0]);
    text(4, 120, '***', 'horizontalalignment', 'center', 'fontunits', 'centimeters', 'fontsize', 0.5);
end
xl = xlabel('Region');
yl = ylabel('Ripple frequency (Hz)');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'xlim', [0.4, size(myLFPChanRegion, 1) + 0.6], 'xtick', 1:size(myLFPChanRegion, 1), 'xticklabel', myLFPChanRegion(:, 2), 'ylim', [80, 120], ...
    'xticklabelrotation', 0, 'tickdir', 'out', 'ticklength', [0.01, 0.01]);
% save figure
LK_print(f, strcat(paths.save, 'DiffRegions_MeanRippleFreqs_20230323'), '-dpng', '-r300');
