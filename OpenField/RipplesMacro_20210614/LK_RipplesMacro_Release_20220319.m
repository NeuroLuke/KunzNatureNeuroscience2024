%==========================================================================
% This script identifies ripples in local field potentials (recorded on
% iEEG macro-electrodes).
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; clc; close all;

% settings
param                           = [];
param.myRNG                     = 777;
%- for subject information
param.info.name                 = 'subjectdata_20210616';
%- type of preprocessing
param.preproc.type              = 'FBS'; % filtered, bipolar, select
%- channels
param.channel.choi              = {'HCleftBipolar'; 'HCrightBipolar'}; % channels of interest
%- for artifact detection
param.artifact.time2Exclude     = 1; % time to exclude around artifacts (in seconds)
param.artifact.amp.idx          = 1; % amplitude
param.artifact.amp.threshFac    = 4;
param.artifact.gra.idx         	= 2; % gradient
param.artifact.gra.threshFac    = 4;
param.artifact.pow.idx          = 3; % power
param.artifact.pow.threshFac    = 4;
param.artifact.pow.lowerFreq    = 1;
param.artifact.pow.upperFreq    = 60;
param.artifact.pow.numFreqs     = 30;
param.artifact.bInspect         = false;
%- for ripple detection
param.ripple.bpfilter           = 'yes'; % whether to use a bandpass filter
param.ripple.bpfilttype         = 'but'; % filter type
param.ripple.bpfreq            	= [80, 140]; % ripple band
param.ripple.envSmoothTime      = 0.02; % time window for smoothing the envelope
param.ripple.threshMinFac       = 2; % minimum threshold
param.ripple.threshMaxFac       = Inf; % maximum threshold
param.ripple.threshPeakMinFac   = 3; % minimum peak threshold
param.ripple.durationMin        = 0.02; % minimum ripple duration
param.ripple.durationMax        = 0.5; % maximum ripple duration
param.ripple.numCyclesMin       = 3; % minimum number of cycles
param.ripple.twoi4Data          = [-3, 3];
param.ripple.twoi4BaselineERP   = [-3, 3];
param.ripple.twoi4FlipERP       = [-0.2, 0.2];
param.ripple.twoi4PlotERP       = [-0.2, 0.2];
param.ripple.twoi4BaselineTF    = [-3, 3];
param.ripple.twoi4PlotTF        = [-0.2, 0.2];
param.ripple.freq4PlotTF        = 2:4:200;
param.ripple.bInspect           = false;
%- for false-positive detection
param.ripple.falsePos.method    = 'wavelet';
param.ripple.falsePos.width     = 7; % number of wavelet cycles
param.ripple.falsePos.pad       = 'nextpow2';
param.ripple.falsePos.output    = 'pow';
param.ripple.falsePos.foi     	= 30:2:190; % frequencies for calculating the power spectrum (Hz)
param.ripple.falsePos.timeRes   = 'all';
%- for detecting surrogate ripples
param.ripple.surro.timeWindow   = 60; % (sec)

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectData_20210616\';
paths.macro     = 'E:\OpenField\MacroPreprocessing_20210627\';
paths.save      = strcat('E:\OpenField\RipplesMacro_20210614\20220320_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '\');
paths.pics      = strcat(paths.save, 'Pics\');
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';
paths.functions = 'E:\OpenField\Functions\';

% add paths
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
ft_warning('off');
mkdir(paths.pics);

% subjects
subjects    = load(strcat(paths.info, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings
save(strcat(paths.save, 'settings'));

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report and reset rng for reproducibility
    fprintf('\nSUBJECT: %s.\n\n', subjects{iSub});
    rng(param.myRNG);
    
    %% data
    
    % load subject information
    subjectdata = load(fullfile(paths.info, subjects{iSub}, param.info.name));
    subjectdata = subjectdata.subjectdata;
    
    % load preprocessed iEEG data
    if strcmp(param.preproc.type, 'FBS')
        eegAll  = load(fullfile(paths.macro, subjects{iSub}, 'eegFBS.mat')); % all electrodes
        eegAll  = eegAll.eegFBS;
        eegGA   = load(fullfile(paths.macro, subjects{iSub}, 'eegFAS.mat')); % grand average
        eegGA   = eegGA.eegFAS;
    elseif strcmp(param.preproc.type, 'BS')
        eegAll  = load(fullfile(paths.macro, subjects{iSub}, 'eegBS.mat')); % all electrodes
        eegAll  = eegAll.eegBS;
        eegGA   = load(fullfile(paths.macro, subjects{iSub}, 'eegAS.mat')); % grand average
        eegGA   = eegGA.eegAS;
    end
    
    %% artifact rejection: ripples in the grand-average signal
    
    % detect epileptic artifacts
    cfg                 = param.artifact;
    artifactsGA         = LK_ArtifactDetection_20220319(cfg, eegGA);
    
    % preprocessing of data for ripple detection
    cfg                 = param.ripple;
    cfg.bArtifact       = artifactsGA.bArtifact; % remove epileptic artifacts
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
    
    %% channels of interest
    
    % channels of interest
    choi        = [];
    if any(strcmp(param.channel.choi, 'HCleftBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCleftBipolar);
    end
    if any(strcmp(param.channel.choi, 'HCrightBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCrightBipolar);
    end
    
    %% loop through channels of interest
    for iChan = 1:numel(choi)
        
        %% check whether results already exist and if so, skip
        
        % create save path
        chanSavePath    = strcat(paths.save, subjects{iSub}, filesep, choi{iChan}, filesep);
        mkdir(chanSavePath);
        
        % continue if results already exist
        chanSaveFile    = strcat(chanSavePath, 'ripples.mat');
        if exist(chanSaveFile, 'file') > 0
            fprintf('\nRipple results already exist in %s, thus skipping.\n\n', chanSaveFile);
            continue;
        end
        
        %% select data from this channel of interest
        
        % select specific channel using fieldtrip
        cfg         = [];
        cfg.channel = choi{iChan};
        eeg         = ft_selectdata(cfg, eegAll);
        
        % report
        fprintf('\n\n------------------------------------------------------\n');
        fprintf('Working on "%s", channel "%s".\n', subjects{iSub}, eeg.label{:});
        fprintf('------------------------------------------------------\n\n');        
        
        %% artifact rejection: epileptic activity
        
        % detect epiletic artifacts in this channel
        cfg      	= param.artifact;
        artifacts   = LK_ArtifactDetection_20220319(cfg, eeg);
        
        %% visual inspection of artifacts using fieldtrip
        
        % inspect artifacts
        if param.artifact.bInspect == true
            cfg                                 = [];
            cfg.preproc.detrend                 = 'yes';
            cfg.linecolor                       = 'k';
            cfg.viewmode                        = 'vertical';
            cfg.blocksize                       = 20; % (seconds)
            cfg.artfctdef.autoartifact.artifact = artifacts.artifactOnsEndsSample; % onsets and ends of artifacts (samples)
            visualArtifacts                     = ft_databrowser(cfg, eeg);
        end
        
        %% ripple detection
        
        % preprocessing of data for ripple detection
        cfg                 = param.ripple;
        cfg.bArtifact       = artifacts.bArtifact | ripplesGA.bRipple | falsePositivesGA.bRipple; % remove epileptic artifacts, grand-average ripples, and grand-average false positives
        data4Ripples        = LK_RippleDetection_Preprocessing_20220319(cfg, eeg);
        
        % detection of ripples
        cfg                 = param.ripple;
        ripples             = LK_RippleDetection_DetectRipples_20220319(cfg, data4Ripples);
        
        % characterization of ripples
        cfg                 = [];
        cfg.idxRipples      = ripples.idxRipples;
        ripples             = LK_RippleDetection_CharacterizeRipples_20220319(cfg, data4Ripples);
                
        % data extraction for ripples
        cfg                 = param.ripple;
        cfg.ripples         = ripples.ripples;
        ripples             = LK_RippleDetection_ExtractData_20220319(cfg, data4Ripples);
        
        %% visual inspection of ripples using fieldtrip
        
        % inspect ripples
        if param.ripple.bInspect == true
            cfg                                     = [];
            cfg.preproc.bpfilter                    = 'yes';
            cfg.preproc.bpfreq                      = param.ripple.bpfreq;
            cfg.linecolor                           = 'k';
            cfg.viewmode                            = 'vertical';
            cfg.blocksize                           = 2; % (seconds)
            cfg.artfctdef.autoartifact.artifact     = [cell2mat({ripples.ripples.startIdx}'), cell2mat({ripples.ripples.endIdx}')]; % onsets and ends of artifacts (samples)
            visualRipples                           = ft_databrowser(cfg, eeg);
        end
        
        %% surrogate ripple detection
        
        % find surrogate ripples
        cfg                 = param.ripple;
        cfg.bArtifact       = artifacts.bArtifact | ripplesGA.bRipple | falsePositivesGA.bRipple; % remove epileptic artifacts, grand-average ripples, and grand-average false positives
        cfg.ripples         = ripples;
        ripplesSurro        = LK_RippleDetection_DetectSurrogates_20220319(cfg, data4Ripples);
        
        % characterization of surrogate ripples
        cfg                 = [];
        cfg.idxRipples      = ripplesSurro.idxSurrogate;
        ripplesSurro        = LK_RippleDetection_CharacterizeRipples_20220319(cfg, data4Ripples);
        
        % data extraction for surrogate ripples
        cfg                 = param.ripple;
        cfg.ripples         = ripplesSurro.ripples;
        ripplesSurro        = LK_RippleDetection_ExtractData_20220319(cfg, data4Ripples);
                
        %% channel-specific figures: time-domain amplitude average and time-frequency power average
        
        % create figure
        if isempty(ripples.eegRipples.trial)
            
            % report
            fprintf('\n==== %s. %s. No ripples detected. ====\n\n', subjects{iSub}, choi{iChan});
        else
                        
            %% figure for ripple ERP
            
            % create figure
            cfg         = [];
            cfg.figEx   = [5, 5, 6, 6];
            cfg.visible = 'off';
            cfg.param   = param;
            cfg.ripples = ripples;
            cfg.title   = {'Ripple ERP', ['(n = ', num2str(numel(ripples.eegRipples.trial)), ', ', ...
                num2str(numel(ripples.eegRipples.trial) / range(eeg.time{1}), '%.2f'), '/s)']};
            f           = LK_PlotRippleERP(cfg);
            
            % save figure
            LK_print(f, strcat(paths.pics, subjects{iSub}, '_Ripples_ERP_', ripples.eegRipples.label{1}), '-dtiff', '-r300');
            
            %% figure for ripple power spectrogram
            
            % create figure
            cfg         = [];
            cfg.figEx   = [5, 5, 6, 6];
            cfg.visible = 'off';
            cfg.param   = param;
            cfg.ripples = ripples;
            cfg.title   = {'Ripple power', ['(n = ', num2str(numel(ripples.eegRipples.trial)), ', ', ...
                num2str(numel(ripples.eegRipples.trial) / range(eeg.time{1}), '%.2f'), '/s)']};
            f           = LK_PlotRipplePowerSpectrum(cfg);
            
            % save figure
            LK_print(f, strcat(paths.pics, subjects{iSub}, '_Ripples_TFPower_', ripples.eegRipples.label{1}), '-dtiff', '-r300');
        end
        
        %% save data

        % save ripple and artifact data
        save(strcat(chanSavePath, 'ripples'), 'ripples');
        save(strcat(chanSavePath, 'data4Ripples'), 'data4Ripples');
        save(strcat(chanSavePath, 'ripplesGA'), 'ripplesGA');
        save(strcat(chanSavePath, 'artifacts'), 'artifacts');
        save(strcat(chanSavePath, 'ripplesSurro'), 'ripplesSurro');
        
        % close all open figures
        close all;
    end
end

%% concatenate the data from all ripples or from all channels - mean ERPs

% preallocate
allMeanEEGRipples       = []; % all average ripple ERPs (averaged across ripples from the same channel)
allMeanEEGRipplesSurro  = [];

% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % path with previously saved data
    sessSavePath    = strcat(paths.save, subjects{iSub}, filesep);
    
    % channels for this subject
    choi            = LK_dir(sessSavePath);
    
    % loop through channels
    for iChan = 1:size(choi, 1)
        
        %% load ripples from this session and channel
        
        % load ripple results from this channel
        ripples                 = load(fullfile(choi(iChan).folder, choi(iChan).name, 'ripples.mat'));
        eegRipples              = ripples.ripples.eegRipples;
        eegRipples.label        = {'choi'}; % rename to channel of interest (choi)
        
        % load surrogate ripple results from this channel
        ripplesSurro            = load(fullfile(choi(iChan).folder, choi(iChan).name, 'ripplesSurro.mat'));
        eegRipplesSurro         = ripplesSurro.ripplesSurro.eegRipples;
        eegRipplesSurro.label   = {'choi'};
        
        %% baseline-corrected average ripple from this channel
        
        % baseline correction
        cfg                 = [];
        cfg.baseline        = param.ripple.twoi4BaselineERP; % e.g., [-3, 3]
        meanEEGRipples      = ft_timelockbaseline(cfg, eegRipples);
        meanEEGRipplesSurro = ft_timelockbaseline(cfg, eegRipplesSurro);
        
        % average
        cfg                 = [];
        meanEEGRipples      = ft_timelockanalysis(cfg, meanEEGRipples);
        meanEEGRipplesSurro = ft_timelockanalysis(cfg, meanEEGRipplesSurro);
        
        % update overall data
        if isempty(allMeanEEGRipples)
            
            % initialize average ripples and average surrogate ripples
            allMeanEEGRipples       = meanEEGRipples;
            allMeanEEGRipplesSurro  = meanEEGRipplesSurro;
        else
            
            % append average ripples and average surrogate ripples
            cfg                     = [];
            allMeanEEGRipples       = ft_appenddata(cfg, allMeanEEGRipples, meanEEGRipples);
            allMeanEEGRipplesSurro  = ft_appenddata(cfg, allMeanEEGRipplesSurro, meanEEGRipplesSurro);
        end
    end
end

%% grand-average ERP for ripples and surrogate ripples - across channels

% loop through real and surrogate ripples
groups  = {'RipplesSurro', [5, 5, 6, 6], [1.7, 1.4, 3.5, 3.5], [-0.5, 0.5], [-10, 25], [-0.04, 0.04], 0.4; ... % type, fig size, axis size, x-limits, y-limits, zoom window, font size
    'RipplesSurro', [5, 5, 4.2, 4], [1.25, 1.2, 2.25, 2.25], [-0.04, 0.04], [-10, 25], [], 0.6; ... % zoom-in
    'Ripples', [5, 5, 6, 6], [1.7, 1.4, 3.5, 3.5], [-0.5, 0.5], [-10, 25], [-0.04, 0.04], 0.4; ...
    'Ripples', [5, 5, 4.2, 4], [1.2, 1.25, 2.25, 2.25], [-0.04, 0.04], [-10, 25], [], 0.6; ... % zoom-in
    'Ripples', [5, 5, 6, 6], [1.7, 1.4, 3.5, 3.5], [-1, 1], [-10, 25], [], 0.4; ... % zoom-out
    'Ripples_MeanOnly', [5, 5, 6, 6], [1.7, 1.4, 3.5, 3.5], [-0.25, 0.25], [-10, 25], [], 0.4}; % for illustrating the coactivity analysis
for iGroup = 1:size(groups, 1)
    
    % select data for this group
    if strcmp(groups{iGroup, 1}, 'Ripples') || strcmp(groups{iGroup, 1}, 'Ripples_MeanOnly')
        selData = allMeanEEGRipples;
    elseif strcmp(groups{iGroup, 1}, 'RipplesSurro')
        selData = allMeanEEGRipplesSurro;
    end
    
    % time and data
    time    = mean(cell2mat(selData.time'), 1);
    data    = cell2mat(selData.trial');
        
    % mean and standard error
    m       = mean(data, 1);
    ste     = LK_ste(data);
    
    %% figure for ERP

    % create figure
    dt              = [];
    dt.groups       = groups;
    dt.iGroup       = iGroup;
    dt.time         = time;
    dt.m            = m;
    dt.ste          = ste;
    dt.numChannels  = size(data, 1);
    f = LK_PlotRippleERPAcrossChannels_20231030(dt);
    
    % save figure
    LK_print(f, strcat(paths.pics, groups{iGroup}, '_ERP_averageAcrossChannels_', num2str(range(groups{iGroup, 4}) * 1000), 'ms_20230320'), '-dtiff', '-r300');

    % save figure data
    if strcmp(groups{iGroup, 1}, 'Ripples') && range(groups{iGroup, 4}) == 1
        save(strcat(paths.pics, 'Fig_2h_left'), 'dt');
    elseif strcmp(groups{iGroup, 1}, 'Ripples') && range(groups{iGroup, 4}) == 0.08
        save(strcat(paths.pics, 'Fig_2h_left_inset'), 'dt');
    end
end

%% concatenate the data from all ripples or from all channels - mean power spectrograms

% preallocate
allMeanPowRipples       = []; % all average ripple power spectrograms (power average across ripples from the same channel)
allMeanPowRipplesSurro  = [];

% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % path with previously saved data
    sessSavePath    = strcat(paths.save, subjects{iSub}, filesep);

    % channels
    choi            = LK_dir(sessSavePath);
    
    % loop through channels
    for iChan = 1:size(choi, 1)
        
        %% load ripples from this subject
        
        % load ripple results from this channel
        ripples                 = load(fullfile(choi(iChan).folder, choi(iChan).name, 'ripples.mat'));
        eegRipples              = ripples.ripples.eegRipples;
        eegRipples.label        = {'choi'}; % rename to channel of interest (choi)
        
        % load surrogate ripple results from this channel
        ripplesSurro            = load(fullfile(choi(iChan).folder, choi(iChan).name, 'ripplesSurro.mat'));
        eegRipplesSurro         = ripplesSurro.ripplesSurro.eegRipples;
        eegRipplesSurro.label   = {'choi'};
                
        %% average power spectrum for ripples from this channel
        
        % ripple TF power
        cfg                     = [];
        cfg.trials              = 'all'; % trials to use
        cfg.method              = 'wavelet';
        cfg.width               = 7;
        cfg.output              = 'pow';
        cfg.foi                 = param.ripple.freq4PlotTF;
        cfg.toi                 = 'all';
        cfg.keeptrials          = 'no';
        cfg.pad                 = 'nextpow2';
        meanPowRipples          = ft_freqanalysis(cfg, eegRipples);
        meanPowRipplesSurro     = ft_freqanalysis(cfg, eegRipplesSurro);
        
        % transformation of power values to relative change
        cfg                     = [];
        cfg.baseline            = param.ripple.twoi4BaselineTF;
        cfg.baselinetype        = 'relchange';
        meanPowRipples          = ft_freqbaseline(cfg, meanPowRipples);
        meanPowRipplesSurro     = ft_freqbaseline(cfg, meanPowRipplesSurro);
        
        % update overall data
        if isempty(allMeanPowRipples)
            
            % initialize average ripple power and average surrogate ripple
            % power
            allMeanPowRipples       = meanPowRipples;
            allMeanPowRipplesSurro  = meanPowRipplesSurro;
        else
            
            % append average ripple power and average surrogate ripple
            % power
            cfg                     = [];
            allMeanPowRipples       = ft_appenddata(cfg, allMeanPowRipples, meanPowRipples);
            allMeanPowRipplesSurro  = ft_appenddata(cfg, allMeanPowRipplesSurro, meanPowRipplesSurro);
        end
    end
end

%% grand-average power spectrogram for ripples and surrogate ripples - across channels

% loop through different data groups
groups  = {'RipplesSurro', [5, 5, 6, 6], [1.5, 1.4, 3.5, 3.5], [-0.1, 0.1], [0, 5]; ... % ripple type, figure size, axis size, x-limits, z-limits
    'Ripples', [5, 5, 6, 6], [1.5, 1.4, 3.5, 3.5], [-0.1, 0.1], [0, 5]; ...
    'ActivityAroundRipples', [5, 5, 6.5, 6], [1.5, 1.4, 3.5, 3.5], [-1, 1], [-0.5, 0.5]};
for iGroup = 1:size(groups, 1)
    
    % select data for this group
    if strcmp(groups{iGroup, 1}, 'Ripples') || strcmp(groups{iGroup, 1}, 'ActivityAroundRipples')
        selData = allMeanPowRipples;
    elseif strcmp(groups{iGroup, 1}, 'RipplesSurro')
        selData = allMeanPowRipplesSurro;
    end
    
    % average across channels
    cfg     = [];
    data    = ft_timelockanalysis(cfg, selData);
    
    % create figure showing power spectrogram
    dt              = [];
    dt.iGroup       = iGroup;
    dt.groups       = groups;
    dt.data         = data;
    dt.param        = param;
    dt.numChannels  = numel(selData.trial);
    f = LK_PlotRipplePowerSpectrumAcrossChannels_20231030(dt);
    
    % save figure
    LK_print(f, strcat(paths.pics, groups{iGroup}, '_TFPower_averageAcrossChannels', '_', num2str(range(groups{iGroup, 4}) * 1000), 'ms_20230320'), '-dtiff', '-r300');

    % save figure data
    if strcmp(groups{iGroup}, 'Ripples')
        save(strcat(paths.pics, 'Fig_2h_right'), 'dt');
    end
end

%% ripple characteristics

% report
fprintf('\nExtracting ripple characteristics.\n');

% preallocate
chanRippleTimepoints        = []; % ripple timepoints for all ripples
chanRippleDurations         = []; % ripple durations for all ripples
chanInterRippleIntervals    = []; % inter-ripple-intervals for all ripples
chanRippleFreqs             = []; % ripple frequencies for all ripples
sessRippleRate              = cell(size(subjects, 1), 1); % ripple rate per channel
sessRippleRateArtFree       = cell(size(subjects, 1), 1); % ripple rate per channel in artifact-free periods
sessArtFrac                 = cell(size(subjects, 1), 1); % artifact fraction per channel

% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % ripple channels
    chans                   = LK_dir(strcat(paths.save, subjects{iSub}));
    
    % loop through channels
    thisArtFrac             = nan(size(chans, 1), 1);
    thisRippleRate          = nan(size(chans, 1), 1);
    thisRippleRateArtFree   = nan(size(chans, 1), 1);
    for iChan = 1:size(chans, 1)
        
        %% artifacts
        
        % load artifact information
        artifacts                   = load(fullfile(chans(iChan).folder, chans(iChan).name, 'artifacts.mat'));
        artifacts                   = artifacts.artifacts;
        
        % fraction of artifacts in this channel
        thisArtFrac(iChan, 1)       = sum(artifacts.bArtifact) / numel(artifacts.bArtifact);
        
        %% ripples
        
        % load ripple time points
        ripples                     = load(fullfile(chans(iChan).folder, chans(iChan).name, 'ripples.mat'));
        ripples                     = ripples.ripples;
                
        % ripple-peak timepoints, durations, and inter-ripple-intervals
        thisRippleTimepoints        = cell2mat({ripples.ripples.peakTime}'); % timepoints are in seconds
        thisRippleDurations         = [ripples.ripples.endTime]' - [ripples.ripples.startTime]'; % ripple durations
        thisIRI                     = diff(thisRippleTimepoints); % inter-ripple intervals
        thisRippleFreqs             = cell2mat({ripples.ripples.freq}');
        
        % collect across channels and subjects
        chanRippleTimepoints       	= cat(1, chanRippleTimepoints, {thisRippleTimepoints});
        chanRippleDurations         = cat(1, chanRippleDurations, {thisRippleDurations});
        chanInterRippleIntervals   	= cat(1, chanInterRippleIntervals, {thisIRI});
        chanRippleFreqs             = cat(1, chanRippleFreqs, {thisRippleFreqs});
        
        % overall ripple rate
        thisDuration                = numel(ripples.bRipple) / ripples.eegRipples.fsample; % seconds
        thisRippleRate(iChan, 1)    = numel(thisRippleTimepoints) / thisDuration;
        
        % overall ripple rate, artifact-free periods
        thisDuration                    = sum(~artifacts.bArtifact) / ripples.eegRipples.fsample; % seconds
        thisRippleRateArtFree(iChan, 1) = numel(thisRippleTimepoints) / thisDuration;
    end
    
    % collect across subjects
    sessArtFrac{iSub, 1}            = thisArtFrac;
    sessRippleRate{iSub, 1}         = thisRippleRate;
    sessRippleRateArtFree{iSub, 1}  = thisRippleRateArtFree;
end

%% figure with histogram for ripple rate per channel

% unfold across channels: ripple rates and artifact fractions
chanRippleRate          = cell2mat(sessRippleRate);
chanRippleRateArtFree   = cell2mat(sessRippleRateArtFree);

% report
LK_ReportMeanAndSEM_20220322('Ripple rate per channel', chanRippleRate);
LK_ReportMeanAndSEM_20220322('Ripple rate per channel, only artifact-free periods', chanRippleRateArtFree);

% histogram: ripple rates
cfg                 = [];
cfg.data            = chanRippleRate;
cfg.binEdges        = 0:0.01:0.3;
cfg.xlabel          = 'Ripple rate (Hz)';
cfg.ylabel          = 'Count';
cfg.title           = {'Ripple rate', ['(', LK_IndicateThousands(size(chanRippleRate, 1)), ' channels)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [1.75, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.pics, 'chanRippleRate_20220321'), '-dtiff', '-r300');

%% ripple rates during first vs. second sessions, unpaired test across channels

% identify first sessions vs. second/third sessions
bNotFirst   = cellfun(@(x) contains(x(end), {'b', 'c'}), subjects(:, 1));
bFirst      = ~bNotFirst;

% compare ripple rates between first vs. subsequent sessions
[~, p, ~, stats]    = ttest2(cell2mat(sessRippleRate(bFirst)), cell2mat(sessRippleRate(bNotFirst)));
fprintf('\nTesting whether ripple rates change between subsequent sessions:\n');
fprintf('Number of channels in first and second sessions, respectively: %d and %d.\n', ...
    size(cell2mat(sessRippleRate(bFirst)), 1), size(cell2mat(sessRippleRate(bNotFirst)), 1));
fprintf('Are the ripple rates higher in first than in later sessions (unpaired test)? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% create figure
x   = [1, 2];
d   = {cell2mat(sessRippleRate(bFirst)); cell2mat(sessRippleRate(bNotFirst))};
n   = {sum(~isnan(d{1})), sum(~isnan(d{2}))};
m   = {mean(d{1}), mean(d{2})};
ste = {LK_ste(d{1}), LK_ste(d{2})};
f = figure('units', 'centimeters', 'position', [2, 2, 5, 6]);
axes('units', 'centimeters', 'position', [1.5, 1.4, 3, 3.5]);
hold on;
for iM = 1:numel(m)
    bar(x(iM), m{iM}, 'EdgeColor', [0, 0, 0], 'FaceColor', [0.9, 0.9, 0.9]);
    plot([x(iM), x(iM)], [m{iM} - ste{iM}, m{iM} + ste{iM}], 'Color', [0, 0, 0], 'LineWidth', 2);
    plot(x(iM) + rand(size(d{iM}, 1), 1) .* 0.5 - 0.25, d{iM}, '.', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 3);
end
xl = xlabel('Sessions');
yl = ylabel('Ripple rate (Hz)');
tl = title({'Ripple rate', sprintf('(%d vs %d channels)', n{1}, n{2})});
set(gca, 'xlim', [min(x) - 0.8, max(x) + 0.8], 'xtick', x, 'xticklabel', {'First', 'Later'}, 'ylim', [0, 0.3], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
LK_print(f, strcat(paths.pics, 'chanRippleRate_ComparisonBetweenSessions_20220705'), '-dtiff', '-r300');

%% figure for inter-ripple intervals

% unfold across ripples
ripplesInterRippleIntervals = cell2mat(chanInterRippleIntervals);

% report inter-ripple intervals
fprintf('\nNumber of ripples: %d.\n', size(cell2mat(chanRippleTimepoints), 1));
LK_ReportMeanAndSEM_20220322('\nRipple inter-ripple interval (s)', ripplesInterRippleIntervals);
fprintf('Number of inter-ripple intervals <15 ms: %d.\n', sum(ripplesInterRippleIntervals < 0.015));

% histogram: inter-ripple intervals
cfg                 = [];
cfg.data            = ripplesInterRippleIntervals;
cfg.binEdges        = 0:0.5:10;
cfg.xlabel          = 'Inter-ripple interval (s)';
cfg.ylabel          = 'Count';
cfg.title           = {'Inter-ripple intervals', ['(', LK_IndicateThousands(size(ripplesInterRippleIntervals, 1)), ' intervals)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [1.75, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.pics, 'ripplesInterRippleIntervals_20220321'), '-dtiff', '-r300');

%% figure for ripple durations

% unfold across ripples
ripplesRippleDurations  = cell2mat(chanRippleDurations);

% report grand-average mean
fprintf('\nNumber of ripples: %d.\n', size(ripplesRippleDurations, 1));
LK_ReportMeanAndSEM_20220322('\nRipple duration (ms)', ripplesRippleDurations * 1000);
fprintf('Number of ripple durations <= 20 ms: %d.\n', sum(ripplesRippleDurations <= 0.020));

% histogram: ripple durations
cfg                 = [];
cfg.data            = ripplesRippleDurations .* 1000; % convert to milliseconds
cfg.binEdges        = 0:5:150;
cfg.xlabel          = 'Ripple duration (ms)';
cfg.ylabel          = 'Count';
cfg.title           = {'Ripple durations', ['(', LK_IndicateThousands(size(ripplesRippleDurations, 1)), ' ripples)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [1.75, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.pics, 'ripplesRippleDurations_20220321'), '-dtiff', '-r300');

%% figure for ripple frequencies

% average ripple frequency per channel
ripplesRippleFreqs  = cell2mat(chanRippleFreqs);

% report grand-average mean
LK_ReportMeanAndSEM_20220322('\nRipple frequency (Hz)', ripplesRippleFreqs);

% histogram: ripple frequencies
cfg                 = [];
cfg.data            = ripplesRippleFreqs;
cfg.binEdges        = 80:2:140;
cfg.xlabel          = 'Ripple frequency (Hz)';
cfg.ylabel          = 'Count';
cfg.title           = {'Ripple frequencies', ['(', LK_IndicateThousands(size(ripplesRippleFreqs, 1)), ' ripples)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [2.05, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.pics, 'chanMeanRippleFreqs_20220321'), '-dtiff', '-r300');

%% figure for ripple frequencies, separately for each channel

% sort by median ripple frequency
[~, I]                  = sort(cellfun(@median, chanRippleFreqs));
sortedChanRippleFreqs   = chanRippleFreqs(I);

% assign channel index to each ripple frequency
sortedChanIdx           = cell(size(sortedChanRippleFreqs));
for iChan = 1:size(sortedChanRippleFreqs, 1)
    sortedChanIdx{iChan, 1}   = iChan * ones(size(sortedChanRippleFreqs{iChan}, 1), 1);
end

% statistics: one-way ANOVA to test whether ripple frequencies differ
% between channels
[pRF, tblRF, statsRF]   = anova1(cell2mat(sortedChanRippleFreqs), cell2mat(sortedChanIdx), 'off');
cRF                     = multcompare(statsRF);

% create boxplot figure for ripple frequencies per channel
f = figure('units', 'centimeters', 'position', [2, 2, 18, 6]);
axes('units', 'centimeters', 'position', [1.7, 1.1, 16, 4.5]);
hold on;
for iChan = 1:size(sortedChanRippleFreqs, 1)
    % minimum and maximum ripple frequency for this channel
    plot([iChan, iChan], [min(sortedChanRippleFreqs{iChan}), max(sortedChanRippleFreqs{iChan})], '.', 'Color', [0.6, 0.6, 0.6]);
    % 25th and 75th percentile
    ci = prctile(sortedChanRippleFreqs{iChan}, [25, 75]);
    plot([iChan, iChan], [ci(1), ci(2)], '-', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1);
    % median
    plot(iChan, median(sortedChanRippleFreqs{iChan}), 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
end
xl = xlabel('Channel', 'units', 'normalized', 'position', [0.5, -0.1]);
yl = ylabel('Ripple frequency (Hz)');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
set(gca, 'xlim', [0, numel(sortedChanRippleFreqs) + 1], 'xtick', [1, numel(sortedChanRippleFreqs)], 'ylim', [80, 140], 'box', 'off', 'tickdir', 'out', 'ticklength', [0.01, 0.01]);
LK_print(f, strcat(paths.pics, 'chanRippleFreqs_20230317'), '-dpng', '-r300');

%% figure for relationships between channel-wise ripple rate, ripple duration, and ripple frequency

% channel-wise ripple duration and frequency
chanRippleDuration  = cellfun(@mean, chanRippleDurations);
chanRippleFreq      = cellfun(@mean, chanRippleFreqs);

% comparisons
comparisons     = {'chanRippleRate2ChanRippleDuration', 'poly1', 'Ripple rate (Hz)', 'Ripple duration (ms)'; ...
    'chanRippleRate2ChanRippleFreq', 'power1', 'Ripple rate (Hz)', 'Ripple frequency (Hz)'; ...
    'chanRippleDuration2ChanRippleFreq', 'power1', 'Ripple duration (ms)', 'Ripple frequency (Hz)'};

% loop through comparisons
for iComp = 1:size(comparisons, 1)
    
    % select data
    if strcmp(comparisons{iComp, 1}, 'chanRippleRate2ChanRippleDuration')
        data1   = chanRippleRate;
        data2   = chanRippleDuration .* 1000; % convert to ms
    elseif strcmp(comparisons{iComp, 1}, 'chanRippleRate2ChanRippleFreq')
        data1   = chanRippleRate;
        data2   = chanRippleFreq;
    elseif strcmp(comparisons{iComp, 1}, 'chanRippleDuration2ChanRippleFreq')
        data1   = chanRippleDuration * 1000; % convert to ms
        data2   = chanRippleFreq;
    end
    
    % correlation between channel-wise ripple rate and ripple frequency
    [rho, pval] = corr(data1, data2, 'type', 'spearman');
    fprintf('Correlation for the comparison "%s": Spearman''s rho = %.3f, p = %.3f.\n', comparisons{iComp, 1}, rho, pval);
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
    axes('units', 'centimeters', 'position', [1.6, 1.5, 4, 4]);
    hold on;
    plot(data1, data2, '.', 'Color', [0, 0, 0], 'MarkerSize', 12);
    % fit curve and plot fitted curve
    [fitobject, gof] = fit(data1, data2, comparisons{iComp, 2});
    tmpAx = get(gca);
    x = min(tmpAx.XLim):0.001:max(tmpAx.XLim);
    plot(x, fitobject(x), '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1); % plot fitted curve
    % enhance axes
    xl = xlabel(comparisons{iComp, 3});
    yl = ylabel(comparisons{iComp, 4});
    t1 = text(0.975, 1, sprintf('\\itR\\rm^{2} = %.3f', gof.rsquare), 'units', 'normalized', 'horizontalalignment', 'right', 'verticalalignment', 'top');
    set([gca, xl, yl, t1], 'fontunits', 'centimeters', 'fontsize', 0.4);
    set(gca, 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    LK_print(f, strcat(paths.pics, comparisons{iComp, 1}), '-dpng', '-r300');
end

%% figure for artifact fraction per channel

% report average artifact fraction
chanArtFrac         = cell2mat(sessArtFrac);
LK_ReportMeanAndSEM_20220322('\nArtifact fraction across channels', chanArtFrac * 100);

% histogram: artifact percentages
cfg                 = [];
cfg.data            = chanArtFrac * 100;
cfg.binEdges        = 0:5:100;
cfg.xlabel          = 'Time with IED (%)';
cfg.ylabel          = 'Count';
cfg.title           = {'Artifact fraction', ['(', LK_IndicateThousands(size(chanArtFrac, 1)), ' channels)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [1.75, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.pics, 'chanArtFrac_20220321'), '-dtiff', '-r300');

% correlation between ripple rates and artifact fractions
[rho, pval]         = corr(chanRippleRate, chanArtFrac, 'type', 'spearman');
fprintf('Correlation between ripple rate and artifact fraction across channels: rho = %.3f, p = %.3f.\n', rho, pval);

% correlation between ripple rates in artifact-free periods and artifact
% fractions
[rho, pval]         = corr(chanRippleRateArtFree, chanArtFrac, 'type', 'spearman');
fprintf('Correlation between ripple rate (in artifact-free periods) and artifact fraction across channels: rho = %.3f, p = %.3f.\n', rho, pval);
