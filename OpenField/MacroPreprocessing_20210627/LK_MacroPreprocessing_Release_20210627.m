%==========================================================================
% This script preprocesses the original macro-electrode iEEG data.
%
% Preprocessing includes
% - selection of usable channels;
% - conversion to bipolar montage (adjacent usable channels only);
% - line noise removal (optional).
%
% References
%   Staresina et al., Nat Neurosci, 2015
%   Norman et al., Science, 2019
%   Sakon et al., PNAS, 2022
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% paths
paths               = [];
paths.info          = 'E:\OpenField\SubjectData_20210616\';
paths.macro         = 'E:\OpenField\Macro_20201215\';
paths.save          = 'E:\OpenField\MacroPreprocessing_20210627\';
paths.pics          = strcat(paths.save, 'Pics\');
paths.functions     = 'E:\OpenField\Functions\';
paths.fieldtrip     = 'E:\fieldtrip\fieldtrip-20210614\';

% add paths
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
mkdir(paths.pics);

% settings
param                       = [];
param.info.name             = 'subjectdata_20210616';
%- for selection of channels
param.select.channel        = 'usable'; % usable | all
%- for filtering
param.filter.bsfilter       = 'yes'; % whether to use bandstop filters
param.filter.bsfilttype     = 'but'; % filter type
param.filter.bsfreq         = [50; 100; 150; 200] + [-2, 2]; % stopbands
%- for power spectrum
param.powSpectrum           = [];
param.powSpectrum.method    = 'wavelet';
param.powSpectrum.foi       = 2:4:250;
param.powSpectrum.width     = 7;
param.powSpectrum.pad       = 'nextpow2';
param.powSpectrum.output    = 'pow';
%- for plotting
param.plot.choi             = {'HAL1-HAL2'; 'HAL2-HAL3'; 'HAR1-HAR2'; 'HAR2-HAR3'; ...
    'HL1-HL2'; 'HL2-HL3'; 'HL3-HL4'; 'HL4-HL5'; 'HL5-HL6'; 'HL6-HL7'; 'HL7-HL8'; 'HL8-HL9'; 'HL9-HL10'}; % channels of interest for plotting the hippocampal data
param.plot.numPanels        = 4; % number of data segments to plot
param.plot.duration         = 2; % duration of data segments to plot

% subjects
subjects = {...
    'FR001'; ...
    'FR002a'; 'FR002b'; ...
    'FR003'; ...
    'FR004a'; 'FR004b'; ...
    'FR005a'; 'FR005b'; ...
    'FR006a'; 'FR006b'; ...
    'FR007'; ...
    'FR008a'; 'FR008b'; ...
    'FR009a'; 'FR009b'; ...
    'FR010'; ...
    'FR011'; ...
    'FR012a'; 'FR012b'; ...
    'FR013'; ...
    'FR014'; ...
    'FR015b'; ...
    'FR016'; ...
    'FR017'; ...
    'FR018'; ...
    'FR019'; ...
    'FR020'; ...
    'FR021'; ...
    'FR022'; ...
    'FR023a'; 'FR023b'; ...
    'FR024'; ...
    'FR025a'; 'FR025b'; 'FR025c'; ...
    'FR026'; ...
    'FR027'; ...
    'FR028a'; 'FR028b'; ...
    'FR029'; ...
    'FR030'; ...
    };

% save settings
save(strcat(paths.save, 'settings'));

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n==== Working on: %s.\n', subjects{iSub});
    
    % skip if processing has been done already
    subjSavePath    = strcat(paths.save, subjects{iSub}, filesep);
    if exist(fullfile(subjSavePath, 'eegFBS.mat'), 'file') > 0
        continue;
    else
        mkdir(subjSavePath);
    end
    
    % load original iEEG data
    eegO    = load(fullfile(paths.macro, subjects{iSub}, 'eeg_original.mat'));
    eegO    = eegO.eeg_original;
    
    %% usable channels
    
    % subject-information file
    subjectdata     = load(strcat(paths.info, subjects{iSub}, filesep, param.info.name));
    subjectdata     = subjectdata.subjectdata;
    
    % select usable channels
    cfg             = [];
    if strcmp(param.select.channel, 'usable')
        cfg.channel = subjectdata.macro.usable;
    elseif strcmp(param.select.channel, 'all')
        cfg.channel = subjectdata.macro.all;
    end
    eegS            = ft_selectdata(cfg, eegO); % select
    
    % select trigger channel
    cfg             = [];
    cfg.channel     = 'Trigger';
    Trigger         = ft_selectdata(cfg, eegO);
        
    %% bipolar montage
    
    % compute transition matrix to convert into bipolar montage
    tra     = LK_TransitionMatrixForBipolarMontage_20210627(eegS.label);
    
    % estimate bipolar montage using fieldtrip
    eegBS 	= ft_apply_montage(eegS, tra); % bipolar, select
        
    %% line noise removal using bandstop filter
    
    % bandstop filtering
    eegFBS  = ft_preprocessing(param.filter, eegBS); % filtered, bipolar, select
    
    %% grand average across all select channels (for control analyses) and line noise removal
    
    % compute grand average across all select channels
    eegAS           = [];
    eegAS.label     = {'GrandAverage'};
    eegAS.fsample   = eegS.fsample;
    eegAS.time      = eegS.time;
    eegAS.trial{1}  = mean(eegS.trial{1}, 1); % average, select
    
    % bandstop filtering of the grand average
    eegFAS          = ft_preprocessing(param.filter, eegAS); % filtered, average, select
        
    %% save preprocessed data
    
    % save data
    LK_save(strcat(subjSavePath, 'Trigger'), Trigger, '-v7.3'); % original trigger signal
    LK_save(strcat(subjSavePath, 'eegBS'), eegBS, '-v7.3');
    LK_save(strcat(subjSavePath, 'eegFBS'), eegFBS, '-v7.3');
    LK_save(strcat(subjSavePath, 'eegAS'), eegAS, '-v7.3');
    LK_save(strcat(subjSavePath, 'eegFAS'), eegFAS, '-v7.3');
    
    %% figure: re-referencing matrix for the bipolar montage
    
    % figure for re-referencing matrix
    f = LK_PlotTransitionMatrix(tra);
    LK_print(f, strcat(paths.pics, subjects{iSub}, '_BipolarMontage'), '-dtiff', '-r150');
    close(f);
    
    %% figures: power spectra and some hippocampal data
    
    % loop through channels of interest
    for iChan = 1:size(param.plot.choi, 1)
        
        % skip if the subject does not have this channel of interest
        if sum(strcmp(eegFBS.label, param.plot.choi{iChan})) ~= 1
            fprintf('\tNo %s channel, thus skipping.\n', param.plot.choi{iChan});
            continue;
        end
        
        %% figure: power spectrogram of the original signal and the filtered signal
        
        % time-frequency analysis to extract power
        cfg             = param.powSpectrum;
        cfg.channel     = param.plot.choi{iChan};
        cfg.toi         = min(eegFBS.time{1}):0.5:max(eegFBS.time{1});
        
        % unfiltered signal
        tfUnfiltered    = ft_freqanalysis(cfg, eegBS);
        
        % filtered signal
        tfFiltered      = ft_freqanalysis(cfg, eegFBS);
        
        %% figure: power spectra of unfiltered and filtered data
        
        % create figure
        f = figure('units', 'centimeters', 'position', [1, 1, 16, 16], 'visible', 'off');
        hold on;
        % power spectrum for unfiltered data
        thisPow     = squeeze(mean(log(tfUnfiltered.powspctrm), 3, 'omitnan'));
        plot(cfg.foi, thisPow);
        % power spectrum for filtered data
        thisPow     = squeeze(mean(log(tfFiltered.powspctrm), 3, 'omitnan'));
        plot(cfg.foi, thisPow);
        % enhance axes
        legend('Unfiltered', 'Filtered');
        xl = xlabel('f (Hz)');
        yl = ylabel('Average power over time (log(power))');
        tl = title('Power spectrum');
        set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
        % save figure
        LK_print(f, strcat(paths.pics, subjects{iSub}, '_', param.plot.choi{iChan}, '_PowerSpectra'), '-dtiff', '-r150');
        close(f);
        
        %% show some data snippets
        
        % select n time periods
        toi = transpose(round(linspace(eegFBS.time{1}(1), eegFBS.time{1}(end), param.plot.numPanels + 1)));
        toi = toi(1:param.plot.numPanels, 1);
        
        % create figure
        f = figure('units', 'centimeters', 'position', [1, 1, 16, 24], 'visible', 'off');
        
        % loop through subpanels
        for iPanel = 1:param.plot.numPanels
            
            % select filtered data
            cfg                 = [];
            cfg.channel         = param.plot.choi{iChan};
            cfg.latency         = [toi(iPanel, 1), toi(iPanel, 1) + param.plot.duration];
            eeg4FigFiltered     = ft_selectdata(cfg, eegFBS);
            
            % select unfiltered data for comparison
            eeg4FigUnfiltered   = ft_selectdata(cfg, eegBS);
        
            % plot data segment
            axes('units', 'centimeters', 'position', [1.5, (iPanel - 1) * 6 + 1.5, 14, 4]);
            hold on;
            plot(eeg4FigUnfiltered.time{1}, eeg4FigUnfiltered.trial{1}, '-', 'Color', [0.5, 0.5, 0.5]);
            plot(eeg4FigFiltered.time{1}, eeg4FigFiltered.trial{1}, 'k', 'Color', [0, 0, 0]);
            tmpAx = get(gca);
            set(gca, ...
                'xlim', [min(cfg.latency), max(cfg.latency)], 'xtick', min(cfg.latency):max(cfg.latency), ...
                'ylim', round([min(tmpAx.YLim), max(tmpAx.YLim)]), 'ytick', round([min(tmpAx.YLim), max(tmpAx.YLim)]));
            legend('Unfiltered', 'Filtered');
            xl = xlabel('Time (s)');
            yl = ylabel('Voltage (\muV)');
            set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.3);
        end
        
        % save figure
        LK_print(f, strcat(paths.pics, subjects{iSub}, '_', param.plot.choi{iChan}, '_DataSnippets'), '-dtiff', '-r150');
        close(f);
    end
end

% report
fprintf('\n\n==== Processing completed.\n');
