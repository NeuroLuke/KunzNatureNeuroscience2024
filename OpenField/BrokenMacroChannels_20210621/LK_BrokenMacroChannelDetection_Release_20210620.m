%==========================================================================
% This script detects channels with large artifacts in order to flag them
% in the "subjectdata.mat" file.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;

% paths
paths           = [];
paths.data      = 'E:\OpenField\Macro_20201215\';
paths.info      = 'E:\OpenField\SubjectData_20210616\';
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';

% add fieldtrip
addpath(paths.fieldtrip);
ft_defaults;

% settings
param           = [];
param.infoName  = 'subjectdata_20210616';

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

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report and decide whether to skip this subject
    fprintf('\n=== SUBJECT: %s.\n', subjects{iSub});
    bSkip   = input('Do you want to skip this subject? (1 = yes) ');
    if bSkip == 1
        continue;
    end
    
    % load channel information
    subjectdata     = load(strcat(paths.info, subjects{iSub}, filesep, param.infoName));
    subjectdata     = subjectdata.subjectdata;
        
    % iEEG data
    eegOrig         = load(strcat(paths.data, subjects{iSub}, '\eeg_original.mat'));
    eegOrig         = eegOrig.eeg_original;
    
    % high- and low-pass filter
    cfg             = [];
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 0.5;
    cfg.hpfiltord   = 5;
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 40;
    eegPreprocessed = ft_preprocessing(cfg, eegOrig);
    
    %% loop through channels
    for iChan = 1:size(eegOrig.label, 1)
        
        % report
        fprintf('\nSubject: %s. Channel: %s.\n', subjects{iSub}, eegOrig.label{iChan});
        
        % select data from this channel using fieldtrip
        cfg                     = [];
        cfg.channel             = eegOrig.label{iChan};
        chanEEGOrig             = ft_selectdata(cfg, eegOrig);
        chanEEGPreprocessed     = ft_selectdata(cfg, eegPreprocessed);
                
        %% figure
        
        % create figure
        f = figure('units', 'centimeters', 'position', [0.5, 2, 50, 24]);
        
        % original raw data
        ax1 = axes('units', 'centimeters', 'position', [1.5, 12, 48, 10.5]);
        plot(chanEEGOrig.time{1}(1:2:end), chanEEGOrig.trial{1}(1:2:end));
        set(gca, 'xlim', [min(chanEEGOrig.time{1}), max(chanEEGOrig.time{1})], 'xticklabel', []);
        yl = ylabel('Voltage (\muV)');
        set([gca, yl], 'fontunits', 'centimeters', 'fontsize', 0.3);
        title(strrep(sprintf('Subject: %s. Channel: %s.\n', subjects{iSub}, chanEEGOrig.label{1}), '_', '\_'));
        
        % filtered data
        ax2 = axes('units', 'centimeters', 'position', [1.5, 1, 48, 10.5]);
        plot(chanEEGPreprocessed.time{1}(1:2:end), chanEEGPreprocessed.trial{1}(1:2:end));
        set(gca, 'xlim', [min(chanEEGOrig.time{1}), max(chanEEGOrig.time{1})]);
        xl = xlabel('Time (s)');
        yl = ylabel('Voltage (\muV)');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.3);
        
        % link axes and draw
        linkaxes([ax1, ax2], 'x');
        drawnow;        
        
        % decide whether this channel should be excluded
        input('Do you want to exclude this channel based on potential artifacts? ');
        close(f);        
    end
end
