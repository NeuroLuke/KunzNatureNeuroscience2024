%==========================================================================
% This script examines the percentage of IEDs during the OpenField task and
% plots example IEDs.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param               = [];
%- subject information
param.info.name     = 'subjectdata_20210616';
%- preprocessing
param.preproc.type  = 'FBS'; % filtered, bipolar, select
%- channels
param.channel.choi  = {'HCleftBipolar'; 'HCrightBipolar'}; % channels of interest
%- artifacts
param.IED.type      = 'artifacts'; % epileptic artifacts
param.IED.fileName  = strcat(param.IED.type, '.mat');
param.IED.bPlot     = false;
%- ripples
param.ripple.folder = 'RipplesMacro_20210614';
param.ripple.specs  = '20220320_FBS_80to140Hz_tMF2';
%- for reproducibility
param.myRNG         = 7777;
rng(param.myRNG);

% paths
paths               = [];
paths.fieldtrip     = 'E:\fieldtrip\fieldtrip-20210614\';
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.ripple        = strcat('E:\OpenField\', param.ripple.folder, '\', param.ripple.specs, '\');
paths.save          = strcat('E:\OpenField\IEDsMacro_20220707\20220707_', param.ripple.folder, '_', param.ripple.specs, '_', param.IED.type, filesep);
mkdir(paths.save);

% add fieldtrip
addpath(paths.fieldtrip);
ft_defaults;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% preallocations

% save settings
save(strcat(paths.save, 'settings'));

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
    
    %% channels of interest
    
    % channels of interest in this subject
    choi        = [];
    if any(strcmp(param.channel.choi, 'HCleftBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCleftBipolar);
    end
    if any(strcmp(param.channel.choi, 'HCrightBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCrightBipolar);
    end
    
    %% loop through channels of interest
    for iChan = 1:numel(choi)
        
        %% load data
        
        % load LFP data
        data4Ripples    = load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, 'data4Ripples.mat'));
        data4Ripples    = data4Ripples.data4Ripples;
        
        % load artifact information
        IEDs            = load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, param.IED.fileName));
        IEDs            = IEDs.artifacts;
        
        % report
        fprintf('\nChannel: %s. Total number of artifacts: %d.\n', choi{iChan}, size(IEDs.artifactOnsEndsTime, 1));
        
        %% figures with example IEDs
        
        % loop through example artifacts
        if param.IED.bPlot == true
            for iA = 1:20:size(IEDs.artifactOnsEndsTime, 1)
                
                % create figure
                cfg                 = [];
                cfg.visible         = 'off';
                cfg.artifactTime    = [IEDs.artifactOnsEndsTime(iA, 1), IEDs.artifactOnsEndsTime(iA, 2)]; % time period of this artifact
                cfg.time            = data4Ripples.time{1};
                cfg.data            = data4Ripples.trial{1};
                cfg.idx             = [iSub, iChan, iA];
                f                   = LK_PlotIED_20220323(cfg);
                
                % save figure
                LK_print(f, strcat(paths.save, 'IEDExample_', subjects{iSub}, '_', choi{iChan}, '_', num2str(iA)), '-dtiff', '-r300');
                close(f);
            end
        end
        
        %% total percentage of artifacts and overall artifact rate
        
        % total duration of the experiment
        totExpDur       = range(data4Ripples.time{1});
        
        % data fraction with any sort of artifact
        totAF           = sum(isnan(data4Ripples.dataEB{1})) / numel(data4Ripples.dataEB{1});
        
        % data fraction with IEDs
        totFracIEDs     = sum(IEDs.bArtifact) / numel(IEDs.bArtifact);
        
        % duration per IED
        durationPerIED  = range(IEDs.artifactOnsEndsTime, 2);
        
        %% results
        
        % results for this channel
        thisRes                 = [];
        thisRes.subject         = subjects{iSub};
        thisRes.chan            = choi{iChan};
        thisRes.idx             = [iSub, iChan];
        
        % overall artifact characteristics
        thisRes.totExpDur       = totExpDur;
        thisRes.totAF           = totAF;
        thisRes.totFracIEDs     = totFracIEDs;
        thisRes.durationPerIED  = durationPerIED;
        
        % collect results across all sessions
        allRes              	= cat(1, allRes, thisRes);
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r   = load(strcat(paths.save, 'results.mat'));

%% total fraction of artifacts (any sort) per channel

% total fraction of artifacts for each channel
chanTotAF           = cell2mat({r.allRes.totAF}');

% report average artifact fraction
LK_ReportMeanAndSEM_20220322('\nArtifact fraction across channels', chanTotAF * 100);

%% total fraction of IEDs per channel

% total fraction of IEDs for each channel
chanTotFracIEDs     = cell2mat({r.allRes.totFracIEDs}');

% report average IED fraction
LK_ReportMeanAndSEM_20220322('\nIED fraction across channels', chanTotFracIEDs * 100);

% histogram: IED percentages
cfg                 = [];
cfg.data            = chanTotFracIEDs * 100;
cfg.binEdges        = 0:5:100;
cfg.xlabel          = 'Time with IED (%)';
cfg.ylabel          = 'Count';
cfg.title           = {'IED prevalence', ['(', LK_IndicateThousands(size(chanTotFracIEDs, 1)), ' channels)']};
cfg.fontSize        = 0.4;
cfg.figPosition     = [2, 2, 6, 6];
cfg.axPosition      = [1.75, 1.4, 3.5, 3.5];
f                   = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'chanTotFracIEDs_20220707'), '-dtiff', '-r300');
