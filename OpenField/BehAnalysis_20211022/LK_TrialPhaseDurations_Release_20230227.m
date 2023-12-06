%==========================================================================
% This script examines the durations of the different trial phases.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% paths
paths           = [];
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.beh       = 'E:\OpenField\BehPreprocessing_20210707\';
paths.save      = 'E:\OpenField\BehAnalysis_20211022\Durations\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% settings
param               = [];
param.maxNumTrials  = 160;

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% report
fprintf('\nAnalysis of trial-phase durations.\n');
fprintf('Number of subjects: %d.\n', size(subjects, 1));

%% conditions

% specify the conditions
conditions  = {...
    'ITI', rgb('purple'); ...
    'Cue', rgb('red'); ...
    'Retrieval', rgb('orange'); ...
    'Feedback', rgb('limegreen'); ...
    'Encoding', rgb('darkgreen')};

% preallocate
allDur      = cell(size(subjects, 1), size(conditions, 1)); % duration of each trial phase
expDur      = nan(size(subjects, 1), 1); % experiment duration per session
numTrials   = nan(size(subjects, 1), 1); % number of trials per session

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nSubject: %s.\n', subjects{iSub});
    
    %% load data
    
    % load trial information
    behFile  	= dir(strcat(paths.beh, subjects{iSub, 1}, '\trialInfo.mat'));
    trialInfo   = load(fullfile(behFile.folder, behFile.name));
    trialInfo   = trialInfo.trialInfoMacrotime; % in MACROTIME
    fprintf('\tOriginal number of trials: %d.\n', size(trialInfo, 1));
    
    % remove trials beyond 160
    if size(trialInfo, 1) > param.maxNumTrials
        trialInfo   = trialInfo(1:param.maxNumTrials, :);
    end
    fprintf('\tCorrected number of trials: %d.\n', size(trialInfo, 1));
        
    %% durations of the different conditions
    
    % loop through conditions
    for iCond = 1:size(conditions, 1)
        if strcmp(conditions{iCond, 1}, 'ITI')
            allDur{iSub, iCond}     = trialInfo.Cue - trialInfo.ITI;
        elseif strcmp(conditions{iCond, 1}, 'Cue')
            allDur{iSub, iCond}     = trialInfo.Retrieval - trialInfo.Cue;
        elseif strcmp(conditions{iCond, 1}, 'Retrieval')
            allDur{iSub, iCond}     = trialInfo.Feedback - trialInfo.Retrieval;
        elseif strcmp(conditions{iCond, 1}, 'Feedback')
            allDur{iSub, iCond}     = trialInfo.Reencoding - trialInfo.Feedback;
        elseif strcmp(conditions{iCond, 1}, 'Encoding')
            allDur{iSub, iCond}     = trialInfo.Grab - trialInfo.Reencoding;
        end
    end
    
    %% experiment duration and number of trials
    
    % experiment duration: duration between first ITI onset and last grab
    % onset
    expDur(iSub, 1)     = (max(trialInfo.Grab) - min(trialInfo.ITI)) / 60; % (minutes)
    numTrials(iSub, 1)  = size(trialInfo, 1);
end

%% experiment duration

% report
fprintf('\nExperiment duration.\n');
fprintf('Range: %.3f - %.3f minutes.\n', min(expDur), max(expDur));
fprintf('Mean +/- SEM: %.3f +/- %.3f minutes.\n', mean(expDur), std(expDur) ./ sqrt(size(expDur, 1)));

% figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.6, 1.5, 4, 4]);
hold on;
plot(sort(expDur), '.', 'Color', [0, 0, 0]);
LK_yline(mean(expDur), '-', [0.75, 0.75, 0.75], 'bottom');
xl = xlabel('Session');
yl = ylabel('Duration (min)');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'xlim', [0, numel(expDur) + 1], 'ylim', [0, max(expDur)], 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
LK_print(f, strcat(paths.save, 'allSess_durations'), '-dtiff', '-r300');

%% number of trials

% number of trials per session
fprintf('\nNumber of trials.\n');
fprintf('Range: %d - %d trials.\n', min(numTrials), max(numTrials));
LK_ReportMeanAndSEM_20220322('Number of trials per session', numTrials);

%% durations of trial phases

% figure settings
myXLim = [0, 30];
myYLim = [0, 0.3];

% loop through conditions
for iCond = 1:size(conditions, 1)
    
    % select data for specific condition
    dataAll = cell2mat(allDur(:, iCond));
    
    % report
    fprintf('Durations for "%s".\n', conditions{iCond, 1});
    fprintf('Average duration: %.3f +/- %.3f, n = %.0f.\n', mean(dataAll, 'omitnan'), LK_ste(dataAll), sum(~isnan(dataAll)));
    fprintf('Minimum duration: %.3f. Maximum duration: %.3f.\n', min(dataAll), max(dataAll));
    
    % create figure
    dt              = [];
    dt.dataAll      = dataAll;
    dt.myXLim       = myXLim;
    dt.myYLim       = myYLim;
    dt.condition    = conditions{iCond, 2};
    f = LK_PlotTrialPeriodDurations_20231005(dt);

    % save figure
    LK_print(f, strcat(paths.save, 'allTrials_durations', conditions{iCond, 1}), '-dtiff', '-r300');

    % save figure data
    if strcmp(conditions{iCond}, 'Retrieval')
        save(strcat(paths.save, 'Fig_1b_leftInset'), 'dt');
    elseif strcmp(conditions{iCond}, 'Encoding')
        save(strcat(paths.save, 'Fig_1b_rightInset'), 'dt');
    end

    % figure: duration of each epoch
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
    axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 4]);
    hold on;
    plot(sort(dataAll), '.', 'Color', [0, 0, 0]);
    yline(mean(dataAll, 'omitnan'), '-', 'Color', conditions{iCond, 2}, 'LineWidth', 2);
    xl = xlabel('Trial');
    yl = ylabel('Duration (s)');
    tl = title(conditions{iCond, 1});
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    set(gca, 'xlim', [0, numel(dataAll) + 1], 'ylim', [0, max(dataAll)], 'yscale', 'log', 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    if strcmp(conditions{iCond, 1}, 'Retrieval') || strcmp(conditions{iCond}, 'Encoding')
        set(gca, 'ytick', [1, 10, 100]);
    end
    LK_print(f, strcat(paths.save, 'allTrials_durations', conditions{iCond, 1}, '_sorted'), '-dtiff', '-r300');
end
