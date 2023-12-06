%==========================================================================
% This script examines the memory performance of the subjects.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.save          = 'E:\OpenField\BehAnalysis_20211022\MemoryPerformance\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% subjects
subjects            = load(strcat(paths.subjects, 'subjects.mat'));
subjects            = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% settings
param                   = [];
param.timeEdges         = 0:0.05:1; % normalized time edges
param.timeCenters       = movmean(param.timeEdges, 2, 'endpoints', 'discard'); % normalized time centers
param.maxNumTrials      = 160; % maximum number of trials
param.numTrialsPerBlock = 20; % number of trials per block
param.maxNumSess        = 3; % maximum number of sessions that a subject performed
param.myRNG             = 777;
rng(param.myRNG); % set random seed

% possible random locations within arena
dt                  = [];
dt.maxR             = 5000; % maximum radius
dt.minR             = 0; % minimum radius
dt.N                = 10000001; % number of locations to create
dt.centerX          = 0; % arena center x
dt.centerY          = 0; % arena center y
locsInArena         = LK_RandomPointsInCircle(dt);

% preallocate
allMemPerf          = cell(size(subjects, 1), 1);
allNormTime         = cell(size(subjects, 1), 1);
allMemPerfNormTime  = nan(size(subjects, 1), numel(param.timeCenters));
allMemPerfEarlyLate = nan(size(subjects, 1), 2); % early, late
allObjects          = cell(size(subjects, 1), 1);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nSubject: %s.\n', subjects{iSub});
    
    %% load data
    
    % load trial information
    behFile  	= dir(strcat(paths.beh, subjects{iSub, 1}, '\trialInfo.mat'));
    trials      = load(fullfile(behFile.folder, behFile.name));
    trials    	= trials.trialInfoMacrotime;
    fprintf('\tOriginal number of trials: %d.\n', size(trials, 1));
    
    % restrict to 160 trials (after trial #160, the session is over)
    if size(trials, 1) > param.maxNumTrials
        trials  = trials(1:param.maxNumTrials, :);
    end
    fprintf('\tCorrected number of trials: %d.\n', size(trials, 1));
    fprintf('\tNumber of trials with a response: %d.\n', sum(~isnan(trials.DropError)));
    
    %% drop error and memory performance
    
    % drop error and memory performance per trial
    dropError       = sqrt((trials.xCorrect - trials.xResponse) .^ 2 + (trials.yCorrect - trials.yResponse) .^ 2); % drop error
    dropErrorSurro  = pdist2([trials.xCorrect, trials.yCorrect], locsInArena); % surrogate drop error
    memPerf         = sum(dropError < dropErrorSurro, 2) ./ sum(~isnan(dropErrorSurro), 2); % memory performance
    
    % normalized time (use feedback timepoints as these are the timepoints
    % when the subjects give their response)
    normTime            = (trials.Feedback - min(trials.Feedback)) ./ range(trials.Feedback);
    
    % average memory performance for each normalized time window
    memPerfNormTime     = nan(numel(param.timeCenters), 1);
    for iTime = 1:numel(param.timeCenters)
        if iTime < numel(param.timeCenters)
            logIdx                  = normTime >= param.timeEdges(iTime) & normTime < param.timeEdges(iTime + 1);
        else
            logIdx                  = normTime >= param.timeEdges(iTime) & normTime <= param.timeEdges(iTime + 1); % last bin
        end
        memPerfNormTime(iTime, 1)   = mean(memPerf(logIdx), 'omitnan');
    end
    
    % average memory performance early vs. late during the task
    bHalf1                          = (1:numel(memPerf)) <= (numel(memPerf) / 2);
    memPerfEarlyLate                = [mean(memPerf(bHalf1), 'omitnan'); mean(memPerf(~bHalf1), 'omitnan')];
    
    % collect across subjects
    allMemPerf{iSub, 1}             = memPerf; % memory performance per trial
    allNormTime{iSub, 1}            = normTime; % normalized time per trial
    allMemPerfNormTime(iSub, :)     = memPerfNormTime; % memory performance per normalized time
    allMemPerfEarlyLate(iSub, :)    = memPerfEarlyLate; % memory performance during the first and second half of all trials
    allObjects{iSub, 1}             = trials.Object; % object names
end

%% memory performance during the first vs. second half of the data (separately for all, first, and later sessions)

% session groups
groups          = {'allSessions'; 'firstSessions'; 'laterSessions'};

% identify first sessions
bFirstSession   = ~cellfun(@(x) contains(x(end), {'b', 'c'}), subjects(:, 1));

% loop through different groups
for iGroup = 1:size(groups, 1)
    
    % select sessions based on group
    if strcmp(groups{iGroup}, 'allSessions')
        bMaskSession = true(size(bFirstSession));
    elseif strcmp(groups{iGroup}, 'firstSessions')
        bMaskSession = bFirstSession;
    elseif strcmp(groups{iGroup}, 'laterSessions')
        bMaskSession = ~bFirstSession;
    end
    
    % report
    fprintf('\nAnalysis of memory performance in: %s (n = %d).\n', groups{iGroup}, sum(bMaskSession));
    
    % compare memory accuracy in the first vs. second data half
    [~, pHalves, ~, statsHalves]    = ttest(allMemPerfEarlyLate(bMaskSession, 1), allMemPerfEarlyLate(bMaskSession, 2));
    fprintf('Testing memory accuracy in the first vs. second data half: t(%d) = %.3f, p = %.3f.\n', statsHalves.df, statsHalves.tstat, pHalves);
    
    % compare memory accuracy in the first vs. second data half,
    % microwire-subjects only
    bMicrowire                      = strcmp(subjects(:, 2), 'Microwire');
    [~, p, ~, stats]                = ttest(allMemPerfEarlyLate(bMaskSession & bMicrowire, 1), allMemPerfEarlyLate(bMaskSession & bMicrowire, 2));
    fprintf('Testing memory accuracy in the first vs. second data half, microwire subjects only: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
    
    % figure
    dt                      = [];
    dt.subjects             = subjects;
    dt.subjects(:, 1)       = {''}; % anonymize
    dt.bMaskSession         = bMaskSession;
    dt.allMemPerfEarlyLate  = allMemPerfEarlyLate;
    dt.pHalves              = pHalves;
    f = LK_PlotMemoryPerformanceHalves_20231004(dt);

    % save figure
    print(f, strcat(paths.save, 'allMemPerfEarlyLate_', groups{iGroup}, '_20230227'), '-dtiff', '-r450');

    % save figure data
    if strcmp(groups{iGroup}, 'allSessions')
        save(strcat(paths.save, 'Fig_1d_left'), 'dt');
    end
end

%% stats: memory performance over normalized time

% mean and standard error across sessions
m   = mean(allMemPerfNormTime, 1, 'omitnan');
ste = LK_ste(allMemPerfNormTime);

% correlation between normalized time and memory performance for the mean
[rhoMean, pvalMean] = corr(param.timeCenters', m');
fprintf('\nCorrelation between normalized time and mean memory performance: rho = %.3f, p = %.3f.\n', rhoMean, pvalMean);

% correlation between normalized time and memory performance, session-wise
allRho  = nan(size(subjects, 1), 1);
for iSub = 1:size(subjects, 1)
    allRho(iSub, 1) = corr(param.timeCenters', allMemPerfNormTime(iSub, :)', 'rows', 'complete');
end

% t-test of the session-wise correlation values against 0
[~, p, ~, stats]    = ttest(allRho);
fprintf('T-test of the time-memory correlation values against 0: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

%% figure: memory performance over normalized time

% figure showing the development of memory performance over normalized time
dt          = [];
dt.param    = param;
dt.m        = m;
dt.ste      = ste;
dt.pvalMean = pvalMean;
f = LK_PlotMemoryPerformanceNormTime_20231004(dt);

% save figure
LK_print(f, strcat(paths.save, 'allMemPerfNormTime_20230227'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_1d_right'), 'dt');

%% figure: memory performance per blocks of trials

% loop through sessions
blockMemPerf    = nan(size(allMemPerf, 1), ceil(max(cellfun(@numel, allMemPerf)) / param.numTrialsPerBlock));
for iSub = 1:size(allMemPerf, 1)
    
    % memory performances in this session
    thisMemPerf = allMemPerf{iSub};
    blockIdx    = ceil(transpose(1:numel(thisMemPerf)) / param.numTrialsPerBlock);
    
    % loop through blocks
    for iBlock = 1:max(blockIdx)
        blockMemPerf(iSub, iBlock)  = mean(thisMemPerf(blockIdx == iBlock), 'omitnan');
    end
end

% test for significant differences in memory performance per block between
% first and later sessions
[~, pBlockMemFirstVSLater, ~, statsBlockMemFirstVSLater] = ttest2(blockMemPerf(bFirstSession, :), blockMemPerf(~bFirstSession, :));
pBlockMemFirstVSLater = pBlockMemFirstVSLater * numel(pBlockMemFirstVSLater); % Bonferroni correction for multiple comparisons
fprintf('\nBlock-wise comparison of memory performance between first and later sessions:\n');
disp(statsBlockMemFirstVSLater);
disp(pBlockMemFirstVSLater);

% session groups
groups  = {'allSessions', 'All sessions'; 'firstSessions', 'First sessions'; 'laterSessions', 'Later sessions'};

% loop through the session groups
for iGroup = 1:size(groups, 1)
    
    % select sessions based on group
    if strcmp(groups{iGroup, 1}, 'allSessions')
        bMaskSession = true(size(bFirstSession));
    elseif strcmp(groups{iGroup, 1}, 'firstSessions')
        bMaskSession = bFirstSession;
    elseif strcmp(groups{iGroup, 1}, 'laterSessions')
        bMaskSession = ~bFirstSession;
    end
    
    % select data
    thisData = blockMemPerf(bMaskSession, :);
    
    % report
    fprintf('\nAnalysis of memory performance in: %s (n = %d).\n', groups{iGroup, 2}, sum(bMaskSession));
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
    ax = axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 4]);
    if strcmp(groups{iGroup, 1}, 'laterSessions')
        set(ax, 'position', [0.3, 1.4, 4, 4], 'YAxisLocation', 'right');
    end
    x = 1:size(thisData, 2);
    m = mean(thisData, 1, 'omitnan');
    ste = LK_ste(thisData);
    n = sum(~isnan(thisData));
    hold on;
    patch([x, fliplr(x)], [m - ste, fliplr(m + ste)], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, m, '-', 'Color', [0, 0, 0], 'LineWidth', 1);
    if strcmp(groups{iGroup, 1}, 'firstSessions') || strcmp(groups{iGroup, 1}, 'laterSessions')
        text(x(pBlockMemFirstVSLater < 0.05), repmat(0.925, 1, sum(pBlockMemFirstVSLater < 0.05)), '\_', 'Color', [0, 0, 0], 'FontSize', 12, 'horizontalalignment', 'center');
    end
    t1 = text(0.95, 0, sprintf('\\itn\\rm = %d', size(thisData, 1)), 'units', 'normalized', 'horizontalalignment', 'right', 'verticalalignment', 'bottom'); % sample size
    xl = xlabel('Block');
    yl = ylabel('Memory performance');
    tl = title(groups{iGroup, 2});
    set([ax, t1, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    set(ax, 'xlim', [min(x), max(x)] + [-1, 1], 'ylim', [0.6, 1], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    % save figure
    LK_print(f, strcat(paths.save, 'allMemPerfInBlocks_', groups{iGroup, 1}, '_20230227'), '-dtiff', '-r300');
end

%% memory performance during first vs. second task sessions (paired analysis)

% unique subjects
bMask                   = cellfun(@(x) contains(x(end), {'a', 'b', 'c'}), subjects(:, 1));
strippedSubjects        = subjects(:, 1);
strippedSubjects(bMask) = cellfun(@(x) x(1:end - 1), strippedSubjects(bMask), 'UniformOutput', false);
uniqueSubjects          = unique(strippedSubjects);

% memory performance per task session, organized by subject
memPerfPerSess  = nan(size(uniqueSubjects, 1), param.maxNumSess);
for iSub = 1:size(allMemPerf, 1)
    
    % subject and session index
    subIdx = contains(uniqueSubjects, strippedSubjects{iSub, 1});
    if subjects{iSub}(end) == 'b'
        sessIdx = 2;
    elseif subjects{iSub}(end) == 'c'
        sessIdx = 3;
    else
        sessIdx = 1;
    end
    
    % average memory performance for this subject and session
    memPerfPerSess(subIdx, sessIdx) = mean(allMemPerf{iSub}, 'omitnan');
end

% stats: memory performance between first and second sessions in those
% subjects that performed both first and second sessions
b1And2 = ~isnan(memPerfPerSess(:, 1)) & ~isnan(memPerfPerSess(:, 2));
[~, p1And2, ~, stats1And2] = ttest(memPerfPerSess(b1And2, 1), memPerfPerSess(b1And2, 2));
fprintf('\nAre the memory-performance values different between 1st and 2nd sessions? t(%d) = %.3f, p = %.3f.\n', stats1And2.df, stats1And2.tstat, p1And2);

% figure: memory performance between first and second session
f = figure('units', 'centimeters', 'position', [5, 5, 5, 6]);
ax = axes('units', 'centimeters', 'position', [1.5, 1.4, 3, 4]);
x = [1, 2];
hold on;
for iUSub = 1:size(memPerfPerSess, 1)
    if b1And2(iUSub, 1) == 1 % if the subject performed a first and second session
        plot(x, [memPerfPerSess(iUSub, 1), memPerfPerSess(iUSub, 2)], '-', 'Color', [0, 0, 0]);
    end
end
t1 = text(1, 0, sprintf('\\itn\\rm = %d', sum(b1And2)), 'units', 'normalized', 'horizontalalignment', 'right', 'verticalalignment', 'bottom'); % sample size
xl = xlabel('Session');
yl = ylabel('Memory performance');
if p1And2 < 0.001
    tl = title('\itP\rm < 0.001');
else
    tl = title('');
end
set([gca, t1, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
set(gca, 'xlim', x + [-0.1, 0.1], 'xtick', x, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
% save figure
LK_print(f, strcat(paths.save, 'memPerfPerSess_1And2_20230227'), '-dtiff', '-r300');
