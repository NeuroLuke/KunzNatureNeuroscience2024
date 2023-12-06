%==========================================================================
% This script identifies the trial halves that lead to the largest
% difference in memory performance between the early and the late data
% half.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.save          = 'E:\OpenField\BehAnalysis_20211022\MemoryImprovement\';
mkdir(paths.save);

% settings
param               = [];
param.maxNumTrials  = 160; % maximum number of trials
param.smoothType    = 'movmean';
param.smoothFacAll  = 11; % smoothing factor for the object-unspecific analysis
param.smoothFacObj  = 3; % smoothing factor for the object-specific analysis
param.myRNG         = 111;
rng(param.myRNG);

% subjects
subjects            = load(strcat(paths.subjects, 'subjects.mat'));
subjects            = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% random locations within arena for normalization
dt                  = [];
dt.maxR             = 5000; % maximum radius
dt.minR             = 0; % minimum radius
dt.N                = 10000001; % number of locations to create
dt.centerX          = 0; % arena center x
dt.centerY          = 0; % arena center y
locsInArena         = LK_RandomPointsInCircle(dt);

% preallocate
idxTrial4BestSep        = nan(size(subjects, 1), 1);
idxTrial4BestSepObjWise = nan(size(subjects, 1), 1);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nSubject: %s.\n', subjects{iSub});
    
    %% load data
    
    % load trial information
    behFile = dir(strcat(paths.beh, subjects{iSub, 1}, '\trialInfo.mat'));
    trials  = load(fullfile(behFile.folder, behFile.name));
    trials  = trials.trialInfoMacrotime;
    fprintf('Original number of trials: %d.\n', size(trials, 1));
    
    % restrict to maximum number of trials
    if size(trials, 1) > param.maxNumTrials
        trials  = trials(1:param.maxNumTrials, :);
    end
    fprintf('Adjusted number of trials: %d.\n', size(trials, 1));
    
    %% drop error and memory performance
    
    % drop error and memory performance per trial
    dropError       = sqrt((trials.xCorrect - trials.xResponse) .^ 2 + (trials.yCorrect - trials.yResponse) .^ 2);
    dropErrorSurro  = pdist2([trials.xCorrect, trials.yCorrect], locsInArena);
    trials.memPerf  = sum(dropError < dropErrorSurro, 2) ./ sum(~isnan(dropErrorSurro), 2);
    
    %% timepoint of largest memory improvement (overall)
    
    % smooth the memory performances to get rid of outliers (irrespective
    % of object identity)
    thisMemPerf = smoothdata(trials.memPerf, 1, param.smoothType, param.smoothFacAll);
    
    % loop through possible separations of trials
    statsPerSep = nan(size(trials, 1), 1);
    for iTrial = 1:size(trials, 1) - 1
    
        % t-test for this separation of trials (high t-values represent a
        % large improvement)
        [~, ~, ~, stats]        = ttest2(thisMemPerf(iTrial + 1:end), thisMemPerf(1:iTrial)); % 2nd half > 1st half
        statsPerSep(iTrial, 1)  = stats.tstat;
    end
        
    % do not consider the first or last trial as possible best trial
    statsPerSep([1; numel(statsPerSep)])    = nan;
    
    % last trial of half 1 leading to the strongest separation of early-bad
    % versus late-good memory
    [maxVal, maxIdx]            = max(statsPerSep);
    idxTrial4BestSep(iSub, 1)   = trials.TrialIdx(maxIdx);
    
    % if there is actually no learning, assign all trials to the early half
    if maxVal < 0
        idxTrial4BestSep(iSub, 1)   = Inf;
    end
    
    %% timepoint of largest memory improvement (object-wise)
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 30, 12], 'visible', 'off');
    tiledlayout(2, 4, 'Padding', 'compact');
    
    % loop through objects
    for iObj = min(trials.Object):max(trials.Object)
        
        % trials with this object
        thisTrials  = trials(trials.Object == iObj, :);
        thisMemPerf = trials.memPerf(trials.Object == iObj, :);
        
        % smooth the memory performance to get rid of outliers
        % (object-wise)
        thisMemPerf = smoothdata(thisMemPerf, 1, param.smoothType, param.smoothFacObj);
        
        % loop through possible separations of trials
        statsPerSepObjWise  = nan(size(thisTrials, 1), 1);
        for iTrial = 1:size(thisTrials, 1) - 1
        
            % t-test for this separation of trials
            [~, ~, ~, stats]                = ttest2(thisMemPerf(iTrial + 1:end), thisMemPerf(1:iTrial)); % 2nd half > 1st half
                        
            % store the stats for each separation
            statsPerSepObjWise(iTrial, 1)   = stats.tstat;
        end
                
        % do not consider the first or last trial as possible best trial
        statsPerSepObjWise([1; numel(statsPerSepObjWise)])  = nan;
        
        % last trial of half 1 leading to the strongest memory improvement
        [maxVal, maxIdx]                        = max(statsPerSepObjWise);
        idxTrial4BestSepObjWise(iSub, iObj + 1) = thisTrials.TrialIdx(maxIdx);
        
        % if there is actually no learning, assign all trials to the early
        % half
        if maxVal < 0
            idxTrial4BestSepObjWise(iSub, iObj + 1) = Inf;
        end
        
        % add to figure
        nexttile;
        hold on;
        % left y-axis
        yyaxis left;
        plot(thisTrials.TrialIdx, thisMemPerf, 'Color', [0, 0, 0]);
        if idxTrial4BestSepObjWise(iSub, iObj + 1) ~= Inf
            xline(idxTrial4BestSepObjWise(iSub, iObj + 1), ':', 'LineWidth', 0.5, 'Color', [1, 0, 0, 0.1]);
        end
        xl = xlabel('Trial index', 'units', 'normalized', 'position', [0.5, -0.06]);
        yl = ylabel('Memory performance');
        set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
        % right y-axis
        yyaxis right;
        plot(thisTrials.TrialIdx, statsPerSepObjWise, 'Color', [0.5, 0.5, 0.5]);
        yl = ylabel('Statistic (\itt\rm)');
        tl = title(sprintf('Object #%d: trial #%d', iObj + 1, idxTrial4BestSepObjWise(iSub, iObj + 1)));
        set([gca, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
        ax = gca;
        ax.YAxis(1).Color = [0, 0, 0];
        ax.YAxis(2).Color = [0.5, 0.5, 0.5];
        set(gca, 'xlim', [min(ax.XLim), max(ax.XLim)], 'xtick', [min(ax.XLim), max(ax.XLim)], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    end
    
    % save figure
    LK_print(f, strcat(paths.save, subjects{iSub}, '_MemoryImprovementPerSeparationObjWise'), '-dpng', '-r450');
        
    % close open figures
    close all;
end

%% save output

% convert data into table
idxTrial4BestSep        = table(subjects(:, 1), idxTrial4BestSep, 'VariableNames', {'Subject', 'ObjAll'});
idxTrial4BestSepObjWise = table(subjects(:, 1), idxTrial4BestSepObjWise(:, 1), idxTrial4BestSepObjWise(:, 2), idxTrial4BestSepObjWise(:, 3), ...
    idxTrial4BestSepObjWise(:, 4), idxTrial4BestSepObjWise(:, 5), idxTrial4BestSepObjWise(:, 6), idxTrial4BestSepObjWise(:, 7), idxTrial4BestSepObjWise(:, 8), ...
    'VariableNames', {'Subject', 'Obj0', 'Obj1', 'Obj2', 'Obj3', 'Obj4', 'Obj5', 'Obj6', 'Obj7'});

% save tables
save(strcat(paths.save, 'idxTrial4BestSep'), 'idxTrial4BestSep');
save(strcat(paths.save, 'idxTrial4BestSepObjWise'), 'idxTrial4BestSepObjWise');
