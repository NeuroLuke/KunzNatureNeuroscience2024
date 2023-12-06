%==========================================================================
% This script plots navigation paths examines the speed profile during the 
% trials.
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
paths.save          = 'E:\OpenField\BehAnalysis_20211022\NavigationPath\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% settings
param               = [];
param.trialPhases   = {...
    'ITI', 1, rgb('purple'); ... % name, index, color
    'Cue', 2, rgb('red'); ...
    'Retrieval', 3, rgb('orange'); ...
    'Feedback', 4, rgb('limegreen'); ...
    'Reencoding', 5, rgb('darkgreen')};
param.speedCutoff   = 0.001; % speed cutoff
param.speedYLim     = [-100, 1000]; % for plotting
param.maxNumTrials  = 160;
param.bPlot         = false;

% subjects
subjects            = load(strcat(paths.subjects, 'subjects.mat'));
subjects            = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% preallocate results
allRes              = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nSubject: %s.\n', subjects{iSub});
    
    %% load data
    
    % load behavioral information
    behFile  	= dir(strcat(paths.beh, subjects{iSub, 1}, '\behInfoRes10Hz.mat')); % resampled behavioral information
    beh         = load(fullfile(behFile.folder, behFile.name));
    beh         = beh.behInfoRes;
    
    % speed profile
    speed                       = smoothdata(beh.speed, 1, 'movmean', 1);
    speed(beh.trialPhase <= 2)  = 0; % no movement during ITI and cue
    bMove                       = speed > param.speedCutoff; % movement, if above a particular speed cutoff
    
    % load trial information
    trialFile   = dir(strcat(paths.beh, subjects{iSub, 1}, '\trialInfo.mat')); % trial information in behavioral time
    trials      = load(fullfile(trialFile.folder, trialFile.name));
    trials      = trials.trialInfo;
    
    % restrict trial number
    if size(trials, 1) > param.maxNumTrials
        trials  = trials(1:param.maxNumTrials, :);
    end
    
    %% navigation path and speed profile
    
    % loop through trials
    for iTrial = 1:size(trials, 1)
        
        % timepoints from this trial
        bThisTrial          = beh.trialIdx == iTrial;
        
        % speed profile of this trial
        thisRes             = [];
        thisRes.idx         = [iSub, iTrial];
        thisRes.time        = beh.time(bThisTrial); % time during this trial
        thisRes.speed       = speed(bThisTrial); % speed during this trial
        thisRes.bMove       = bMove(bThisTrial); % movement during this trial
        thisRes.ITI         = trials.ITI(iTrial); % time of ITI start
        thisRes.Cue         = trials.Cue(iTrial);
        thisRes.Retrieval   = trials.Retrieval(iTrial);
        thisRes.Feedback    = trials.Feedback(iTrial);
        thisRes.Reencoding  = trials.Reencoding(iTrial);
        thisRes.Grab        = trials.Grab(iTrial);
        
        % collect speed profile across subjects and trials
        allRes              = cat(1, allRes, thisRes);
        
        %% figure: navigation path
        
        % do the plotting for only a fraction of the trials
        if (param.bPlot == false) || (mod(iTrial, 10) ~= 1)
            continue;
        end
        fprintf('\tTrial: %d.\n', iTrial);
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'visible', 'off');
        axes('units', 'centimeters', 'position', [1.25, 1.25, 4, 4]);
        hold on;
        % boundary
        alpha = 0:0.001:(2*pi);
        plot(cos(alpha) * 5000, sin(alpha) * 5000, 'k-', 'LineWidth', 2);
        % loop through trial phases
        for iTP = 3:size(param.trialPhases, 1)
            bThisPhase  = beh.trialPhase == param.trialPhases{iTP, 2};
            bMask       = bThisTrial & bThisPhase;
            plot(beh.x(bMask), beh.y(bMask), '.', 'Color', param.trialPhases{iTP, 3});
        end
        % enhance axes
        xl = text(0, 5750, 'x', 'HorizontalAlignment', 'center');
        yl = text(-5750, 0, 'y', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        % enhance axes
        set(gca, 'ydir', 'reverse');
        set([xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
        axis equal off;
        % save figure
        LK_print(f, strcat(paths.save, subjects{iSub}, '_Trial', num2str(iTrial), '_NavigationPath_20230407'), '-dpng', '-r300');
        
        %% figure: speed profile
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 16, 4], 'visible', 'off');
        axes('units', 'centimeters', 'position', [1.5, 1, 14, 2.5]);
        hold on;
        % indicate movement periods
        patch([beh.time(bThisTrial); flipud(beh.time(bThisTrial))], ...
            [range(param.speedYLim) .* bMove(bThisTrial) + param.speedYLim(1); param.speedYLim(1) .* ones(size(bMove(bThisTrial)))], [0.9, 0.9, 0.9], 'EdgeColor', 'none');
        % loop through trial phases
        for iTP = 1:size(param.trialPhases, 1)
            % speed during this trial phase
            bThisPhase  = beh.trialPhase == param.trialPhases{iTP, 2};
            bMask       = bThisTrial & bThisPhase;
            plot(beh.time(bMask), speed(bMask), '-', 'Color', param.trialPhases{iTP, 3});
            % indicate phase onsets
            plot(trials{iTrial, param.trialPhases(iTP)}, param.speedYLim(2), 'p', 'Color', param.trialPhases{iTP, 3}, 'MarkerFaceColor', param.trialPhases{iTP, 3});
        end
        % enhance axes
        xl = xlabel('Time (s)', 'units', 'normalized', 'position', [0.5, -0.13]);
        yl = ylabel({'Speed', '(vu/s)'}, 'units', 'normalized', 'position', [-0.025, 0.5]);
        set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
        set(gca, 'xlim', [min(beh.time(bThisTrial)), max(beh.time(bThisTrial))], 'xtick', [min(beh.time(bThisTrial)), max(beh.time(bThisTrial))], 'xticklabel', [round(min(beh.time(bThisTrial))), round(max(beh.time(bThisTrial)))], ...
            'ylim', param.speedYLim, 'ytick', [0, param.speedYLim(2)], 'tickdir', 'out', 'ticklength', [0.01, 0.01], 'layer', 'top');
        % save figure
        LK_print(f, strcat(paths.save, subjects{iSub}, '_Trial', num2str(iTrial), '_SpeedProfile_20230407'), '-dpng', '-r300');
        
        %% close figures
        close all;
    end
end

%% speed profile across all trials

% normalized time
normT       = 0:0.02:1;

% preallocate resampled speed for each trial and phase
resSpeed    = cell(size(allRes, 1), size(param.trialPhases, 1));
fracMove    = nan(size(allRes, 1), size(param.trialPhases, 1));

% normalize the time line within each trial phase
for iTrial = 1:size(allRes, 1)
    
    % loop through trial phases
    for iPhase = 1:size(param.trialPhases, 1)
        
        % start and end time point of this trial phase
        if strcmp(param.trialPhases{iPhase}, 'ITI')
            startEnd    = [allRes(iTrial).ITI, allRes(iTrial).Cue];
        elseif strcmp(param.trialPhases{iPhase}, 'Cue')
            startEnd    = [allRes(iTrial).Cue, allRes(iTrial).Retrieval];
        elseif strcmp(param.trialPhases{iPhase}, 'Retrieval')
            startEnd    = [allRes(iTrial).Retrieval, allRes(iTrial).Feedback];
        elseif strcmp(param.trialPhases{iPhase}, 'Feedback')
            startEnd    = [allRes(iTrial).Feedback, allRes(iTrial).Reencoding];
        elseif strcmp(param.trialPhases{iPhase}, 'Reencoding')
            startEnd    = [allRes(iTrial).Reencoding, allRes(iTrial).Grab];
        end
        
        % data from this phase
        bPhase  = allRes(iTrial).time >= min(startEnd) & allRes(iTrial).time <= max(startEnd);
        
        % interpolate onto new time axis
        if sum(bPhase) > 1
            resTime                     = linspace(min(startEnd), max(startEnd), numel(normT));
            resSpeed{iTrial, iPhase}    = interp1(allRes(iTrial).time(bPhase), allRes(iTrial).speed(bPhase), resTime);
        else
            resSpeed{iTrial, iPhase}    = nan(size(normT));
        end
        
        % fraction of this trial phase during which the subject moves
        fracMove(iTrial, iPhase) = sum(allRes(iTrial).bMove(bPhase)) / sum(bPhase);
    end
end

%% figure: speed profiles per trial phase and movement probability per trial phase

% loop through trial phases
for iPhase = 1:size(param.trialPhases, 1)
    
    % resampled speed profiles from all trials for this trial phase
    allResSpeed = cell2mat(resSpeed(:, iPhase));
    
    % metrics
    med     = median(allResSpeed, 1, 'omitnan');
    p25     = prctile(allResSpeed, 25, 1);
    p75     = prctile(allResSpeed, 75, 1);
    minima  = min(allResSpeed, [], 1);
    maxima  = max(allResSpeed, [], 1);
    mMove   = mean(allResSpeed > param.speedCutoff, 1, 'omitnan'); % movement probability
    steMove = LK_ste(allResSpeed > param.speedCutoff);
    
    % create figure: speed profile and movement probability
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 8]);
    
    % speed profile
    ax1 = axes('units', 'centimeters', 'position', [1.7, 3.4, 4, 4]);
    hold on;
    plot(normT, minima, '.', 'Color', [0.8, 0.8, 0.8], 'MarkerSize', 2);
    plot(normT, maxima', '.', 'Color', [0.8, 0.8, 0.8], 'MarkerSize', 2);
    plot([normT; normT], [p25; p75], '-', 'Color', [0.6, 0.6, 0.6]);
    plot(normT, med, '.', 'Color', param.trialPhases{iPhase, 3}, 'MarkerSize', 10);
    yl = ylabel('Speed (vu/s)', 'units', 'normalized', 'position', [-0.3, 0.5]);
    tl = title(param.trialPhases{iPhase, 1}, 'Color', param.trialPhases{iPhase, 3});
    set([gca, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
    set(gca, 'xtick', [], 'ylim', [0, 1200], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    
    % movement probability
    ax2 = axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 1.25]);
    hold on;
    patch([normT, fliplr(normT)], [mMove - steMove, fliplr(mMove + steMove)], param.trialPhases{iPhase, 3}, 'EdgeColor', 'none');
    plot(normT, mMove, '-', 'Color', param.trialPhases{iPhase, 3}, 'LineWidth', 1);
    xl = xlabel('Normalized time');
    yl = ylabel({'Movement', 'probability'}, 'units', 'normalized', 'position', [-0.163, 0.5]);
    set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
    set(gca, 'ylim', [0, 1], 'ytick', [0, 1], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    linkaxes([ax1, ax2], 'x');
    
    % save figure
    LK_print(f, strcat(paths.save, param.trialPhases{iPhase, 1}, '_AllTrials_SpeedProfile_20230407'), '-dpng', '-r300');
end

%% figure: percent of trials during which the subject moves

% settings for this figure
xEdges  = 0:5:100;

% loop through trial phases
for iPhase = 1:size(param.trialPhases, 1)
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 5, 4]);
    axes('units', 'centimeters', 'position', [1.4, 1.6, 3, 2]);
    hold on;
    histogram(100 * fracMove(:, iPhase), xEdges, 'normalization', 'probability', ...
        'edgecolor', param.trialPhases{iPhase, 3}, 'facecolor', param.trialPhases{iPhase, 3}, 'facealpha', 0.2, 'linewidth', 1);
    xline(mean(100 * fracMove(:, iPhase), 'omitnan'), '-', 'Color', [0, 0, 0], 'LineWidth', 2);
    xl = xlabel('Movement time (%)');
    yl = ylabel('% trials', 'units', 'normalized', 'position', [-0.25, 0.5]);
    ax = get(gca);
    set(gca, 'xlim', [min(xEdges), max(xEdges)], 'ylim', round(ax.YLim, 2), 'ytick', round(ax.YLim, 2), 'yticklabel', {round(ax.YLim, 2) * 100}, ...
        'tickdir', 'out', 'ticklength', [0.04, 0.04], 'Layer', 'Top');
    set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
    
    % save figure
    LK_print(f, strcat(paths.save, param.trialPhases{iPhase, 1}, '_AllTrials_FracMove_20230407'), '-dpng', '-r300');
end
