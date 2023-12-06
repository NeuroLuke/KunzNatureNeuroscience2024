%==========================================================================
% This script performs behavioral and ripple analyses regarding the initial
% encoding period.
%
% It estimates the duration of the initial encoding for each subject and
% performs analyses with regard to ripple rates during the initial encoding
% events.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions'));

% paths
paths         	= [];
paths.subjects 	= 'E:\OpenField\SubjectData_20210616\';
paths.beh     	= 'E:\OpenField\Beh_20201215\';
paths.trial    	= 'E:\OpenField\BehPreprocessing_20210707\';
paths.ripples 	= 'E:\OpenField\RipplesMacro_20210614\20220320_FBS_80to140Hz_tMF2\';
paths.save     	= 'E:\OpenField\InitialEncoding_20230226\20230226\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% settings
param                       = [];
param.arenaR                = 5000; % radius of the arena
param.ripple.timeRes        = 0.001; % time resolution of the ripple-PSTH
param.ripple.timeEdges      = -5:param.ripple.timeRes:5; % time edges for the ripple-PSTH
param.ripple.timeCenters    = movmean(param.ripple.timeEdges, 2, 'endpoints', 'discard'); % time centers for the ripple-PSTH
param.ripple.smoothType     = 'Gaussian'; % type of smoothing
param.ripple.smoothTime     = 0.5; % time period of smoothing kernel

% subjects
subjects     	= load(strcat(paths.subjects, 'subjects.mat'));
subjects        = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% preallocate
allSessRes      = [];
allRippleRes    = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\n==== SUBJECT: %s.\n', subjects{iSub});
    
    %% data from initial encoding
    
    % file that contains the initial-encoding information
    fileName    = strcat(paths.save, subjects{iSub}, '_InitialEncoding.mat');
    
    % get data from initial encoding
    if exist(fileName, 'file')
        
        % load the previously extracted information
        enc     = load(fileName);
        enc     = enc.enc;
    else
        
        % load original behavioral data
        b           = load(strcat(paths.beh, subjects{iSub}, '\data.mat'));
        origData    = b.data;
        
        % initialize encoding data
        enc = table('Size', [0, 7], 'VariableNames', {'Time', 'PlayerX', 'PlayerY', 'PlayerYaw', 'Object', 'ObjX', 'ObjY'}, ...
            'VariableTypes', {'double', 'double', 'double', 'double', 'char', 'double', 'double'});
        
        % initialize original values
        thisTime        = nan;
        thisPlayerX     = nan;
        thisPlayerY     = nan;
        thisPlayerYaw   = nan;
        thisObject      = 'none';
        thisObjX        = nan;
        thisObjY        = nan;
        
        % reset status of beginning data reading
        bStarted        = false;
        
        % extract information from the logfile
        for iLine = 1:size(origData, 1)
            
            % information from this line
            thisLine    = origData{iLine};
            thisInfo    = split(thisLine, ' ');
            
            % check whether initial encoding has started
            if bStarted == false && contains(thisLine, 'ScriptLog: 1')
                bStarted = true;
                fprintf('\tInitial encoding started: %s.\n', thisLine);
            elseif bStarted == true
            else
                continue;
            end
            
            % information about time
            thisTime    = str2double(thisInfo{3});
            
            % information about player location
            if contains(thisLine, 'Location')
                thisPlayerX = str2double(thisInfo{6});
                thisPlayerY = str2double(thisInfo{8});
            end
            
            % information about player orientation
            if contains(thisLine, 'YawUnit')
                thisPlayerYaw   = str2double(thisInfo{10});
            end
            
            % information about object (start)
            if contains(thisLine, 'Show')
                thisObject  = replace(replace(thisInfo{5}, 't_', ''), 'Neuro2.Obj.', '');
                thisObjX    = str2double(thisInfo{7});
                thisObjY    = str2double(thisInfo{9});
            end
            
            % collect information
            thisEnc = cell2table({thisTime, thisPlayerX, thisPlayerY, thisPlayerYaw, {thisObject}, thisObjX, thisObjY}, ...
                'VariableNames', {'Time', 'PlayerX', 'PlayerY', 'PlayerYaw', 'Object', 'ObjX', 'ObjY'});
            enc = cat(1, enc, thisEnc);
            
            % finish if testing period has started
            if contains(thisLine, 'Begin Phase2')
                fprintf('\tPhase 2 reached: %s.\n', thisLine);
                break;
            end
        end
        
        % save information about initial encoding
        save(fileName, 'enc');
    end
    
    %% output: duration of the initial encoding period
    
    % use the first player movement as start
    distances   = sqrt(diff([0; enc.PlayerX]) .^ 2 + diff([0; enc.PlayerY]) .^ 2);
    bFirstMove  = find(distances > 0, 1);
    
    % results from this session
    thisSessRes                     = [];
    thisSessRes.subject             = subjects{iSub};
    thisSessRes.duration            = max(enc.Time) - min(enc.Time); % in seconds
    thisSessRes.durationFirstMove   = max(enc.Time) - enc.Time(bFirstMove); % in seconds
    
    % concatenate the results from different session
    allSessRes                      = cat(1, allSessRes, thisSessRes);
    
    %% figure: navigation path and object locations during initial encoding
    
    % colors
    uniqueObjects   = unique(enc.Object(~strcmp(enc.Object, 'none')), 'stable');
    colors          = flipud(distinguishable_colors(numel(uniqueObjects)));
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8.5, 6], 'visible', 'off');
    axes('units', 'centimeters', 'position', [0.5, 0.75, 5, 5]);
    set(gca, 'xlim', [-param.arenaR, param.arenaR], 'ylim', [-param.arenaR, param.arenaR]);
    axis off;
    hold on;
    plot(param.arenaR .* cos(0:0.001:(2*pi)), param.arenaR .* sin(0:0.001:(2*pi)), '-', 'Color', [0, 0, 0], 'LineWidth', 4);
    p = nan(size(uniqueObjects, 1), 1);
    for iO = 1:size(uniqueObjects, 1)
        logIdx = strcmp(enc.Object, uniqueObjects{iO});
        plot(enc.PlayerX(logIdx), enc.PlayerY(logIdx), '-', 'Color', colors(iO, :), 'LineWidth', 1); % player path
        p(iO) = plot(enc.ObjX(logIdx), enc.ObjY(logIdx), 'x', 'Color', colors(iO, :), 'MarkerSize', 10, 'LineWidth', 3); % object location
        hLine = findobj(p(iO), 'type', 'line');
        set(hLine, 'MarkerSize', 10);
    end
    lg = legend(p, {'Object 1', 'Object 2', 'Object 3', 'Object 4', 'Object 5', 'Object 6', 'Object 7', 'Object 8'}, ...
        'location', [0.74, 0.275, 0.2, 0.5], 'units', 'normalized', 'box', 'off', 'FontSize', 10); % legend
    text(-0.05, -0.05, ['Subject: ', num2str(iSub)], 'fontunits', 'centimeters', 'fontsize', 0.2, 'units', 'normalized'); % indicate subject number
    xl = text(0, -5500, 'x', 'horizontalalignment', 'center');
    yl = text(-5750, 0, 'y', 'verticalalignment', 'middle');
    set([xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
    LK_print(f, strcat(paths.save, subjects{iSub}, '_InitialEncodingPath_20230226'), '-dtiff', '-r450');
    
    %% ripple effects during the initial encoding period
    
    % estimate the iEEG timepoints of the behavioral events during the
    % initial encoding period
    trialFile       = dir(strcat(paths.trial, subjects{iSub, 1}, '\trialInfo.mat'));
    trials          = load(fullfile(trialFile.folder, trialFile.name));
    coef            = polyfit(trials.trialInfo.Cue, trials.trialInfoMacrotime.Cue, 1); % regression between behavioral times and behavioral times in iEEG time
    enc.Macrotime   = enc.Time .* coef(1) + coef(2);
    
    % ripple channels
    rippleChans     = LK_dir(strcat(paths.ripples, subjects{iSub}));
    
    % loop through ripple channels
    for iRippleChan = 1:size(rippleChans, 1)
        
        %% ripple rates during initial encoding and the main task
        
        % load ripples
        r       = load(fullfile(rippleChans(iRippleChan).folder, rippleChans(iRippleChan).name, 'ripples.mat'));
        ripples = r.ripples.ripples;
        
        % ripple timepoints
        rippleTimepoints    = [ripples.peakTime]';
        
        % mean ripple rate during initial encoding (IE)
        numRipplesIE    = sum(rippleTimepoints >= min(enc.Macrotime) & rippleTimepoints <= max(enc.Macrotime));
        durationIE      = max(enc.Macrotime) - min(enc.Macrotime);
        rippleRateIE    = numRipplesIE / durationIE;
        
        % mean ripple rate the the main task (MT)
        numRipplesMT    = sum(rippleTimepoints >= min(trials.trialInfoMacrotime.Cue) & rippleTimepoints <= max(trials.trialInfoMacrotime.Cue));
        durationMT      = max(trials.trialInfoMacrotime.Cue) - min(trials.trialInfoMacrotime.Cue);
        rippleRateMT    = numRipplesMT / durationMT;
        
        %% time-resolved ripple rates around the initial-encoding events
        
        % ripple rates time-locked to the "grab" events during initial
        % encoding
        uniqueObjects   = unique(enc.Object(~strcmp(enc.Object, 'none')), 'stable');
        grabTime        = nan(size(uniqueObjects, 1), 1); % "grab" timepoints in original behavioral time
        grabMacrotime   = nan(size(uniqueObjects, 1), 1); % "grab" timepoints in iEEG time
        for iO = 1:size(uniqueObjects, 1)
            bThisObject             = strcmp(enc.Object, uniqueObjects{iO});
            grabTime(iO, 1)         = max(enc.Time(bThisObject));
            grabMacrotime(iO, 1)    = max(enc.Macrotime(bThisObject));
        end
        
        % compute grab-locked ripple PSTHs
        grabTimeEdges   = grabMacrotime + param.ripple.timeEdges;
        grabTimeCenters = grabMacrotime + param.ripple.timeCenters;
        numRipplesPSTH  = nan(size(grabTimeCenters));
        durationsPSTH   = nan(size(grabTimeCenters));
        for iO = 1:size(grabTimeEdges, 1)
            numRipplesPSTH(iO, :)   = histcounts(rippleTimepoints, grabTimeEdges(iO, :));
            durationsPSTH(iO, :)    = round(diff(grabTimeEdges(iO, :), [], 2), 4);
        end
        rippleRatePSTH      = numRipplesPSTH ./ durationsPSTH;
        meanRippleRatePSTH  = mean(rippleRatePSTH, 1); % average across initial-encoding events
        
        %% output: ripple rates during initial encoding
        
        % results for this channel
        thisRippleRes                       = [];
        thisRippleRes.subject            	= subjects{iSub};
        thisRippleRes.rippleChan          	= rippleChans(iRippleChan).name;
        thisRippleRes.idx                  	= [iSub, iRippleChan];
        
        % mean ripple rates
        thisRippleRes.rippleRateIE       	= rippleRateIE;
        thisRippleRes.rippleRateMT      	= rippleRateMT;
        
        % time-resolved ripple rates locked to grab
        thisRippleRes.rippleRatePSTH    	= rippleRatePSTH;
        thisRippleRes.meanRippleRatePSTH  	= meanRippleRatePSTH;
        
        % collect all results
        allRippleRes  = cat(1, allRippleRes, thisRippleRes);
    end
end

%% duration of all initial encoding periods

% report
fprintf('\nResults:\n');

% report durations of the initial encoding periods
allDur          = [allSessRes.duration]' ./ 60; % entire duration (minutes)
allDurFirstMove = [allSessRes.durationFirstMove]' ./ 60; % duration since first movement (minutes)
LK_ReportMeanAndSEM_20220322('Duration of initial encoding periods (min)', allDur);
LK_ReportMeanAndSEM_20220322('Duration of initial encoding periods since first movement (min)', allDurFirstMove);

% create
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
ax = axes('units', 'centimeters', 'position', [1.5, 1.5, 4, 4]);
hold on;
plot(sort(allDur), '.', 'Color', [0, 0, 0]);
LK_yline(mean(allDur), '-', [0.75, 0.75, 0.75], 'bottom');
xl = xlabel('Session');
yl = ylabel('Duration (min)');
set([ax, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(ax, 'ytick', 0:2:10, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
LK_print(f, strcat(paths.save, 'InitialEncodingDurations_20230315'), '-dtiff', '-r300');

%% ripple rates during initial encoding vs. main task

% mean ripple rate during initial encoding and the main task
allRippleRateIE = [allRippleRes.rippleRateIE]';
allRippleRateMT = [allRippleRes.rippleRateMT]';

% compare mean ripple rates
[~, p, ~, stats] = ttest(allRippleRateIE, allRippleRateMT);
fprintf('Are ripple rates during initial encoding higher than during the actual task? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

% correlation between ripple rates from initial encoding vs. main task
[rho, pval] = corr(allRippleRateIE, allRippleRateMT, 'Type', 'Pearson');
fprintf('Are the ripple rates during initial encoding correlated with the ripple rates during the main task? rho = %.3f, p = %.3f.\n', rho, pval);

% figure: correlation between ripple rates during initial encoding vs.
% during the main task
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [2, 2, 3.5, 3.5]);
hold on;
plot(allRippleRateIE, allRippleRateMT, '.', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10);
P = polyfit(allRippleRateIE, allRippleRateMT, 1);
plot(allRippleRateIE, P(1) .* allRippleRateIE + P(2), 'r-');
t1 = text(0.975, 0.025, sprintf('\\itr\\rm = %.3f', rho), 'units', 'normalized', 'horizontalalignment', 'right', 'verticalalignment', 'bottom');
xl = xlabel({'Ripple rate (Hz)', 'Initial encoding'});
yl = ylabel({'Main task', 'Ripple rate (Hz)'});
set([gca, t1, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
LK_print(f, strcat(paths.save, 'CorrRippleRateIEvsMT_20230315'), '-dpng', '-r300');

%% time-resolved ripple rates around the "grab" events

% average grab-related ripple PSTH per channel
allMeanRRPSTH       = cell2mat({allRippleRes.meanRippleRatePSTH}');

% smooth the PSTHs
allMeanRRPSTH       = smoothdata(allMeanRRPSTH, 2, param.ripple.smoothType, param.ripple.smoothTime / param.ripple.timeRes + 1);

% test against average ripple rate during initial encoding
[~, p, ~, stats]    = ttest(allMeanRRPSTH - repmat(allRippleRateIE, 1, size(allMeanRRPSTH, 2)));

% mean and standard error
m   = mean(allMeanRRPSTH, 1);
ste = LK_ste(allMeanRRPSTH);

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.7, 1.5, 4, 4]);
hold on;
patch([param.ripple.timeCenters, fliplr(param.ripple.timeCenters)], [m - ste, fliplr(m + ste)], [0.7, 0.7, 0.7], 'EdgeColor', 'none');
plot(param.ripple.timeCenters, m, 'Color', [0, 0, 0], 'LineWidth', 1);
LK_yline(mean(allRippleRateIE), '-', [1, 0, 0], 'top');
grid minor;
xl = xlabel('Time (s)');
yl = ylabel('Ripple rate (Hz)');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'xlim', [-2, 2], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
LK_print(f, strcat(paths.save, 'ChanRippleRate_AbsTimeFromGrab_InitialEncoding_20230315'), '-dpng', '-r300');
