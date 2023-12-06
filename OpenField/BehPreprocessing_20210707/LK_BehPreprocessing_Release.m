%==========================================================================
% This script preprocesses the behavioral data and aligns them with (i) the
% macrodata timeline and (ii) with the microdata timeline.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));

% paths
paths       = [];
paths.beh   = 'E:\OpenField\Beh_20201215\';
paths.micro = 'E:\OpenField\Micro_20201215\';
paths.macro = 'E:\OpenField\Macro_20201215\';
paths.save  = 'E:\OpenField\BehPreprocessing_20210707\';

% parameters
param                   = [];
param.newTimeRes      	= 0.1; % new temporal resolution (Hz) for resampled behavioral data
param.direc.maxAbsYaw   = 32768; % maximum absolute Unreal yaw value
param.macro.triggerChan = 'Trigger'; % macro trigger channel
param.micro.triggerChan = 'ainp2'; % microwire trigger channel

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
    
    % report
    fprintf('\n\n==== SUBJECT: %s.\n', subjects{iSub});
    
    % save path for this subject
    thisSavePath    = strcat(paths.save, subjects{iSub});
    mkdir(thisSavePath);
    
    % skip this subject if the final output exists already
    if exist(strcat(thisSavePath, '\behInfoRes', num2str(1 / param.newTimeRes), 'Hz.mat'), 'file') > 0
        continue;
    end
    
    %% behavioral triggers
    
    % load behavioral triggers and estimate inter-trigger intervals
    t                   = load(strcat(paths.beh, subjects{iSub}, '\triggers.mat'));
    behTrigTimepoints   = t.triggers;
    behTrigITI          = diff(behTrigTimepoints);
    
    %% macro triggers and their congruence with behavioral triggers
    
    % load macro trigger data
    e           = load(strcat(paths.macro, subjects{iSub}, '\eeg_original.mat'));
    trigIdx     = strcmp(e.eeg_original.label, param.macro.triggerChan);
    macroData   = e.eeg_original.trial{1}(trigIdx, :);
    macroData   = smoothdata(macroData, 2, 'movmean', 10); % remove tiny bumps
    macroTime   = e.eeg_original.time{1}; % first sample has the timepoint 0
    
    % determine macrodata trigger timepoints
    thisThresh          = max(macroData) / 2; % threshold for detecting triggers
    logIdx              = diff(macroData > thisThresh) == 1;
    macroTrigTimepoints = transpose(macroTime(logIdx));
    macroTrigITI        = diff(macroTrigTimepoints);
    
    % correlation between inter-trigger-intervals
    [rho, pval]         = corr(behTrigITI, macroTrigITI);
    
    % report
    fprintf('Correlation between "behtrigITI" and "macrotrigITI": rho = %.6f, p = %.6f.\n', rho, pval);
    if rho < 0.999
        error('The correlation is lower than rho = 0.999');
    end
        
    %% micro triggers and their congruence with behavioral triggers
    
    % load microwire trigger data
    microFile   = fullfile(paths.micro, subjects{iSub}, param.micro.triggerChan, 'datacut.mat');
    
    % load microwire trigger data if they exist
    if exist(microFile, 'file') > 0
        
        % microwire trigger data and time
        m           = load(microFile);
        microData   = m.data;
        microData 	= smoothdata(microData, 2, 'movmean', 11); % remove tiny bumps
        microTime   = ((1:size(microData, 2)) - 1) ./ m.sr; % first sample has the timepoint 0
        
        % determine microwire trigger timepoints
        thisThresh          = max(microData) / 2; % threshold for detecting triggers
        logIdx              = diff(microData > thisThresh) == 1;
        microTrigTimepoints = transpose(microTime(logIdx));
        microTrigITI        = diff(microTrigTimepoints);
        
        % correlation between inter-trigger intervals
        [rho, pval]         = corr(behTrigITI, microTrigITI);
        
        % report
        fprintf('Correlation between "behTrigITI" and "microTrigITI": rho = %.6f, p = %.6f.\n', rho, pval);
        if rho < 0.999
            error('The correlation is lower than rho = 0.999');
        end
    else
        
        % use vector of nans for the trigger timepoints
        microTrigTimepoints = nan(size(behTrigTimepoints));
        microTrigITI        = nan(size(behTrigITI));
        
        % report
        warning('The subject does not have microwire data. Using nans for "microTrigTimepoints" and "microTrigITI".');
    end
        
    %% behavioral data
    
    % load behavioral data
    b           = load(strcat(paths.beh, subjects{iSub}, '\data.mat'));
    origData    = b.data;
    
    %% ==================================================================== process trial information ==> trialInfo
    
    % total number of trials
    numTrials   = sum(contains(origData, 'TrialPhase2Counter ='));
    
    % report
    fprintf('\nExtracting relevant trial information ==> "trialInfo.mat".\n');
    fprintf('Number of trials the subject completed: %d.\n', numTrials);
    
    % create "trials" variable with relevant trial information: trialidx,
    % object, ITI onset, cue onset, retrieval onset, feedback onset,
    % re-encoding onset, grab onset, x-correct, y-correct, x-drop, y-drop,
    % drop error
    trials      = nan(numTrials, 13);
    bTestPhase  = false; % whether the test phase has started
    for iLine = 1:size(origData, 1)
        
        % information from this line
        thisLineText        = origData{iLine, :};
        thisLineSplits      = strsplit(thisLineText, ' ');
        if iLine < size(origData, 1)
            nextLineText    = origData{iLine + 1, :};
            nextLineSplits  = strsplit(nextLineText, ' ');
        end
        
        % look for information
        if regexp(thisLineText, 'TrialPhase2Counter =')
            
            % trial index
            trialIdx                = str2double(thisLineSplits{4});
            trials(trialIdx, 1)     = trialIdx;
            
            % note that the test phase has begun
            bTestPhase              = true;
            
        elseif regexp(thisLineText, ' Cue ')
            
            % object identity
            trials(trialIdx, 2)     = str2double(thisLineSplits{5});
        
        elseif regexp(thisLineText, 'ITI_start')
            
            % timepoint of ITI start
            trials(trialIdx, 3)     = str2double(thisLineSplits{3});

        elseif regexp(thisLineText, 'CUE_start')
            
            % timepoint of cue start
            trials(trialIdx, 4)     = str2double(thisLineSplits{3});
            
        elseif regexp(thisLineText, 'MyExp2Spawner.bReadyDrop True')
            
            % timepoint of retrieval start
            trials(trialIdx, 5)     = str2double(nextLineSplits{3}); % note: timepoint is in the next line
            
        elseif regexp(thisLineText, ' Drop ')
            
            % response location (x/y)
            trials(trialIdx, [11, 12])  = [str2double(thisLineSplits{7}), str2double(thisLineSplits{9})];
            
        elseif regexp(thisLineText, 'DropBeep')
            
            % timepoint of feedback
            trials(trialIdx, 6)     = str2double(thisLineSplits{2});
            
        elseif regexp(thisLineText, ' Show ')
            
            % skip if test phase has not begun yet
            if bTestPhase == false
                continue;
            end
            
            % correct location (x/y)
            trials(trialIdx, [9, 10])   = [str2double(thisLineSplits{7}), str2double(thisLineSplits{9})];
            
        elseif regexp(thisLineText, 'how accurately placed')
            
            % drop error, computed by the paradigm
            trials(trialIdx, 13)    = str2double(thisLineSplits{6});
            
        elseif regexp(thisLineText, 'STOP_SMILEYFB')
            
            % timepoint of re-encoding start
            trials(trialIdx, 7)     = str2double(thisLineSplits{4});
            
        elseif regexp(thisLineText, 'GrabBeep')
            
            % timepoint of grab
            trials(trialIdx, 8)     = str2double(thisLineSplits{2});
        end
    end
    
    % remove trials that do not have a cue onset
    if any(isnan(trials(:, 4)))
        trials  = trials(~isnan(trials(:, 4)), :);
        warning('Removing trials from the "trials" that do not have a cue onset.');
    end
    
    % sanity check that time points are strictly increasing
    tmpTimes    = transpose(trials(:, 3:8));
    if any(diff(tmpTimes(:)) < 0)
        error('Timepoints in "trials" not strictly increasing.');
    end
    
    %% align behavioral time with macro/microwire time
    
    % report
    fprintf('\nConverting the behavioral time points in "trials" to macro/microwire time using a trial-specific conversion procedure.\n');
    
    % trial information with macro/microwire time
    trialsMacrotime = trials;
    trialsMicrotime = trials;
    
    % trial-specific alignment
    for iTrial = 1:size(trials, 1) - 1
        
        % linear regression: behavioral to MACRO time
        P                               = polyfit(behTrigTimepoints(iTrial:iTrial + 1), macroTrigTimepoints(iTrial:iTrial + 1), 1);
        trialsMacrotime(iTrial, 3:8)    = trials(iTrial, 3:8) .* P(1) + P(2);
        if iTrial == (size(trials, 1) - 1)
            trialsMacrotime(iTrial + 1, 3:8)    = trials(iTrial + 1, 3:8) .* P(1) + P(2); % for the last trial, use the conversion of the second to last trial
        end
        
        % linear regression: behavioral to MICRO time
        P                               = polyfit(behTrigTimepoints(iTrial:iTrial + 1), microTrigTimepoints(iTrial:iTrial + 1), 1);
        trialsMicrotime(iTrial, 3:8)    = trials(iTrial, 3:8) .* P(1) + P(2);
        if iTrial == (size(trials, 1) - 1)
            trialsMicrotime(iTrial + 1, 3:8)    = trials(iTrial + 1, 3:8) .* P(1) + P(2);
        end
    end
    
    %% convert into "trialInfo" table and save
    
    % variable names
    varNames    = {'TrialIdx', 'Object', 'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab', 'xCorrect', 'yCorrect', 'xResponse', 'yResponse', 'DropError'};
    
    % trial information with behavioral, macro, and micro time
    trialInfo           = array2table(trials, 'VariableNames', varNames);
    trialInfoMacrotime  = array2table(trialsMacrotime, 'VariableNames', varNames);
    trialInfoMicrotime  = array2table(trialsMicrotime, 'VariableNames', varNames);
        
    % save trial information
    save(strcat(thisSavePath, '\trialInfo'), 'trialInfo', 'trialInfoMacrotime', 'trialInfoMicrotime');
    
    %% ==================================================================== process behavioral information ==> behInfo
    
    % report
    fprintf('\nExtracting relevant behavioral information ==> "behInfo.mat".\n');
    
    % preallocate output
    behData  	= nan(size(origData, 1), 7); % time, x, y, z, yawUnreal, trial idx, trial phase
    
    % initialize trial index and trial phase as nans
    trialIdx    = nan;
    trialPhase  = nan;
    
    % loop through text lines
    for iLine = 1:size(origData, 1) - 1
        
        % information from this and the next line
        thisLineText = origData{iLine, :};
        nextLineText = origData{iLine + 1, :}; % needed for yawUnreal
        
        % update current trial index
        if regexp(thisLineText, ' TrialPhase2Counter ')
            splits    	= strsplit(thisLineText, ' ');
            trialIdx   	= str2double(splits{4});
        end
        
        % update current trial phase
        if regexp(thisLineText, 'ITI_start')
            trialPhase 	= 1; % ITI
        elseif regexp(thisLineText, 'CUE_start')
            trialPhase 	= 2; % cue
        elseif regexp(thisLineText, 'bReadyDrop True')
            trialPhase 	= 3; % retrieval
        elseif regexp(thisLineText, 'DropBeep')
            trialPhase 	= 4; % feedback
        elseif regexp(thisLineText, 'STOP_SMILEYFB')
            trialPhase 	= 5; % re-encoding
        elseif regexp(thisLineText, 'GrabBeep')
            trialPhase 	= 6; % grab
        end
        
        % collect all information whenever subject is at a specific
        % location
        if regexp(thisLineText, 'Location X')
            
            % time and xyz-location information
            splits              = strsplit(thisLineText, ' ');
            behData(iLine, 1:4) = [str2double(splits{3}), str2double(splits{6}), str2double(splits{8}), str2double(splits{10})];
            
            % heading information (from the next line)
            splits              = strsplit(nextLineText, ' ');
            behData(iLine, 5)   = str2double(splits{10});
            if ~contains(nextLineText, 'Yaw')
                error('Next line does not contain "Yaw" information.'); % sanity-check
            end
            
            % trial-index information
            behData(iLine, 6)   = trialIdx;
            
            % trial-phase information
            behData(iLine, 7)   = trialPhase;
        end
    end
    
    % remove lines that do not have a time stamp
    fprintf('Size of "behData" before cutting lines without a time stamp: \t%d x %d.\n', size(behData));
    behData = behData(~isnan(behData(:, 1)), :);
    fprintf('Size of "behData" after cutting lines without a time stamp: \t%d x %d.\n', size(behData));
    
    %% align behavioral time with macro/microwire time
    
    % report
    fprintf('\nConverting the behavioral time points in "behData" to macro/microwire time using a trial-specific conversion procedure.\n');
    
    % behavioral information with macro/microwire time
    behDataMacrotime    = behData;
    behDataMicrotime    = behData;
    
    % trial-specific alignment
    for iTrial = 1:size(trials, 1) - 1
        
        % behavioral timepoints from this inter-trigger interval
        logIdx     	= behData(:, 1) >= behTrigTimepoints(iTrial) & behData(:, 1) <= behTrigTimepoints(iTrial + 1);
        if iTrial == 1
            logIdx 	= behData(:, 1) <= behTrigTimepoints(iTrial + 1);
        elseif iTrial == (size(trials, 1) - 1)
            logIdx  = behData(:, 1) >= behTrigTimepoints(iTrial);
        end
        
        % conversion of behavioral to MACRO time
        P                           = polyfit(behTrigTimepoints(iTrial:iTrial + 1), macroTrigTimepoints(iTrial:iTrial + 1), 1);
        behDataMacrotime(logIdx, 1) = behData(logIdx, 1) .* P(1) + P(2);
        
        % conversion of behavioral to MICRO time
        P                           = polyfit(behTrigTimepoints(iTrial:iTrial + 1), microTrigTimepoints(iTrial:iTrial + 1), 1);
        behDataMicrotime(logIdx, 1) = behData(logIdx, 1) .* P(1) + P(2);
    end
    
    % check validity of behavioral data in macrotime
    if sum(diff(behDataMacrotime(:, 1)) < 0) > 0 || any(isnan(behDataMacrotime(:, 1)))
        error('There is a problem with "behDataMacrotime".');
    end
    
    % check validity of behavioral data in microtime
    if exist(microFile, 'file') > 0
        if sum(diff(behDataMicrotime(:, 1)) < 0) > 0 || any(isnan(behDataMicrotime(:, 1)))
            error('There is a problem with "behDataMicrotime".');
        end
    end
    
    %% convert into "behInfo" table and add additional variables
    
    % variable names
    varNames    = {'time', 'x', 'y', 'z', 'yawUnreal', 'trialIdx', 'trialPhase'};
    
    % loop through different data groups
    groups      = {'behInfo'; 'behInfoMacrotime'; 'behInfoMicrotime'};
    for iGroup = 1:size(groups, 1)
        
        % select data and convert into table
        if strcmp(groups{iGroup}, 'behInfo')
            thisInfo    = array2table(behData, 'VariableNames', varNames);
        elseif strcmp(groups{iGroup}, 'behInfoMacrotime')
            thisInfo    = array2table(behDataMacrotime, 'VariableNames', varNames);
        elseif strcmp(groups{iGroup}, 'behInfoMicrotime')
            thisInfo    = array2table(behDataMicrotime, 'VariableNames', varNames);
        end
        
        % durations
        thisInfo.durations  = [diff(thisInfo.time); median(diff(thisInfo.time), 'omitnan')]; % add a median duration at the end
        
        % translational distances and speed
        thisInfo.distances  = [sqrt(diff(thisInfo.x) .^ 2 + diff(thisInfo.y) .^ 2); 0];
        thisInfo.speed  	= thisInfo.distances ./ thisInfo.durations;
        
        % angular yaws, distances, and speed
        thisInfo.yaw            = (thisInfo.yawUnreal ./ param.direc.maxAbsYaw) .* pi; % in radians
        thisInfo.yawDistances   = abs([angdiff(thisInfo.yaw); 0]);
        thisInfo.yawSpeed       = thisInfo.yawDistances ./ thisInfo.durations;
        
        % distribute data
        if strcmp(groups{iGroup}, 'behInfo')
            behInfo             = thisInfo;
        elseif strcmp(groups{iGroup}, 'behInfoMacrotime')
            behInfoMacrotime    = thisInfo;
        elseif strcmp(groups{iGroup}, 'behInfoMicrotime')
            behInfoMicrotime    = thisInfo;
        end
    end
    
    %% convert onto new timeline (e.g., 10 Hz)
    
    % loop through different data groups
    groups  = {'behInfoRes'; 'behInfoMacrotimeRes'; 'behInfoMicrotimeRes'};
    for iGroup = 1:size(groups, 1)
        
        % select data
        if strcmp(groups{iGroup}, 'behInfoRes')
            thisInfo    = behInfo;
        elseif strcmp(groups{iGroup}, 'behInfoMacrotimeRes')
            thisInfo    = behInfoMacrotime;
        elseif strcmp(groups{iGroup}, 'behInfoMicrotimeRes')
            thisInfo    = behInfoMicrotime;
        end
        
        % new timeline
        newTimeMin  = LK_floor(min(thisInfo.time), 1);
        newTimeMax  = LK_ceil(max(thisInfo.time), 1);
        newTime  	= transpose(newTimeMin:param.newTimeRes:newTimeMax);
        fprintf('New timeline: minimum = %.3f (from %.3f), maximum = %.3f (from %.3f).\n', min(newTime), min(thisInfo.time), max(newTime), max(thisInfo.time));
        
        % initialize the resampled data
        thisInfoRes = array2table(nan(size(newTime, 1), size(thisInfo, 2)), 'VariableNames', thisInfo.Properties.VariableNames);
        
        % time
        thisInfoRes.time    = newTime;
        
        % loop through timepoints and estimate resampled data
        for iTime = 1:size(thisInfoRes, 1) - 1
            
            % original timepoints that fall into this new time window
            bTWOI   = thisInfo.time >= thisInfoRes.time(iTime) & thisInfo.time < thisInfoRes.time(iTime + 1);
            
            % continue if there is no data for this time window
            if sum(bTWOI) == 0
                continue;
            end
            
            % x, y, and z
            thisInfoRes.x(iTime)    = mean(thisInfo.x(bTWOI), 'omitnan');
            thisInfoRes.y(iTime)    = mean(thisInfo.y(bTWOI), 'omitnan');
            thisInfoRes.z(iTime)    = mean(thisInfo.z(bTWOI), 'omitnan');
            
            % yaw
            tmpYaw                          = thisInfo.yawUnreal(bTWOI) ./ param.direc.maxAbsYaw .* pi; % convert Unreal units to [-pi, pi]
            tmpYaw                          = circ_mean(tmpYaw(~isnan(tmpYaw))); % remove nans before calculating the circular mean
            thisInfoRes.yawUnreal(iTime)    = tmpYaw ./ pi .* param.direc.maxAbsYaw; % convert from [-pi, pi] to Unreal units
                        
            % trial index and trial phase
            thisInfoRes.trialIdx(iTime)     = mode(thisInfo.trialIdx(bTWOI));
            thisInfoRes.trialPhase(iTime)   = mode(thisInfo.trialPhase(bTWOI));
        end
        
        % replace possible nans with preceding entry
        for iTime = 2:size(thisInfoRes, 1) % skip first line
            for iVar = 1:size(thisInfoRes, 2)
                if isnan(thisInfoRes{iTime, iVar})
                    thisInfoRes{iTime, iVar}    = thisInfoRes{iTime - 1, iVar};
                end
            end
        end
        
        % durations
        thisInfoRes.durations   = [diff(thisInfoRes.time); median(diff(thisInfoRes.time), 'omitnan')]; % add a median duration at the end
        
        % translational distances and speed
        thisInfoRes.distances   = [sqrt(diff(thisInfoRes.x) .^ 2 + diff(thisInfoRes.y) .^ 2); 0]; % add a zero distance at the end
        thisInfoRes.speed       = thisInfoRes.distances ./ thisInfoRes.durations;
        
        % angular distances and speed
        thisInfoRes.yaw             = thisInfoRes.yawUnreal ./ param.direc.maxAbsYaw .* pi; % convert Unreal units to [-pi, pi]
        thisInfoRes.yawDistances    = abs([angdiff(thisInfoRes.yaw); 0]); % angular distances
        thisInfoRes.yawSpeed        = thisInfoRes.yawDistances ./ thisInfoRes.durations; % angular speed      
                
        % distribute data
        if strcmp(groups{iGroup}, 'behInfoRes')
            behInfoRes         	= thisInfoRes;
        elseif strcmp(groups{iGroup}, 'behInfoMacrotimeRes')
            behInfoMacrotimeRes = thisInfoRes;
        elseif strcmp(groups{iGroup}, 'behInfoMicrotimeRes')
            behInfoMicrotimeRes = thisInfoRes;
        end
    end
    
    %% save output
    
    % save behavioral information, original timeline
    save(strcat(thisSavePath, '\behInfo'), 'behInfo', 'behInfoMacrotime', 'behInfoMicrotime');
    
    % save behavioral information, new timeline
    save(strcat(thisSavePath, '\behInfoRes', num2str(1 / param.newTimeRes), 'Hz'), 'behInfoRes', 'behInfoMacrotimeRes', 'behInfoMicrotimeRes');
        
    %% close figures
    close all;
end
