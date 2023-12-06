%==========================================================================
% This script aligns macrodata triggers and microwire triggers to get the
% spike times in macrotime.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% paths
paths           = [];
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.spikes    = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.macro     = 'E:\OpenField\MacroPreprocessing_20210627\';

% settings
param                   = [];
param.macro.trigFile    = 'Trigger.mat';
param.micro.trigChan    = 'ainp2';
param.micro.trigFile    = 'datacut.mat';
param.spikes.file       = 'times_datacut.mat';
param.output.name     	= 'cluster_class_macrotime'; % name of the output file

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects 	= subjects.subjects(strcmp(subjects.subjects(:, 2), 'Microwire'));
fprintf('Converting microtime into macrotime. Number of microwire sessions: %d.\n', size(subjects, 1));

% preallocate control output
allRes      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\nSUBJECT: %s.\n\n', subjects{iSub});
    
    %% load macro and microwire trigger data
    
    % data for macro triggers
    macroTrig       = load(strcat(paths.macro, subjects{iSub}, filesep, param.macro.trigFile));
    macroTrigTime   = macroTrig.Trigger.time{1}; % first sample has the timepoint 0
    macroTrigData   = macroTrig.Trigger.trial{1};
    macroTrigData   = smoothdata(macroTrigData, 2, 'movmean', 10); % smooth to get rid of tiny bumps
    
    % macro triggers and inter-trigger-intervals
    thisThresh          = max(macroTrigData) / 2; % threshold for detecting triggers
    logIdx              = diff(macroTrigData > thisThresh) == 1;
    macroTrigTimepoints = transpose(macroTrigTime(logIdx));
    macroTrigITI        = diff(macroTrigTimepoints);
    
    % data for micro triggers
    microTrig       = load(strcat(paths.spikes, subjects{iSub}, filesep, param.micro.trigChan, filesep, param.micro.trigFile));
    microTrigData   = microTrig.data;
    microTrigData   = smoothdata(microTrigData, 2, 'movmean', 10); % smooth to get rid of tiny bumps
    microTrigTime   = ((1:size(microTrigData, 2)) - 1) ./ microTrig.sr; % first sample has the timepoint 0
    
    % micro triggers and inter-trigger-intervals
    thisThresh          = max(microTrigData) / 2; % threshold for detecting triggers
    logIdx              = diff(microTrigData > thisThresh) == 1;
    microTrigTimepoints = transpose(microTrigTime(logIdx));
    microTrigITI        = diff(microTrigTimepoints);
    
    % correlation between inter-trigger-intervals
    [rho, pval]         = corr(macroTrigITI, microTrigITI);
    if rho < 0.999
        error('The correlation is lower than rho = 0.999');
    end
    
    %% spike analysis
    
    % channels
    chans   = dir(strcat(paths.spikes, subjects{iSub}, '\chan*'));
    [~, I]  = sort(cellfun(@str2double, replace({chans.name}', lettersPattern, '')));
    chans   = chans(I);
    
    %% loop through channels
    for iChan = 1:size(chans, 1)
        
        % load wave-clus output for this channel
        wcFile  = strcat(chans(iChan).folder, filesep, chans(iChan).name, filesep, param.spikes.file);
        if exist(wcFile, 'file') > 0
            
            % load information about cluster class and spike time
            t   = load(wcFile);
        else
            % report
            warning('Wave-clus did not extract spikes from this channel.\n');
            
            % save empty matrix as output if no data available
            cluster_class_macrotime     = [];
            save(fullfile(chans(iChan).folder, chans(iChan).name, param.output.name), param.output.name);
            
            continue;
        end
        
        %% loop through spikes and assign macrotime
        
        % cluster, microtime (msec), macrotime (sec)
        cluster_class_macrotime     = nan(size(t.cluster_class, 1), 3);
        for iSpike = 1:size(t.cluster_class, 1)
            
            % microwire time of this spike in seconds(!)
            thisSpikeTime   = t.cluster_class(iSpike, 2) / 1000; % convert (msec) into (sec)
            
            % find micro trigger that is directly before the spike
            timeDiff        = microTrigTimepoints - thisSpikeTime;
            [~, trigIdx]    = max(timeDiff(timeDiff < 0));
            
            % continue, if spike time not between first cue and last cue
            if isempty(trigIdx) || trigIdx == length(microTrigTimepoints)
                continue;
            end
            
            % relative timing of spike within the two neighboring micro
            % triggers
            relTiming           = (thisSpikeTime - microTrigTimepoints(trigIdx)) / (microTrigTimepoints(trigIdx + 1) - microTrigTimepoints(trigIdx));
            
            % corresponding macro time
            thisSpikeMacroTime  = macroTrigTimepoints(trigIdx) + relTiming * (macroTrigTimepoints(trigIdx + 1) - macroTrigTimepoints(trigIdx));
            
            % collect output
            cluster_class_macrotime(iSpike, 1:2)    = t.cluster_class(iSpike, :); % cluster class and microwire timing
            cluster_class_macrotime(iSpike, 3)      = thisSpikeMacroTime; % macrodata timing
        end
        
        %% save output
        
        % save output file
        outputFile  = strcat(chans(iChan).folder, filesep, chans(iChan).name, filesep, param.output.name, '.mat');
        save(outputFile, param.output.name);
        
        %% bookkeeping
        
        % number of units available for analyses
        uniqueNames = unique(cluster_class_macrotime(~isnan(cluster_class_macrotime(:, 1)), 1));
        numUnits    = sum(uniqueNames ~= 0); % ignore "rest" clusters
        
        % original number of units
        uniqueNamesOrig = unique(t.cluster_class(~isnan(t.cluster_class(:, 1)), 1));
        numUnitsOrig    = sum(uniqueNamesOrig ~= 0);
        
        % percentage of spikes between the first and last cue
        fprintf('\t\tBetween first and last cue: %d spikes (%.2f%% of all spikes).\n', sum(~isnan(cluster_class_macrotime(:, 1))), ...
            100 * sum(~isnan(cluster_class_macrotime(:, 1))) / size(cluster_class_macrotime, 1));        
        percSpikesWithinTOI = 100 * sum(~isnan(cluster_class_macrotime(:, 1))) / size(cluster_class_macrotime, 1);
        
        %% collect information across units
        
        % results from this unit
        thisRes                         = [];
        thisRes.idx                     = [iSub, iChan];
        thisRes.numUnits                = numUnits;
        thisRes.numUnitsOrig            = numUnitsOrig;        
        thisRes.percSpikesWithinTOI     = percSpikesWithinTOI;
        
        % all results
        allRes                          = cat(1, allRes, thisRes);
    end
end

%% save results
save(strcat(paths.spikes, 'allRes'), 'allRes');

%% load results
r = load(strcat(paths.spikes, 'allRes.mat'));

%% report

% all indices
allIdx  = cell2mat({r.allRes.idx}');

% total number of units for analyses
allNumUnits = cell2mat({r.allRes.numUnits}');
totNumUnits = sum(allNumUnits);
fprintf('Total number of units: %d.\n', totNumUnits);

% total number of units originally
allNumUnitsOrig     = cell2mat({r.allRes.numUnitsOrig}');
totNumUnitsOrig     = sum(allNumUnitsOrig);
fprintf('Total number of units originally: %d.\n', totNumUnitsOrig);

% percentage of spikes between first and last cue
allPerc     = cell2mat({r.allRes.percSpikesWithinTOI}');
fprintf('Minimum percentage of spikes kept after transformation into macrotime: %.2f.\n', min(allPerc));
fprintf('Maximum percentage of spikes kept after transformation into macrotime: %.2f.\n', max(allPerc));
