%==========================================================================
% This script examines the coactivity of object cells and place cells
% during hippocampal ripples.
%
% Coactivity is quantified via a z-score (Sosa et al., Neuron, 2020) or a
% Pearson correlation.
%
% The coactivity maps of associative cell pairs are tested against
% (a) chance (zero or surrogates);
% (b) baseline coactivity maps (before/after the actual coactivity map);
% (c) coactivity maps of non-associative cell pairs.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all; clc;

% settings
%- random number generator
param                       = [];
param.myRNG                 = 999;
rng(param.myRNG);
%- subjects
param.subject.fileName      = 'subjectdata_20210616.mat';
%- behavior
param.beh.maxNumTrials      = 160; % maximum number of trials
param.beh.speedCutoff       = 0.001; % speed cutoff
param.beh.memFile           = 'idxTrial4BestSepObjWise.mat'; % object-wise separation of trials regarding before vs. after learning
%- microwires
param.micro.c4aFileName     = 'LK_Cluster4Analysis_20210521.mat';
param.micro.ccmFileName     = 'cluster_class_macrotime.mat';
%- ripples
param.ripple.folder         = 'RipplesMacro_20210614';
param.ripple.specs          = '20220320_FBS_80to140Hz_tMF2';
param.ripple.type           = 'ripples';
param.ripple.fileName       = 'ripples.mat';
param.ripple.timeRes        = 0.001; % temporal resolution for estimating firing rates during ripples
param.ripple.timeEdges   	= -2:param.ripple.timeRes:2; % time-bin edges around ripples for estimating firing rates
param.ripple.timeCenters    = movmean(param.ripple.timeEdges, 2, 'endpoints', 'discard');
%- coactivity
param.coac.type             = 'z'; % zscore ('z') or Pearson correlation ('r')
param.coac.twoiRes          = 0.1; % temporal resolution of the sliding time-window-of-interest
param.coac.twoiShift        = 0.005; % temporal shift of the sliding time-window-of-interest
param.coac.basePre          = [-0.25, 0.25] - 0.5; % time window for the pre-baseline coactivity map
param.coac.basePost         = [-0.25, 0.25] + 0.5; % time window for the post-baseline coactivity map
param.coac.target           = [-0.25, 0.25]; % time window for the target coactivity map
param.coac.twoiOns          = (min(param.coac.basePre) - param.coac.twoiRes / 2):param.coac.twoiShift:(max(param.coac.basePost) - param.coac.twoiRes / 2);
param.coac.twoi             = round(transpose([param.coac.twoiOns; param.coac.twoiOns + param.coac.twoiRes]), 4);
param.coac.time             = transpose(mean(param.coac.twoi, 2));
if strcmp(param.coac.type, 'z') % zscore
    param.coac.clim.ret     = [-0.4, 0.4]; % c-limits for 2D-coactivity plots during retrieval
    param.coac.clim.enc     = [-0.6, 0.6]; % c-limits for 2D-coactivity plots during encoding
elseif strcmp(param.coac.type, 'r') % correlation
    param.coac.clim.ret     = [-0.2, 0.2];
    param.coac.clim.enc     = [-0.3, 0.3];
end
param.coac.bPlotIndividuals = true; % whether to plot individual coactivity maps
param.coac.numSurrogates    = 2001; % number of surrogates

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.beh           = 'E:\OpenField\BehPreprocessing_20210707\';
paths.memory        = 'E:\OpenField\BehAnalysis_20211022\MemoryImprovement\';
paths.spikes        = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.objCells      = 'E:\OpenField\CueObjAnalysis_20210812\20211116\';
paths.placeCells    = 'E:\OpenField\PlaceCellAnalysis_20210511\20211116_Loc25x25_facThreshPF0.75\';
paths.ripple        = strcat('E:\OpenField\', param.ripple.folder, '\', param.ripple.specs, '\');
paths.functions     = 'E:\OpenField\Functions\';
paths.fieldtrip     = 'E:\fieldtrip\fieldtrip-20210614\';
paths.save          = strcat('E:\OpenField\NeuronalCoactivityDuringRipples_20211025\20230424_', ...
    param.ripple.folder, '_', param.ripple.specs, '\', ...
    param.ripple.type, '_', param.coac.type, '_tR', num2str(param.coac.twoiRes * 1000), '_tS', num2str(param.coac.twoiShift * 1000), '\');
paths.pics          = strcat(paths.save, 'Pics\');

% add functions
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
ft_warning off;

% subjects
subjects  	= load(strcat(paths.subjects, 'subjects.mat'));
subjects  	= subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings

% create save path
if ~exist(paths.pics, 'dir')
    mkdir(paths.pics);
end

% save settings
save(strcat(paths.save, 'settings'));

%% cell types

% report
fprintf('\nLoading cell classifications.\n');

% object cells
objRes      = load(strcat(paths.objCells, 'results.mat'), 'allRes');
objCells    = load(strcat(paths.objCells, 'cellClassification.mat'));
bObjCell    = strcmp(objCells.cellClassification.cellType, 'ObjectCell');
fprintf('Number of object cells: %d (out of %d).\n', sum(bObjCell), size(bObjCell, 1));

% place cells
placeRes    = load(strcat(paths.placeCells, 'results.mat'), 'allRes', 'place');
placeCells  = load(strcat(paths.placeCells, 'cellClassification.mat'));
bPlaceCell  = strcmp(placeCells.cellClassification.cellType, 'PlaceCell');
fprintf('Number of place cells: %d (out of %d).\n', sum(bPlaceCell), size(bPlaceCell, 1));

%% random locations to convert drop errors into memory-performance values

% random locations
cfg         = [];
cfg.maxR    = 5000;
cfg.minR    = 0;
cfg.N       = 10000001;
cfg.centerX = 0;
cfg.centerY = 0;
randLocs    = LK_RandomPointsInCircle(cfg);

%% memory improvement

% load information about time of memory improvement
memImprovement  = load(strcat(paths.memory, param.beh.memFile));

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByInspection  = [];
exByNoRipples   = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\n\n==========================================================\n');
    fprintf('SUBJECT: %s.\n', subjects{iSub});
    
    % skip this subject, if it is not a microwire subject
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        fprintf(2, 'Skipping this subject because it does not have microwires.\n');
        continue;
    end
    
    %% behavioral data: trial information
    
    % report
    fprintf('\nLoading and processing the behavioral data.\n');
    
    % load trial information
    trialInfo           = load(strcat(paths.beh, subjects{iSub}, '\trialInfo.mat'));
    trialInfo           = trialInfo.trialInfoMacrotime; % in macrotime
    fprintf('\tOriginal number of cue periods: %d.\n', sum(~isnan(trialInfo.Cue)));
    
    % restrict to maximum number of trials
    if size(trialInfo, 1) > param.beh.maxNumTrials
        trialInfo       = trialInfo(1:param.beh.maxNumTrials, :);
    end
    fprintf('\tCorrected number of cue periods: %d.\n', sum(~isnan(trialInfo.Cue)));
    
    % memory performance per trial
    dropError           = sqrt((trialInfo.xResponse - trialInfo.xCorrect) .^ 2 + (trialInfo.yResponse - trialInfo.yCorrect) .^ 2);
    dropErrorSurro      = pdist2([trialInfo.xCorrect, trialInfo.yCorrect], randLocs);
    trialInfo.memPerf   = sum(dropError < dropErrorSurro, 2) ./ sum(~isnan(dropErrorSurro), 2);
    
    %% behavioral data: navigation information
    
    % load behavioral information to identify movement periods
    behInfo             = load(strcat(paths.beh, subjects{iSub}, '\behInfoRes10Hz.mat'));
    behInfo             = behInfo.behInfoMacrotimeRes;
    
    % identify movement periods
    bMovement           = behInfo.speed > param.beh.speedCutoff;
    bStart              = diff([bMovement(1); bMovement]) == 1;
    bEnd                = diff([bMovement; bMovement(end)]) == -1;
    movementStartEnd    = [behInfo.time(bStart), behInfo.time(bEnd)];
    fprintf('\tNumber of movement periods: %d.\n', size(movementStartEnd, 1));
    
    % trials before/after memory improvement
    if contains(param.beh.memFile, 'idxTrial4BestSepObjWise')
        
        % trials before and after memory improvement, object-wise
        bB4Learn = nan(size(trialInfo, 1), 1); % boolean indicating whether a trial is before or after learning
        for iObj = min(trialInfo.Object):max(trialInfo.Object)
            
            % separating trial for this object
            bC  = strcmp(memImprovement.idxTrial4BestSepObjWise.Properties.VariableNames, ['Obj', num2str(iObj)]); % relevant column
            bR  = strcmp(memImprovement.idxTrial4BestSepObjWise.Subject, subjects{iSub}); % relevant row
            thisIdxTrial4BestSep    = memImprovement.idxTrial4BestSepObjWise{bR, bC};
            
            % separation for this object
            bThisObj            = trialInfo.Object == iObj;
            bThisB4Learn        = trialInfo.TrialIdx(bThisObj) <= thisIdxTrial4BestSep; % trials before learning (includes the separating trial)
            bB4Learn(bThisObj)  = bThisB4Learn;
        end
        trialInfo.bB4Learn      = bB4Learn;
    else
        error("Name of memory-improvement file not valid.");
    end
    
    %% ripple data
    
    % ripple channels
    fprintf('\nLoading and processing the ripple data.\n');
    rippleChans = LK_dir(strcat(paths.ripple, subjects{iSub}, '\*'));
    
    % double-check whether the subject has hippocampal channels
    if isempty(rippleChans)
        error('Subject does not have any ripple channels.');
    end
    
    % assemble the ripple data from all ripple channels
    rippleData      = repmat(struct, size(rippleChans, 1), 1);
    for iRippleChan = 1:size(rippleChans, 1)
                
        % ripples
        r               = load(fullfile(rippleChans(iRippleChan).folder, rippleChans(iRippleChan).name, param.ripple.fileName));
        ripples         = r.ripples.ripples; % includes information about the ripple peak timepoint
        
        % remove ripples before the first or after the last cue
        bRipples2Keep   = [ripples.startTime]' >= min(trialInfo.Cue) & [ripples.endTime]' <= max(trialInfo.Cue);
        ripples     	= ripples(bRipples2Keep);
        
        % report
        fprintf('\tRipple channel: %s. Keeping === %d === of %d ripples.\n', rippleChans(iRippleChan).name, sum(bRipples2Keep), size(bRipples2Keep, 1));
        
        % ripple timepoints
        rippleTimepoints   	= [ripples.peakTime]';
                
        % assign information to each ripple
        trialIdxPerRipple   = nan(size(rippleTimepoints, 1), 1); % trial index associated with each ripple
        trialPhasePerRipple = strings(size(rippleTimepoints, 1), 1); % trial phase associated with each ripple (with/without movement)
        objNamePerRipple    = nan(size(rippleTimepoints, 1), 1); % object (0-7) associated with each ripple
        objLocPerRipple     = nan(size(rippleTimepoints, 1), 2); % correct object location associated with each ripple
        respLocPerRipple    = nan(size(rippleTimepoints, 1), 2); % response location associated with each ripple
        bB4LearnPerRipple   = nan(size(rippleTimepoints, 1), 1); % learning status associated with each ripple
        for iR = 1:size(rippleTimepoints, 1)
            
            % timepoint of this ripple
            thisRippleTimepoint = rippleTimepoints(iR, 1);
            
            % assign information to this ripple
            for iTrial = 1:size(trialInfo, 1)
                
                % flag whether the ripple was contained in this trial
                bWithinTrial    = false;
                
                % test whether this ripple belongs to one of the trial
                % phases: cue, retrieval, feedback, reencoding, or nextITI
                if thisRippleTimepoint >= trialInfo.Cue(iTrial, 1) && thisRippleTimepoint < trialInfo.Retrieval(iTrial, 1) % cue
                    trialPhasePerRipple(iR, 1)  = 'C';
                    bWithinTrial                = true;
                
                elseif thisRippleTimepoint >= trialInfo.Retrieval(iTrial, 1) && thisRippleTimepoint < trialInfo.Feedback(iTrial, 1) % retrieval
                    trialPhasePerRipple(iR, 1)  = 'R';
                    bWithinTrial                = true;
                    
                    % movement versus non-movement periods
                    if any(thisRippleTimepoint >= movementStartEnd(:, 1) & thisRippleTimepoint <= movementStartEnd(:, 2))
                        trialPhasePerRipple(iR, 1)  = 'R_Movement';
                    else
                        trialPhasePerRipple(iR, 1)  = 'R_NoMovement';
                    end
                    
                elseif thisRippleTimepoint >= trialInfo.Feedback(iTrial, 1) && thisRippleTimepoint < trialInfo.Reencoding(iTrial, 1) % feedback
                    trialPhasePerRipple(iR, 1)  = 'F';
                    bWithinTrial                = true;
                    
                    % movement versus non-movement periods
                    if any(thisRippleTimepoint >= movementStartEnd(:, 1) & thisRippleTimepoint <= movementStartEnd(:, 2))
                        trialPhasePerRipple(iR, 1)  = 'F_Movement';
                    else
                        trialPhasePerRipple(iR, 1)  = 'F_NoMovement';
                    end
                
                elseif thisRippleTimepoint >= trialInfo.Reencoding(iTrial, 1) && thisRippleTimepoint < trialInfo.Grab(iTrial, 1) % re-encoding
                    trialPhasePerRipple(iR, 1)  = 'E';
                    bWithinTrial                = true;
                    
                    % movement versus non-movement periods
                    if any(thisRippleTimepoint >= movementStartEnd(:, 1) & thisRippleTimepoint <= movementStartEnd(:, 2))
                        trialPhasePerRipple(iR, 1)  = 'E_Movement';
                    else
                        trialPhasePerRipple(iR, 1)  = 'E_NoMovement';
                    end
                    
                elseif iTrial < size(trialInfo, 1) && thisRippleTimepoint >= trialInfo.Grab(iTrial, 1) && thisRippleTimepoint < trialInfo.Cue(iTrial + 1, 1) % next-ITI
                    trialPhasePerRipple(iR, 1)  = 'NI';
                    bWithinTrial              	= true;
                end
                
                % if the ripple was part of any period of this trial, store
                % various information
                if bWithinTrial == true
                    trialIdxPerRipple(iR, 1)    = trialInfo.TrialIdx(iTrial, 1); % trial index
                    objNamePerRipple(iR, 1)     = trialInfo.Object(iTrial, 1); % object name
                    objLocPerRipple(iR, :)      = [trialInfo.xCorrect(iTrial, 1), trialInfo.yCorrect(iTrial, 1)]; % object location
                    respLocPerRipple(iR, :)     = [trialInfo.xResponse(iTrial, 1), trialInfo.yResponse(iTrial, 1)]; % response location
                    bB4LearnPerRipple(iR, 1)    = trialInfo.bB4Learn(iTrial, 1); % before learning
                end
            end
        end
        
        % cutoff for early vs. late ripples
        halfCutoff      = median(trialIdxPerRipple, 'omitnan');
                
        % collect the information from all ripple channels
        rippleData(iRippleChan).rippleChanName      = rippleChans(iRippleChan).name; % name of this ripple channel
        rippleData(iRippleChan).rippleTimepoints    = rippleTimepoints; % ripple timepoints
        rippleData(iRippleChan).trialIdxPerRipple   = trialIdxPerRipple; % trial index assigned to each ripple
        rippleData(iRippleChan).bHalf1PerRipple     = trialIdxPerRipple < halfCutoff; % early (1) vs late (0) ripples
        rippleData(iRippleChan).trialPhasePerRipple = trialPhasePerRipple; % trial phase assigned to each ripple
        rippleData(iRippleChan).objNamePerRipple    = objNamePerRipple; % object name assigned to each ripple
        rippleData(iRippleChan).objLocPerRipple     = objLocPerRipple; % correct object location associated with each ripple
        rippleData(iRippleChan).respLocPerRipple    = respLocPerRipple; % response location associated with each ripple
        rippleData(iRippleChan).bB4LearnPerRipple   = bB4LearnPerRipple; % learning status associated with each ripple
    end
    
    %% single-neuron activity during the ripples
    
    % microwire channels and meta information
    spikeChans  = LK_GetChans(strcat(paths.spikes, subjects{iSub}, '\chan*'));
    subjectdata = load(fullfile(paths.subjects, subjects{iSub}, param.subject.fileName));
    
    % loop through microwire channels
    for iSpikeChan = 1:size(spikeChans, 1)
        
        % report
        fprintf('\tSpike channel: %s.\n', spikeChans(iSpikeChan).name);
        
        % channel index
        spikeChanIdx    = str2double(replace(spikeChans(iSpikeChan).name, lettersPattern, ''));
        
        % wave-clus files with spike times in macrotime
        c4aFile = fullfile(spikeChans(iSpikeChan).folder, spikeChans(iSpikeChan).name, param.micro.c4aFileName); % cluster-for-analysis file
        ccmFile = fullfile(spikeChans(iSpikeChan).folder, spikeChans(iSpikeChan).name, param.micro.ccmFileName); % cluster-class, microtime (msec), macrotime (sec)
        
        % load decision which clusters to analyze
        c4a = load(c4aFile);
        
        % if the decision file is empty, skip this channel
        if isempty(c4a.Cluster4Analysis)
            fprintf('\t... skipping this microwire channel, because the decision file is empty.\n');
            continue;
        end
        
        % load the spike-time information
        ccm                 = load(ccmFile); % cluster class, microtime (MSEC), macrotime (SEC)
        
        % brain region and hemisphere of this wire
        logIdx              = any(cell2mat(subjectdata.subjectdata.micro2macro(:, 1)) == spikeChanIdx, 2);
        spikeChanRegion     = subjectdata.subjectdata.micro2macro{logIdx, 3}; % brain region
        spikeChanHemisphere = subjectdata.subjectdata.micro2macro{logIdx, 4}; % hemisphere
        
        %% loop through clusters
        for iClus = 1:max(ccm.cluster_class_macrotime(:, 1))
            
            % consider visual inspection
            if ~strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'yes')
                exByInspection  = cat(1, exByInspection, [iSub, spikeChanIdx, iClus]);
                continue;
            end
            
            % data from this cluster
            thisCCM             = ccm.cluster_class_macrotime(ccm.cluster_class_macrotime(:, 1) == iClus, :); % cluster class, microtime (MSEC), macrotime (SEC)
            
            %% firing rates around ripples
            % for each combination of a given neuron and a given ripple
            % channel
            
            % loop through ripple channels
            for iRippleChan = 1:size(rippleData, 1)
                
                % macro-timepoints of all ripples of this channel
                rippleTimepoints   	= rippleData(iRippleChan).rippleTimepoints;
                
                % skip if there are no ripples on this ripple channel
                if size(rippleTimepoints, 1) < 1
                    fprintf('\tSkipping === %s (%s) === because there are no ripples on this ripple channel.\n', subjects{iSub}, rippleData(iRippleChan).rippleChanName);
                    exByNoRipples   = cat(1, exByNoRipples, [iSub, spikeChanIdx, iClus, iRippleChan]);
                    continue;
                end
                
                % firing rates around the ripple peak
                rippleFR    = nan(size(rippleTimepoints, 1), numel(param.ripple.timeCenters)); % ripples x time-centers
                for iR = 1:size(rippleTimepoints, 1)
                    thisTimeEdges   = rippleTimepoints(iR, 1) + param.ripple.timeEdges; % time bins around this ripple
                    rippleFR(iR, :) = histcounts(thisCCM(:, 3), thisTimeEdges) ./ round(diff(thisTimeEdges, [], 2), 5); % firing rates
                end
                
                %% results for this neuron-ripple channel-combination
                
                % unit information
                thisRes                     = [];
                thisRes.subject             = subjects{iSub};
                thisRes.idx                 = [iSub, spikeChanIdx, iClus];
                thisRes.spikeChanRegion     = spikeChanRegion;
                thisRes.spikeChanHemisphere = spikeChanHemisphere;
                
                % ripple information
                thisRes.rippleChanName      = rippleData(iRippleChan).rippleChanName;
                thisRes.rippleTimepoints    = rippleData(iRippleChan).rippleTimepoints;
                thisRes.trialIdxPerRipple   = rippleData(iRippleChan).trialIdxPerRipple;
                thisRes.bHalf1PerRipple     = rippleData(iRippleChan).bHalf1PerRipple;
                thisRes.trialPhasePerRipple = rippleData(iRippleChan).trialPhasePerRipple;
                thisRes.objNamePerRipple    = rippleData(iRippleChan).objNamePerRipple;
                thisRes.objLocPerRipple     = rippleData(iRippleChan).objLocPerRipple;
                thisRes.respLocPerRipple    = rippleData(iRippleChan).respLocPerRipple;
                thisRes.bB4LearnPerRipple   = rippleData(iRippleChan).bB4LearnPerRipple;
                
                % firing rates during ripples
                thisRes.rippleFR            = rippleFR;
                
                % behavior
                thisRes.trialInfo           = trialInfo;
                
                % assemble across different units
                allRes                      = cat(1, allRes, thisRes);
                
                %% close open figures
                close all;
            end
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load results

% load previous results
fprintf('\n==== Results.\n');
r = load(strcat(paths.save, 'results.mat'), 'allRes', 'param', 'paths', 'subjects', 'objRes', 'objCells', 'placeRes', 'placeCells');

%% general information

% cellular indices
allIdx      = cell2mat({r.allRes.idx}'); % all unit indices of all unit-ripple channel-combinations
uniqueIdx   = unique(allIdx, 'rows', 'stable'); % unique unit indices of all unit-ripple channel-combinations

% labels of the ripple channels
allRippleChanName       = {r.allRes.rippleChanName}'; % ripple channel name for each cellular index
uniqueRippleChanName    = unique(allRippleChanName, 'stable'); % unique labels of the ripple channels

% brain regions of the cells
allSpikeChanRegion      = {r.allRes.spikeChanRegion}';
uniqueSpikeChanRegion   = unique(allSpikeChanRegion); % unique brain regions of the spike channels

% hemisphere of all cells
allSpikeChanHemisphere  = {r.allRes.spikeChanHemisphere}';

% information whether a neuron and its associated ripples were recorded
% from the same hemisphere
allBSameHemiSpikeRipple = false(size(r.allRes, 1), 1);
for iRes = 1:size(r.allRes, 1)
    if (strcmp(allSpikeChanHemisphere{iRes}, 'left') && contains(allRippleChanName{iRes}, 'L')) || ...
            (strcmp(allSpikeChanHemisphere{iRes}, 'right') && contains(allRippleChanName{iRes}, 'R'))
        allBSameHemiSpikeRipple(iRes, 1)    = true;
    end
end

% report
fprintf('Number of unit-ripple combinations: %d.\n', size(allIdx, 1));
fprintf('Number of unique units: %d.\n', size(uniqueIdx, 1));
fprintf('Number of units recorded from the same hemisphere as their hippocampal ripples: %d.\n', sum(allBSameHemiSpikeRipple));
fprintf('Number of units not recorded from the same hemisphere as their hippocampal ripples: %d.\n', sum(~allBSameHemiSpikeRipple));
fprintf('Names of the ripple channels: %s.\n', uniqueRippleChanName{:});

%% assign cell types to neuron-ripple channel-combinations

% for each cell from the unit-ripple analysis, identify whether it is an
% object cell or a place cell
allCellTypes    = repmat(struct(), size(allIdx, 1), 1);
for iIdx = 1:size(allIdx, 1)
    
    % check whether it is an object cell
    tmpObjIdx                       = all(allIdx(iIdx, :) == [r.objCells.cellClassification.session, r.objCells.cellClassification.channel, r.objCells.cellClassification.cluster], 2);
    allCellTypes(iIdx).objectCell   = r.objCells.cellClassification.cellType{tmpObjIdx};
    
    % check whether it is a place cell
    tmpPlaceIdx                     = all(allIdx(iIdx, :) == [r.placeCells.cellClassification.session, r.placeCells.cellClassification.channel, r.placeCells.cellClassification.cluster], 2);
    allCellTypes(iIdx).placeCell    = r.placeCells.cellClassification.cellType{tmpPlaceIdx};
    
    % sanity check whether the indexing for object cells and place cells is
    % identical
    if sum(tmpObjIdx == tmpPlaceIdx) ~= size(tmpObjIdx, 1)
        error('"tmpObjIdx" and "tmpPlaceIdx" are different.');
    end
end

%% cellular activity during the ripples, separately for different ripple-locked time windows (discretized: yes/no)

% cells x time-windows-during-ripples
allBActiveDuringRipples = cell(size(r.allRes, 1), size(r.param.coac.twoi, 1)); % in each cell, you store whether the neuron is active during each of the ripples
for iRes = 1:size(r.allRes, 1) % loop through cells
    for iTwoi = 1:size(r.param.coac.twoi, 1) % loop through time windows
        
        % time window of interest (relative to ripple peak)
        thisTwoi    = r.param.coac.twoi(iTwoi, :);
        bThisTwoi   = r.param.ripple.timeCenters >= min(thisTwoi) & r.param.ripple.timeCenters < max(thisTwoi);
        
        % check whether this cell was active during this time window
        thisTwoiFR                              = r.allRes(iRes).rippleFR(:, bThisTwoi); % ripples x timepoints
        allBActiveDuringRipples{iRes, iTwoi}    = any(thisTwoiFR, 2); % ripples x 1
    end
end

%% coactivity between object cells and place cells during ripples

% compute and evaluate coactivity
cfg                             = [];
cfg.r                           = r;
cfg.allIdx                      = allIdx; % unit indices
cfg.allCellTypes                = allCellTypes; % cell types (object, place)
cfg.allSpikeChanRegion          = allSpikeChanRegion; % spike-channel regions
cfg.allSpikeChanHemisphere      = allSpikeChanHemisphere; % spike-channel hemispheres
cfg.allRippleChanName           = allRippleChanName; % ripple-channel names
cfg.allBSameHemiSpikeRipple     = allBSameHemiSpikeRipple; % whether spikes and ripples were recorded from the same hemisphere
cfg.allBActiveDuringRipples     = allBActiveDuringRipples; % cellular activity during ripples
outCoac                         = LK_CoactivityDuringRipples_Function_20230317(cfg);

%% coactivity between object cells and place cells during ripples, first sessions only

% restrict data to first sessions
allSubject                      = {r.allRes.subject}';
bMaskSess                       = cellfun(@(x) ~contains(x(end), {'b', 'c'}), allSubject);

% compute and evaluate coactivity
cfg                             = [];
cfg.r                           = r;
cfg.r.paths.save                = strcat(r.paths.save, 'Sess1\');
cfg.r.paths.pics                = strcat(cfg.r.paths.save, 'Pics\');
cfg.r.allRes                    = r.allRes(bMaskSess);
cfg.allIdx                      = allIdx(bMaskSess, :);
cfg.allCellTypes                = allCellTypes(bMaskSess);
cfg.allSpikeChanRegion          = allSpikeChanRegion(bMaskSess);
cfg.allSpikeChanHemisphere      = allSpikeChanHemisphere(bMaskSess);
cfg.allRippleChanName           = allRippleChanName(bMaskSess);
cfg.allBSameHemiSpikeRipple     = allBSameHemiSpikeRipple(bMaskSess);
cfg.allBActiveDuringRipples     = allBActiveDuringRipples(bMaskSess, :); % cellular activity during ripples
outCoacSess1                 	= LK_CoactivityDuringRipples_Function_20230317(cfg);

%% coactivity between object cells and place cells during ripples, Pearson correlation

% compute and evaluate coactivity
cfg                             = [];
cfg.r                           = r;
cfg.r.param.coac.type           = 'r';
cfg.r.param.coac.clim.ret       = [-0.2, 0.2]; % c-limits for 2D-coactivity plots during retrieval
cfg.r.param.coac.clim.enc       = [-0.3, 0.3]; % c-limits for 2D-coactivity plots during encoding
cfg.r.paths.save                = strrep(r.paths.save, '_z_', '_r_');
cfg.r.paths.pics                = strcat(cfg.r.paths.save, 'Pics\');
cfg.allIdx                      = allIdx;
cfg.allCellTypes                = allCellTypes;
cfg.allSpikeChanRegion          = allSpikeChanRegion;
cfg.allSpikeChanHemisphere      = allSpikeChanHemisphere;
cfg.allRippleChanName           = allRippleChanName;
cfg.allBSameHemiSpikeRipple     = allBSameHemiSpikeRipple;
cfg.allBActiveDuringRipples     = allBActiveDuringRipples; % cellular activity during ripples
outCoacPearson                  = LK_CoactivityDuringRipples_Function_20230317(cfg);
