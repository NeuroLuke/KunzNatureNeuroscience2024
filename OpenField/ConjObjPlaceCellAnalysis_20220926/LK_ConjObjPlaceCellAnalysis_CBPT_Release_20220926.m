%==========================================================================
% This script analyzes conjunctive object-place cells using a procedure
% that involves cluster-based permutation testing.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                	= [];
param.numSurrogates    	= 1001; % number of surrogates
param.myRNG            	= 666;
param.subjectFileName   = 'subjectdata_20210616.mat';
param.clusterName     	= 'LK_Cluster4Analysis_20210521';
%- for behavior
beh                     = [];
beh.newTimeRes         	= 0.1;
beh.naviCutoff        	= 0.001;
beh.naviSmooth      	= 2;
beh.rotCutoff         	= 0.001;
beh.rotSmooth         	= 2;
beh.minNumObs           = 2; % minimum number of observations per bin
%- for location
place                	= [];
place.res             	= 25; % spatial resolution
place.xEdges          	= linspace(-5000, 5000, place.res + 1);
place.xCenters         	= movmean(place.xEdges, 2, 'endpoints', 'discard');
place.yEdges          	= linspace(-5000, 5000, place.res + 1);
place.yCenters        	= movmean(place.yEdges, 2, 'endpoints', 'discard');
place.idxTemplate       = flipud(reshape(1:place.res ^ 2, place.res, place.res)); % template for location indexing
place.facThreshPF       = 0.75;
place.smoothSettings    = {'gaussian', 5, 1.5};

% paths
paths         	= [];
paths.subjects 	= 'E:\OpenField\SubjectData_20210616\';
paths.spikes  	= 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.beh     	= 'E:\OpenField\BehPreprocessing_20210707\';
paths.cueObj    = 'E:\OpenField\CueObjAnalysis_20210812\20211116\'; % object cells
paths.place     = 'E:\OpenField\PlaceCellAnalysis_20210511\20211116_Loc25x25_facThreshPF0.75\';
paths.save    	= strcat('E:\OpenField\ConjObjPlaceCellAnalysis_20220926\20220926_', ...
    'Loc', num2str(numel(place.xCenters)), 'x', num2str(numel(place.yCenters)), ...
    '_facThreshPF', num2str(place.facThreshPF), '\');
mkdir(paths.save);

% subjects
subjects     	= load(strcat(paths.subjects, 'subjects.mat'));
subjects    	= subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% save settings
save(strcat(paths.save, 'settings'));

%% preallocations

% main results for each cell
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByVisInsp 	= []; % excluded by visual inspection

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % skip if it is not a microwire subject
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        continue;
    end
    
    % reset rng for reproducibility
    rng(param.myRNG);
    
    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
        
    %% behavioral processing: masking for movements and rotations
    
    % report
    fprintf('\nProcessing of behavioral data (using microtime).\n');
    
    % load trial and behavioral information in microtime
    trialInfo       = load(strcat(paths.beh, subjects{iSub}, '\trialInfo.mat'), 'trialInfoMicrotime');
    trialInfo       = trialInfo.trialInfoMicrotime;
    behInfo         = load(strcat(paths.beh, subjects{iSub}, '\behInfoRes', num2str(1/beh.newTimeRes), 'Hz.mat'), 'behInfoMicrotimeRes'); % resampled data
    behInfo         = behInfo.behInfoMicrotimeRes;
    
    % restrict behavior to time period between first and last cue
    logIdx          = behInfo.time >= min(trialInfo.Cue) & behInfo.time <= max(trialInfo.Cue);
    behInfo         = behInfo(logIdx, :);
    fprintf('The duration between the first and the last cue is %.1f s.\n', range(trialInfo.Cue));
    fprintf('Percentage of "behInfo" between the first and the last cue: %.3f%%.\n', 100 * sum(logIdx) / size(logIdx, 1));
    
    % behavioral mask
    bMove           = smoothdata(behInfo.speed > beh.naviCutoff, 1, 'movmean', beh.naviSmooth / beh.newTimeRes + 1) > 0;
    bRotate         = smoothdata(behInfo.yawSpeed > beh.rotCutoff, 1, 'movmean', beh.rotSmooth / beh.newTimeRes + 1) > 0;
    bActive         = behInfo.trialPhase > 2;
    bMask4Analysis  = (bMove | bRotate) & bActive;
    
    % report
    fprintf('Percentage of the "active time" that is included into the analysis after masking for movements and rotations: %.3f%%.\n', ...
        sum(behInfo.durations(bMask4Analysis)) / sum(behInfo.durations(bActive)) * 100);
    
    %% behavioral processing: predictors
    
    % binned x/y locations
    xBins   = discretize(behInfo.x, place.xEdges);
    yBins 	= discretize(behInfo.y, place.yEdges);
    xyBins 	= (xBins - 1) .* numel(place.yCenters) + yBins;
    
    % object identity for each time bin
    objBins = nan(size(behInfo, 1), 1);
    for iTrial = 1:max(trialInfo.TrialIdx)
        bThisTrial          = behInfo.trialIdx == iTrial;
        objBins(bThisTrial) = trialInfo.Object(iTrial);
    end
    
    % sanity-check
    if sum(isnan(objBins)) > 0
        error('There is a NaN in the predictor "objBins".');
    end
    
    %% behavioral processing: sufficient sampling across the entire time
    
    % identify timepoints that are part of sufficient sampling
    cfg                 = [];
    cfg.data            = xyBins(bMask4Analysis); % time x predictors
    cfg.factorNames     = {'Loc'}; % predictor names
    cfg.minNumObs       = beh.minNumObs;
    bAddMask4Analysis   = LK_MaskBehavioralSampling(cfg); % additional mask for analysis
    
    % extend the behavioral mask to ensure sufficient sampling
    bMask4Analysis(bMask4Analysis)  = bAddMask4Analysis;
    
    % report
    fprintf('Percentage of the "active time" that is included in the analysis after ensuring sufficient sampling: %.3f%%.\n', ...
        sum(behInfo.durations(bMask4Analysis)) / sum(behInfo.durations(bActive)) * 100);
    
    %% channels with microwire data
    
    % available channels
    chans = LK_GetChans(fullfile(paths.spikes, subjects{iSub}, 'chan*'));
    
    % report
    fprintf('\nProcessing of neural data (%d channels from "%s" to "%s").\n', size(chans, 1), chans(1).name, chans(end).name);
    
    %% loop through wires
    for iChan = 1:size(chans, 1)
        
        %% index and brain region of this microwire channel
        
        % report
        fprintf('\tChannel: %s.\n', chans(iChan).name);
        
        % original channel index
        chanIdx     = str2double(replace(chans(iChan).name, lettersPattern, ''));
        
        % subject information
        s           = load(strcat(paths.subjects, subjects{iSub}, filesep, param.subjectFileName));
        logIdx  	= any(cell2mat(s.subjectdata.micro2macro(:, 1)) == chanIdx, 2);
        chanRegion  = s.subjectdata.micro2macro{logIdx, 3}; % brain region
        
        %% wave-clus output
        
        % load wave-clus output
        wcFile  = fullfile(chans(iChan).folder, chans(iChan).name, 'times_datacut.mat');
        if exist(wcFile, 'file') > 0
            times_datacut   = load(wcFile);
        else
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC          = cat(1, exByWC, [iSub, chanIdx]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on inspection)
        c4a = load(strcat(chans(iChan).folder, filesep, chans(iChan).name, filesep, param.clusterName, '.mat'));
        
        %% loop through clusters
        for iClus = 1:max(times_datacut.cluster_class(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            tic
            
            % continue if you decided that this cluster is insufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'no')
                fprintf('\t\t- You decided not to analyze this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, chanIdx, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster = times_datacut.cluster_class(times_datacut.cluster_class(:, 1) == iClus, :); % cluster-class, MICROTIME (MSEC)
            thisSpike   = times_datacut.spikes(times_datacut.cluster_class(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            % convert spike times to seconds
            thisCluster(:, 2)   = thisCluster(:, 2) ./ 1000;
            
            %% firing rates
            
            % number of spikes and firing rate per timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 2), behInfo.time), 0]); % add 0 spikes for the last time bin
            FR          = numSpikes ./ behInfo.durations;
            
            %% place tuning, separately for the eight different objects
            
            % loop through the objects
            for iObj = min(objBins):max(objBins)
                
                % reduce the behavior and the firing rates to the analysis
                % time window
                rxyBins     = xyBins(bMask4Analysis & objBins == iObj);
                robjBins    = objBins(bMask4Analysis & objBins == iObj);
                rFR         = FR(bMask4Analysis & objBins == iObj);
                rNumSpikes  = numSpikes(bMask4Analysis & objBins == iObj);
                rDurations  = behInfo.durations(bMask4Analysis & objBins == iObj);
                
                % surrogate firing rates (obtained by circularly shifting
                % the empirical firing rates)
                randShifts  = datasample(1:numel(rFR), param.numSurrogates, 'replace', true);
                rFRSurro    = nan(size(rFR, 1), param.numSurrogates);
                for iSurro = 1:param.numSurrogates
                    rFRSurro(:, iSurro) = circshift(rFR, randShifts(iSurro));
                end
                
                %% empirical place tuning
                
                % identify one place field per cell and its strength
                cfg         = [];
                cfg.FR      = rFR; % firing rates
                cfg.xyBins  = rxyBins; % place information
                cfg.place   = place; % place settings
                placeTest   = LK_PlaceCellAnalysis_placeFieldTest_20210814(cfg);
                
                %% surrogate place tuning
                
                % loop through surrogates
                tPlaceSurro     = nan(param.numSurrogates, 1);
                strengthPFSurro = nan(param.numSurrogates, 1);
                for iSurro = 1:param.numSurrogates
                    
                    % identify one surrogate place field per cell and its
                    % strength
                    cfg             = [];
                    cfg.FR          = rFRSurro(:, iSurro);
                    cfg.xyBins      = rxyBins;
                    cfg.place       = place;
                    placeTestSurro  = LK_PlaceCellAnalysis_placeFieldTest_20210814(cfg);
                    
                    % collect information across surrogate rounds
                    tPlaceSurro(iSurro, 1)      = placeTestSurro.tPlace;
                    strengthPFSurro(iSurro, 1)  = placeTestSurro.strengthPF;
                end
                
                %% collect information for this unit and object
                
                % basics
                unitRes                 = [];
                unitRes.idx             = [iSub, chanIdx, iClus, iObj];
                unitRes.subjectName     = subjects{iSub, 1};
                unitRes.chanRegion      = chanRegion;
                unitRes.nspk            = sum(rNumSpikes);
                unitRes.totDur          = sum(rDurations);
                
                % place information
                unitRes.FRMap           = placeTest.FRMap;
                unitRes.smFRMap         = placeTest.smFRMap;
                unitRes.threshPF        = placeTest.threshPF;
                unitRes.PF              = placeTest.PF;
                unitRes.strengthPF      = placeTest.strengthPF;
                unitRes.pPlace          = placeTest.pPlace;
                unitRes.tPlace          = placeTest.tPlace;
                unitRes.strengthPFSurro = strengthPFSurro;
                unitRes.tPlaceSurro     = tPlaceSurro;
                                
                % collapse across units
                allRes                  = cat(1, allRes, unitRes);
                
                %% figure: place tuning
                
                % p-value
                thisPT              = 1 - sum(placeTest.tPlace > tPlaceSurro) / sum(~isnan(tPlaceSurro));
                
                % data for figure
                cfg                 = [];
                cfg.visible         = 'off';
                cfg.place           = place;
                cfg.FR              = placeTest.smFRMap;
                cfg.PF              = placeTest.PF;
                cfg.path            = [behInfo.x, behInfo.y]; % entire navigation path
                cfg.bPathOfInterest = objBins == iObj;
                cfg.thisSpike       = thisSpike; % all spikes
                cfg.sr              = times_datacut.par.sr; % sampling rate
                cfg.nspk            = sum(rNumSpikes); % number of spikes used for this analysis (object-specific)
                cfg.figTitle        = {'Place tuning', ['\itP\rm = ', num2str(thisPT, '%.3f'), ' (', chanRegion, ')']};
                if thisPT < 0.001
                    cfg.figTitle    = {'Place tuning', ['\itP\rm < 0.001 (', chanRegion, ')']};
                end
                cfg.idx             = [iSub, chanIdx, iClus, iObj + 1]; % state the object index, not its "name"
                cfg.tPlace          = placeTest.tPlace;
                cfg.tPlaceSurro     = tPlaceSurro;
                f                   = LK_PlotPlaceCell_20220330(cfg);
                
                % save figure
                set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
                print(f, strcat(paths.save, subjects{iSub}, '_', chans(iChan).name, '_clus', num2str(iClus), '_obj', num2str(iObj), '_loc_20220330'), '-dtiff', '-r450');
            end
            
            %% close all open figures
            close all;
            toc
        end
    end
end

%% results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r = load(strcat(paths.save, 'results.mat'));

% remove fieldtrip from the search path
rmpath(genpath('E:\fieldtrip\'));
clc;

%% general

% general information about all cell x object combinations
allIdx          = cell2mat({r.allRes.idx}');
allSubjectName  = {r.allRes.subjectName}';
fprintf('Number of cell x object combinations: %d.\n', size(allIdx, 1));

% brain regions of all cell x object combinations
allRegions      = {r.allRes.chanRegion}';

%% conjunctive cells

% define alpha level
uncorrAlpha     = 0.05; % uncorrected alpha level
corrAlpha       = uncorrAlpha / 8; % correction for eight comparisons

% identify conjunctive cells based on the object-wise t-statistic
allTPlace       = cell2mat({r.allRes.tPlace}');
allTPlaceSurro  = cell2mat({r.allRes.tPlaceSurro})';
allTPlaceRank   = sum(allTPlace > allTPlaceSurro, 2) ./ sum(~isnan(allTPlaceSurro), 2);

% loop through cells and identify information at the cell level
cellIdx         = unique(allIdx(:, 1:3), 'rows', 'stable');
cellRegions     = cell(size(cellIdx, 1), 1);
cellBConjCell   = false(size(cellIdx, 1), 1);
for iCell = 1:size(cellIdx, 1)
    
    % logical indexing of this cell
    bThisCell               = all(allIdx(:, 1:3) == cellIdx(iCell, :), 2);
    
    % cell region
    cellRegions(iCell, 1)   = unique(allRegions(bThisCell));
    
    % identification of conjunctive cells: is the place tuning significant
    % for exactly one object?
    cellBConjCell(iCell, 1) = sum(allTPlaceRank(bThisCell) > (1 - corrAlpha)) == 1;
end

% report
fprintf('Number of conjunctive cells based on a significant t-statistic: %d of %d (%.3f%%, P = %.3f).\n', sum(cellBConjCell), numel(cellBConjCell), ...
    100 * sum(cellBConjCell) / numel(cellBConjCell), myBinomTest(sum(cellBConjCell), numel(cellBConjCell), 0.05));

%% conjunctive cells per brain region

% report
fprintf('\nAnalysis of conjunctive cells in different brain regions.\n');

% unique regions
uniqueRegions   = unique(cellRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% loop through the different regions
for iR = 1:size(uniqueRegions, 1)
    
    % cells from this region
    bThisRegion             = strcmp(cellRegions, uniqueRegions{iR});
    numCellsPerReg(iR, 1)   = sum(bThisRegion);
    fprintf('Region: %s. Number of cells: %d.\n', uniqueRegions{iR}, sum(bThisRegion));
    
    % binomial test
    pBinom  = myBinomTest(sum(cellBConjCell & bThisRegion), sum(bThisRegion), 0.05);
    
    % report
    fprintf('\tNumber of conjunctive cells: %d (%.3f%%, P = %.3f).\n', sum(cellBConjCell & bThisRegion), ...
        100 * sum(cellBConjCell & bThisRegion) / sum(bThisRegion), pBinom);
end

% restrict unique regions to regions with enough units
uniqueRegions       = uniqueRegions(numCellsPerReg > 30);

% create figure for percentage of place cells per brain region
cfg                         = [];
cfg.uniqueRegions           = uniqueRegions;
cfg.allRegions              = cellRegions;
cfg.bCell                   = cellBConjCell;
cfg.ylabel                  = 'Conjunctive cells (%)';
cfg.BonferroniCorrection    = true;
f                           = LK_PlotCellsPerRegion(cfg);
LK_print(f, strcat(paths.save, 'conjunctiveCellsPerRegion_20230925'), '-dtiff', '-r300');

%% load previous place- and object-cell results

% report
fprintf('\nLoading information about object cells and place cells.\n');

% object cells and object results
objRes      = load(strcat(r.paths.cueObj, 'results.mat'), 'allRes');
objCells    = load(strcat(r.paths.cueObj, 'cellClassification.mat'));
bObjCell    = strcmp(objCells.cellClassification.cellType, 'ObjectCell');
fprintf('Number of object cells: %d (out of %d).\n', sum(bObjCell), size(bObjCell, 1));

% place cells
placeRes    = load(strcat(r.paths.place, 'results.mat'), 'allRes', 'place');
placeCells  = load(strcat(r.paths.place, 'cellClassification.mat'));
bPlaceCell  = strcmp(placeCells.cellClassification.cellType, 'PlaceCell');
fprintf('Number of place cells: %d (out of %d).\n', sum(bPlaceCell), size(bPlaceCell, 1));

%% overlap between conjunctive cells, object cells, and place cells

% report
fprintf('\nExamining the overlap between conjunctive cells, object cells and place cells.\n');

% overlap between conjunctive cells, object cells and place cells
fprintf('Number of cells that are conjunctive, object, and place cells: %d.\n', sum(cellBConjCell & bPlaceCell & bObjCell));

% figure: Venn diagram for the overlap between conjunctive, place, and
% object cells
A = [sum(bObjCell), sum(bPlaceCell), sum(cellBConjCell)];
I = [sum(bObjCell & bPlaceCell), sum(bObjCell & cellBConjCell), sum(bPlaceCell & cellBConjCell), sum(bObjCell & bPlaceCell & cellBConjCell)];
f = figure('units', 'centimeters', 'position', [2, 5, 12, 8]);
venn(A, I, 'FaceColor', {rgb('orange'), rgb('blue'), mean([rgb('orange'); rgb('blue')], 1)}, 'LineWidth', 2);
axis equal tight off;
% save figure
LK_print(f, strcat(r.paths.save, 'Venn_ObjCells_vs_PlaceCells_vs_ConjCells_20220927'), '-dtiff', '-r300');

%% similarity of place tuning between the different objects

% loop through cells
simPlaceTuning  = nan(size(cellIdx, 1), 1);
for iCell = 1:size(cellIdx, 1)
    
    % logical indexing of this cell
    bThisCell   = all(allIdx(:, 1:3) == cellIdx(iCell, :), 2);
    
    % place tuning in this cell, separately for the different objects
    thisMaps    = {r.allRes(bThisCell).smFRMap};
    
    % pair-wise correlation between the firing-rate maps
    rho         = nan(numel(thisMaps), numel(thisMaps));
    for i = 1:numel(thisMaps)
        for j = 1:numel(thisMaps)
            rho(i, j)   = corr(thisMaps{i}(:), thisMaps{j}(:), 'rows', 'complete', 'type', 'pearson');
        end
    end
    
    % average similarity
    bOI                         = tril(ones(size(rho)), -1) ~= 0;
    simPlaceTuning(iCell, 1)    = mean(rho(bOI), 1, 'omitnan');
end

%% figure: bar plots showing the place-tuning similarity across objects

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.6, 1.4, 4.3, 4.3]);
hold on;

% loop through cell groups
cellGroups  = {'conjunctiveCells', 'Conj.'; 'objectCells', 'Object'; 'placeCells', 'Place'};
for iCG = 1:size(cellGroups, 1)

    % select data
    if strcmp(cellGroups{iCG}, 'conjunctiveCells')
        thisData    = simPlaceTuning(cellBConjCell);
    elseif strcmp(cellGroups{iCG}, 'objectCells')
        thisData    = simPlaceTuning(bObjCell);
    elseif strcmp(cellGroups{iCG}, 'placeCells')
        thisData    = simPlaceTuning(bPlaceCell);
    end
    
    % mean and standard error
    m   = mean(thisData, 1, 'omitnan');
    ste = LK_ste(thisData);
    
    % statistical test against 0
    [~, p, ~, stats]    = ttest(thisData, 0);
    fprintf('Stable place tuning across objects in %s? t(%d) = %.3f, p = %.3f.\n', cellGroups{iCG}, stats.df, stats.tstat, p);
    
    % plot bar
    b1 = bar(iCG, m, 'FaceColor', [0.9, 0.9, 0.9]);
    s1 = plot([iCG, iCG], [m - ste, m + ste], 'Color', [0, 0, 0], 'LineWidth', 2);
    % indicate significance
    if p < 0.001
        t = text(b1.XData, max(s1.YData), '***');
    elseif p < 0.01
        t = text(b1.XData, max(s1.YData), '**');
    elseif p < 0.05
        t = text(b1.XData, max(s1.YData), '*');
    end
    % enhance text
    if p < 0.05
        set(t, 'fontunits', 'centimeters', 'fontsize', 0.4, ...
            'verticalalignment', 'baseline', 'horizontalalignment', 'center');
    end
end

% stats: is place-tuning similarity higher for place than for conjunctive
% cells?
[~, p, ~, stats]    = ttest2(simPlaceTuning(bPlaceCell), simPlaceTuning(cellBConjCell));
fprintf('More stable tuning across objects in place cells than in conjunctive cells? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
plot([1, 3], [0.0575, 0.0575], 'Color', [0, 0, 0], 'LineWidth', 1);
if p < 0.01
    text(2, 0.057, '**', 'fontunits', 'centimeters', 'fontsize', 0.4, 'verticalalignment', 'baseline', 'horizontalalignment', 'center');
end

% stats: is place-tuning similarity higher for place than for object cells?
[~, p, ~, stats]    = ttest2(simPlaceTuning(bPlaceCell), simPlaceTuning(bObjCell));
fprintf('More stable tuning across objects in place cells than in objects cells? t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
plot([2, 3], [0.0525, 0.0525], 'Color', [0, 0, 0], 'LineWidth', 1);
if p < 0.05
    text(2.5, 0.052, '*', 'fontunits', 'centimeters', 'fontsize', 0.4, 'verticalalignment', 'baseline', 'horizontalalignment', 'center');
end
hold off;

% enhance axes
set(gca, ...
    'xlim', [0.2, size(cellGroups, 1) + 0.8], 'ylim', [0, 0.06], ...
    'xtick', 1:size(cellGroups, 1), 'xticklabel', cellGroups(:, 2), 'tickdir', 'out');
xl = xlabel('cells');
yl = ylabel('Place-tuning similarity (\itr\rm)');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);

% save figure
LK_print(f, strcat(r.paths.save, 'simPlaceTuning_ObjCells_vs_PlaceCells_vs_ConjCells_20220927'), '-dtiff', '-r300');
