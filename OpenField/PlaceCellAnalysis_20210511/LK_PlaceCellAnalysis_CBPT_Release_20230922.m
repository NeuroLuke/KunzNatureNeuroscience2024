%==========================================================================
% This script analyzes place cells using a procedure that involves
% cluster-based permutation testing.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                	= [];
param.myRNG            	= 444;
param.numSurrogates    	= 1001; % number of surrogates
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
place.facThreshPF       = 0.75; % threshold for defining place fields
place.smoothSettings    = {'gaussian', 5, 1.5};

% paths
paths         	= [];
paths.subjects 	= 'E:\OpenField\SubjectData_20210616\';
paths.spikes  	= 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.beh     	= 'E:\OpenField\BehPreprocessing_20210707\';
paths.cueObj    = 'E:\OpenField\CueObjAnalysis_20210812\20211116\'; % object cells
paths.save    	= strcat('E:\OpenField\PlaceCellAnalysis_20210511\20211116_', ...
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
    
    % subject-specific save path
    subjSavePath    = strcat(paths.save, subjects{iSub}, '\');
    mkdir(subjSavePath);
        
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
    
    %% behavioral processing: sufficient sampling
    
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
        chanRegion  = s.subjectdata.micro2macro{logIdx, 3}; % brain region of this channel
        
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
            
            % continue if you decided that this cluster is insufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'no')
                fprintf('\t\t- You decided not to analyze this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, chanIdx, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster = times_datacut.cluster_class(times_datacut.cluster_class(:, 1) == iClus, :); % cluster-class, microtime (msec)
            thisSpike   = times_datacut.spikes(times_datacut.cluster_class(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            % convert spike times to seconds
            thisCluster(:, 2)   = thisCluster(:, 2) ./ 1000;
            
            %% firing rates
            
            % number of spikes and firing rate per timepoint
            numSpikes   = transpose([histcounts(thisCluster(:, 2), behInfo.time), 0]); % 0 spikes for the last time bin
            FR          = numSpikes ./ behInfo.durations;
            
            % reduce the behavior and the firing rates to the analysis time
            % window
            rxyBins     = xyBins(bMask4Analysis);
            rFR         = FR(bMask4Analysis);
            rNumSpikes  = numSpikes(bMask4Analysis);
            rDurations  = behInfo.durations(bMask4Analysis);
            
            % surrogate firing rates (obtained by circularly shifting the
            % empirical firing rates)
            randShifts  = datasample(1:numel(rFR), param.numSurrogates, 'replace', false);
            rFRSurro    = nan(size(rFR, 1), param.numSurrogates);
            for iSurro = 1:param.numSurrogates
                rFRSurro(:, iSurro) = circshift(rFR, randShifts(iSurro));
            end
                        
            %% empirical place tuning
            
            % identify place tuning
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
                
                % identify surrogate place tuning
                cfg             = [];
                cfg.FR          = rFRSurro(:, iSurro);
                cfg.xyBins      = rxyBins;
                cfg.place       = place;
                placeTestSurro  = LK_PlaceCellAnalysis_placeFieldTest_20210814(cfg);
                
                % collect information across surrogate rounds
                tPlaceSurro(iSurro, 1)      = placeTestSurro.tPlace;
                strengthPFSurro(iSurro, 1)  = placeTestSurro.strengthPF;
            end
                        
            %% additional property: temporal stability (correlation between the firing-rate maps from the 1st vs. 2nd data half)
            
            % identify data from 1st half
            b1stHalf    = 1:numel(rFR) <= (numel(rFR) / 2);
            
            % groups
            groups    	= {'1stHalf'; '2ndHalf'};
            smFRMaps    = nan(numel(placeTest.smFRMap), size(groups, 1)); % e.g., size = 625 x 2
            for iG = 1:size(groups, 1)
                
                % select data
                if strcmp(groups{iG}, '1stHalf')
                    bTmpMask        = b1stHalf;
                elseif strcmp(groups{iG}, '2ndHalf')
                    bTmpMask        = ~b1stHalf;
                end
                
                % firing-rate map for this data half
                cfg             = [];
                cfg.FR          = rFR(bTmpMask);
                cfg.xyBins      = rxyBins(bTmpMask);
                cfg.place       = place;
                placeTestHalf   = LK_PlaceCellAnalysis_placeFieldTest_20210814(cfg);
                
                % collect across groups
                smFRMaps(:, iG) = reshape(flipud(placeTestHalf.smFRMap), numel(placeTestHalf.smFRMap), 1); % unfold the firing-rate map
            end
            
            % temporal stability as the Pearson correlation
            temporalStab    = corr(smFRMaps, 'rows', 'complete');
            temporalStab    = mean(temporalStab(tril(temporalStab, -1) ~= 0));
            
            %% collect information for this unit
            
            % basics
            unitRes                 = [];
            unitRes.idx             = [iSub, chanIdx, iClus];
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
            
            % additional properties
            unitRes.temporalStab    = temporalStab;
            
            % collapse across units
            allRes                  = cat(1, allRes, unitRes);
            
            %% figure: place tuning
            
            % p-value
            thisPT              = 1 - sum(placeTest.tPlace > tPlaceSurro) / sum(~isnan(tPlaceSurro));
            
            % data for figure
            dt                  = [];
            dt.visible          = 'off';
            dt.place            = place;
            dt.FR               = placeTest.smFRMap;
            dt.PF               = placeTest.PF;
            dt.path             = [behInfo.x, behInfo.y]; % navigation path
            dt.thisSpike        = thisSpike;
            dt.sr               = times_datacut.par.sr;
            dt.nspk             = sum(numSpikes(bMask4Analysis));
            dt.figTitle         = {'Place tuning', ['\itP\rm = ', num2str(thisPT, '%.3f'), ' (', chanRegion, ')']};
            if thisPT < 0.001
                dt.figTitle     = {'Place tuning', ['\itP\rm < 0.001 (', chanRegion, ')']};
            end
            dt.idx              = [iSub, chanIdx, iClus];
            dt.tPlace           = placeTest.tPlace;
            dt.tPlaceSurro      = tPlaceSurro;
            f                   = LK_PlotPlaceCell_20220330(dt);
            
            % save figure
            set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
            print(f, strcat(paths.save, subjects{iSub}, '_', chans(iChan).name, '_clus', num2str(iClus), '_loc_20220330'), '-dtiff', '-r450');
            
            % save figure data
            if all(dt.idx == [15, 14, 2], 2)
                save(strcat(paths.save, 'Fig_6a_upperLeft'), 'dt');
            elseif all(dt.idx == [16, 7, 2], 2)
                save(strcat(paths.save, 'Fig_6a_upperMiddleLeft'), 'dt');
            elseif all(dt.idx == [17, 19, 1], 2)
                save(strcat(paths.save, 'Fig_6a_upperMiddleRight'), 'dt');
            elseif all(dt.idx == [27, 7, 1], 2)
                save(strcat(paths.save, 'Fig_6a_upperRight'), 'dt');
            elseif all(dt.idx == [28, 43, 1], 2)
                save(strcat(paths.save, 'Fig_6a_lowerLeft'), 'dt');
            elseif all(dt.idx == [31, 9, 1], 2)
                save(strcat(paths.save, 'Fig_6a_lowerMiddleLeft'), 'dt');
            elseif all(dt.idx == [32, 9, 5], 2)
                save(strcat(paths.save, 'Fig_6a_lowerMiddleRight'), 'dt');
            elseif all(dt.idx == [39, 126, 1], 2)
                save(strcat(paths.save, 'Fig_6a_lowerRight'), 'dt');
            end

            %% close all open figures
            close all;
        end
    end
end

%% results
save(strcat(paths.save, 'results'), '-v7.3');

%% load previous results
r = load(strcat(paths.save, 'results.mat'));

%% general

% indices
allIdx          = cell2mat({r.allRes.idx}');
allSubjectName  = {r.allRes.subjectName}';
fprintf('Number of cells: %d.\n', size(allIdx, 1));

% brain regions
allRegions      = {r.allRes.chanRegion}';

%% place cells

% identify place cells based on the t-statistic as compared to the
% surrogate t-statistics
allTPlace       = cell2mat({r.allRes.tPlace}');
allTPlaceSurro  = cell2mat({r.allRes.tPlaceSurro})';
allTPlaceRank   = sum(allTPlace > allTPlaceSurro, 2) ./ sum(~isnan(allTPlaceSurro), 2);
bPlaceCell      = allTPlaceRank > 0.95;

% report
fprintf('Number of place cells based on a significant t-statistic: %d of %d (%.3f%%, P = %.3f).\n', sum(bPlaceCell), numel(bPlaceCell), ...
    100 * sum(bPlaceCell) / numel(bPlaceCell), myBinomTest(sum(bPlaceCell), numel(bPlaceCell), 0.05));

% save cell classification
cellClassification              = cellstr(repmat('NonPlaceCell', size(allIdx, 1), 1));
cellClassification(bPlaceCell)  = {'PlaceCell'};
cellClassification              = table(cellClassification, allSubjectName, allIdx(:, 1), allIdx(:, 2), allIdx(:, 3), ...
    'VariableNames', {'cellType', 'subjectName', 'session', 'channel', 'cluster'});
save(strcat(paths.save, 'cellClassification'), 'cellClassification');

%% place cells per brain region

% report
fprintf('\nAnalysis of place cells in different brain regions.\n');

% unique regions
uniqueRegions   = unique(allRegions);
numCellsPerReg  = nan(size(uniqueRegions, 1), 1);

% loop through the different regions
for iR = 1:size(uniqueRegions, 1)
    
    % cells from this region
    bThisRegion             = strcmp(allRegions, uniqueRegions{iR});
    numCellsPerReg(iR, 1)   = sum(bThisRegion);
    fprintf('Region: %s. Number of cells: %d.\n', uniqueRegions{iR}, sum(bThisRegion));
    
    % binomial test
    pBinom  = myBinomTest(sum(bPlaceCell & bThisRegion), sum(bThisRegion), 0.05);
    
    % report
    fprintf('\tNumber of place cells: %d (%.3f%%, P = %.3f).\n', sum(bPlaceCell & bThisRegion), ...
        100 * sum(bPlaceCell & bThisRegion) / sum(bThisRegion), pBinom);
end

% restrict unique regions to regions with enough units
uniqueRegions   = uniqueRegions(numCellsPerReg > 30);

% create figure for percentage of place cells per brain region
dt                      = [];
dt.uniqueRegions        = uniqueRegions;
dt.allRegions           = allRegions;
dt.bCell                = bPlaceCell;
dt.ylabel               = 'Place cells (%)';
dt.BonferroniCorrection = true;
f                       = LK_PlotCellsPerRegion(dt);

% save figure
LK_print(f, strcat(paths.save, 'placeCellsPerRegion_20230926'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6b'), 'dt');

%% properties: place fields; size of place fields; firing rate within vs. outside the place field

% preallocate and loop through cells
allPF           = nan(numel(r.place.yCenters), numel(r.place.xCenters), size(r.allRes, 1)); % place fields
allSmFRMap      = nan(numel(r.place.yCenters), numel(r.place.xCenters), size(r.allRes, 1)); % firing-rate maps
allPFSize       = nan(size(r.allRes, 1), 1); % place-field sizes
allFRInOut      = nan(size(r.allRes, 1), 2); % firing rates within vs. outside the place field
for iCell = 1:size(r.allRes, 1)
    
    % place fields
    allPF(:, :, iCell)      = r.allRes(iCell).PF; % orientation according to the indexing template
    
    % smoothed firing-rate map
    allSmFRMap(:, :, iCell) = r.allRes(iCell).smFRMap;
    
    % place-field size (as a fraction of the firing-rate map)
    allPFSize(iCell, 1)     = sum(sum(r.allRes(iCell).PF)) / sum(sum(~isnan(r.allRes(iCell).smFRMap)));
    
    % firing rate within vs. outside the place field
    allFRInOut(iCell, :)    = [mean(r.allRes(iCell).smFRMap(r.allRes(iCell).PF == 1), 'omitnan'), ...
        mean(r.allRes(iCell).smFRMap(r.allRes(iCell).PF == 0), 'omitnan')];
end

%% figure: overlay of all place fields

% sum of all place fields
sumPF       = 100 .* sum(allPF(:, :, bPlaceCell), 3) ./ sum(~isnan(allSmFRMap(:, :, bPlaceCell)), 3); % expressed as percentages

% create figure
dt          = [];
dt.place    = r.place;
dt.sumPF    = sumPF;
f = LK_PlotOverlayOfPlaceFields_20231012(dt);

% save figure
print(f, strcat(r.paths.save, 'placeCells_allPlaceFields_20220331'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6c'), 'dt');

%% figure: place-field size

% report on the size of the place fields
LK_ReportMeanAndSEM_20220322('Size of the place fields (%%)', 100 .* allPFSize(bPlaceCell)); % in percent

% create histogram
dt              = [];
dt.figPosition  = [2, 2, 5, 6];
dt.axPosition   = [1.4, 1.75, 3, 3.5];
dt.data         = 100 .* allPFSize(bPlaceCell); % in percent
dt.binEdges     = 100 .* (0:0.01:0.3);
dt.xline        = mean(100 .* allPFSize(bPlaceCell));
dt.xlabel       = 'Place-field size (%)';
dt.ylabel       = 'Count';
dt.fontSize     = 0.4;
f               = LK_HistogramFigure_20220321(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'placeCells_placeFieldSize_20220331'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6d'), 'dt');

%% figure: average firing rates within vs. outside the place fields

% report percentage increase in firing rates
FRDiff  = allFRInOut(:, 1) - allFRInOut(:, 2);
LK_ReportMeanAndSEM_20220322('FR difference within minus outside the place field', FRDiff(bPlaceCell));
FRRatio = allFRInOut(:, 1) ./ allFRInOut(:, 2);
LK_ReportMeanAndSEM_20220322('FR ratio within divided by outside the place field', FRRatio(bPlaceCell));

% report overall percentage increase in firing rates
fprintf('FR ratio within divided by outside the place field, after averaging across place cells: %.6f.\n', mean(allFRInOut(bPlaceCell, 1)) / mean(allFRInOut(bPlaceCell, 2)));

% stats
[~, pFR, ~, statsFR]    = ttest(allFRInOut(bPlaceCell, 1), allFRInOut(bPlaceCell, 2));
fprintf('Are the firing rates inside the place fields different from those outside the place fields? t(%d) = %.3f, p = %.3f.\n', statsFR.df, statsFR.tstat, pFR);

% create figure
dt          = [];
dt.groups   = {'Inside', rgb('darkred'); 'Outside', rgb('darkblue')};
dt.data     = {allFRInOut(bPlaceCell, 1), allFRInOut(bPlaceCell, 2)};
f = LK_PlotFiringRateInOutPlaceField_20231012(dt);

% save figure
LK_print(f, strcat(paths.save, 'placeCells_FRWithinOutsidePF_20230922'), '-dpng', '-r600');

% save figure data
save(strcat(paths.save, 'Fig_6g'), 'dt');

%% property: temporal stability of firing-rate maps

% all temporal stabilities
allTemporalStab     = cell2mat({r.allRes.temporalStab})';

% t-test of temporal stabilities against zero
[~, p, ~, stats]    = ttest(allTemporalStab(bPlaceCell));
fprintf('Temporal stability of place cells against zero: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);
LK_ReportMeanAndSEM_20220322('Temporal-stability values of place cells', allTemporalStab(bPlaceCell));
fprintf('Temporal-stability values of place cells range between %.3f and %.3f.\n', min(allTemporalStab(bPlaceCell)), max(allTemporalStab(bPlaceCell)));

% histogram: temporal stability
dt              = [];
dt.figPosition  = [2, 2, 5, 6];
dt.axPosition   = [1.4, 1.75, 3, 3.5];
dt.data         = allTemporalStab(bPlaceCell);
dt.binEdges     = -0.7:0.1:0.7;
dt.xline        = mean(allTemporalStab(bPlaceCell));
dt.xlabel       = 'Temporal stability (\itr\rm)';
dt.ylabel       = 'Count';
dt.fontSize     = 0.4;
f               = LK_HistogramFigure_20220321(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'placeCells_temporalStab_20220331'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6h'), 'dt');

%% overlap between place cells and object cells

% load object-cell results
objCC       = load(strcat(r.paths.cueObj, 'cellClassification.mat'));
bObjCell    = strcmp(objCC.cellClassification.cellType, 'ObjectCell');
fprintf('\nExamining the overlap between object cells and place cells.\n');
fprintf('Number of object cells: %d (out of %d cells).\n', sum(bObjCell), size(bObjCell, 1));

% double-check whether cellular indices are identical
objIdx      = [objCC.cellClassification.session, objCC.cellClassification.channel, objCC.cellClassification.cluster];
if ~all(all(allIdx == objIdx, 2))
    error('Cellular indices from the object-cell analysis and the place-cell analysis are not identical.');
end

% overlap between place cells and object cells
fprintf('Number of cells that are both place and object cells: %d.\n', sum(bPlaceCell & bObjCell));

% figure: Venn diagram for the overlap between place and object cells
dt      = [];
dt.A    = [sum(bObjCell), sum(bPlaceCell)];
dt.I    = sum(bObjCell & bPlaceCell);
f = LK_PlotVennDiagramTwoSamples_20231012(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'Venn_PlaceCells_vs_ObjCells_20220331'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_5f'), 'dt');

%% test whether place fields touch the borders of the arena

% preallocate
nBinsPF         = nan(size(r.allRes, 1), 1); % number of bins of the place field
nBinsPFAtEdge   = nan(size(r.allRes, 1), 1); % number of bins of the place field at an edge of the firing-rate map
fracPFAtEdge    = nan(size(r.allRes, 1), 1); % fraction of the place field that is at an edge of the firing-rate map

% loop through cells
for iCell = 1:size(r.allRes, 1)
    
    % skip if not a place field
    if bPlaceCell(iCell) == false
        continue;
    end
    
    % firing-rate map and place field of this cell
    thisFRMap   = r.allRes(iCell).smFRMap;
    thisPF      = r.allRes(iCell).PF;
    
    % fill holes
    thisFRMapFilled = imfill(~isnan(thisFRMap), 'holes'); % fills holes in the input image
    thisPFFilled    = imfill(thisPF, 'holes');
    
    % perimeters of the filled place field and the filled firing-rate map
    thisFRMapPerim  = bwperim(thisFRMapFilled); % returns the perimeter
    thisPFPerim     = bwperim(thisPFFilled);
    
    % bins that are edges of both the place field and the firing-rate map
    bBothEdges      = (thisFRMapPerim == 1) & (thisPFPerim == 1);
    
    % percentage of the place field that is at an edge
    nBinsPF(iCell, 1)       = sum(thisPF(:));
    nBinsPFAtEdge(iCell, 1) = sum(bBothEdges(:));
    fracPFAtEdge(iCell, 1)  = nBinsPFAtEdge(iCell, 1) / nBinsPF(iCell, 1);
end

% count how often the place field has some contact to the edges of the
% firing-rate map
bPFAtFREdge = fracPFAtEdge > 0;
fprintf('Number of place fields with at least some contact to the edge of the firing-rate map: %d (of %d).\n', sum(bPFAtFREdge & bPlaceCell), sum(bPlaceCell, 1));
LK_ReportMeanAndSEM_20220322('Fraction of the place fields being at an edge (%%)', 100 .* fracPFAtEdge(bPlaceCell));

%% figure: histogram of place-field fraction touching the firing-rate map edges

% create histogram
dt              = [];
dt.figPosition  = [2, 2, 5, 6];
dt.axPosition   = [1.4, 1.75, 3, 3.5];
dt.data         = 100 .* fracPFAtEdge(bPlaceCell); % in percent
dt.binEdges     = 100 .* (0:0.05:1);
dt.xline        = mean(100 .* fracPFAtEdge(bPlaceCell));
dt.xlabel       = '% of place field at edge    ';
dt.ylabel       = 'Count';
dt.fontSize     = 0.4;
f               = LK_HistogramFigure_20220321(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'placeCells_PFFractionAtEdge_20230403'), '-dpng', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6e'), 'dt');

%% test whether place fields contain more object locations than expected by chance

% reset rng for reproducibility
rng(r.param.myRNG);

% report
fprintf('\nInvestigation of objects within place fields.\n');

% preallocate
numObjInPF      = nan(size(r.allRes, 1), 1); % number of objects within the place field, separately for each cell
numObjInPFSurro = nan(size(r.allRes, 1), r.param.numSurrogates); % surrogate number of objects within the place field

% loop through cells
for iCell = 1:size(r.allRes, 1)
    
    % skip if this is not a place cell
    if bPlaceCell(iCell) == false
        continue;
    end
    fprintf('\tCell index: %d (of %d).\n', iCell, size(r.allRes, 1));
    
    % trial information for object locations
    trialInfo   = load(strcat(r.paths.beh, r.subjects{r.allRes(iCell).idx(1), 1}, '\trialInfo.mat'), 'trialInfo');
    objLocs    	= [trialInfo.trialInfo.xCorrect, trialInfo.trialInfo.yCorrect];
    objLocs   	= unique(objLocs(~isnan(objLocs(:, 1)), :), 'rows', 'stable');
    
    % bin the object locations (like the firing rates)
    xBinsObj    = discretize(objLocs(:, 1), r.place.xEdges);
    yBinsObj    = discretize(objLocs(:, 2), r.place.yEdges);
    xyBinsObj   = (xBinsObj - 1) .* numel(r.place.yCenters) + yBinsObj;
    
    % count how many discretized object locations are part of the place
    % field
    xyBinsPF                = r.place.idxTemplate(r.allRes(iCell).PF);
    numObjInPF(iCell, 1)    = sum(ismember(xyBinsObj, xyBinsPF));
    
    % create surrogates
    for iSurro = 1:r.param.numSurrogates
        
        % draw surrogate place field
        surroPF = LK_DrawConnPixelsFromArea_20230304(~isnan(r.allRes(iCell).smFRMap), sum(r.allRes(iCell).PF(:))); % mask, number of connected pixels to draw

        % count the number of objects in the surrogate place field
        if isempty(surroPF)
            continue;
        else
            xyBinsPFSurro                   = r.place.idxTemplate(logical(surroPF));
            numObjInPFSurro(iCell, iSurro)  = sum(ismember(xyBinsObj, xyBinsPFSurro));
        end
    end
end

%% statistics and figure: number of objects within place fields as compared to surrogates

% report
LK_ReportMeanAndSEM_20220322('Number of objects per place field', numObjInPF(bPlaceCell, :));

% stats: two-sample Kolmogorov-Smirnov test
data1                   = numObjInPF(bPlaceCell, :);
data2                   = numObjInPFSurro(bPlaceCell, :);
[~, pObj, ksstatObj]    = kstest2(data1, data2(:));
fprintf('Is the number of objects per place field increased as compared to surrogates? D = %.3f, p = %.3f.\n', ksstatObj, pObj);

% create figure
dt          = [];
dt.data1    = data1;
dt.data2    = data2;
dt.pObj     = pObj;
f = LK_PlotNumObjInsidePF_20231013(dt);

% save figure
LK_print(f, strcat(r.paths.save, 'placeCells_numObjInPF_20230403'), '-dpng', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_6f'), 'dt');
