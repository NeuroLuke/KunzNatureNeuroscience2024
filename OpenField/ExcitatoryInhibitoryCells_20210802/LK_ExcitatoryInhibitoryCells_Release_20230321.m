%==========================================================================
% This script distinguishes between putatively excitatory and inhibitory
% cells and investigates their overlap with object/place cells.
%
% References:
% - Viskontas et al., Hippocampus, 2007
% - Ison et al., Journal of Neurophysiology, 2011
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param                   = [];
param.clusterName       = 'LK_Cluster4Analysis_20210521';
param.subjectFileName   = 'subjectdata_20210616.mat';

% paths
paths               = [];
paths.subjects      = 'E:\OpenField\SubjectData_20210616\';
paths.spikes        = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.save          = 'E:\OpenField\ExcitatoryInhibitoryCells_20210802\20230320\';
paths.objCells      = 'E:\OpenField\CueObjAnalysis_20210812\20211116\'; % object cells
paths.placeCells    = 'E:\OpenField\PlaceCellAnalysis_20210511\20211116_Loc25x25_facThreshPF0.75\'; % place cells
mkdir(paths.save);

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% preallocate

% results for all units
allRes          = [];

% for controls
exByWC          = [];
exByVisInsp     = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)

    % report
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    
    % skip if it is not a microwire subject
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        fprintf('\nSkipping subject %s because no microwires.\n', subjects{iSub});
        continue;
    end
    
    %% channels with microwire data    
    
    % available channels
    chans   = LK_GetChans(fullfile(paths.spikes, subjects{iSub}, 'chan*'));
    
    % report
    fprintf('\nProcessing of neural data (%d channels from "%s" to "%s").\n', size(chans, 1), chans(1).name, chans(end).name);
    
    % total duration of the recording (to calculate mean firing rate)
    ainp2File   = load(fullfile(paths.spikes, subjects{iSub}, 'ainp2\datacut.mat'));
    totalDur    = numel(ainp2File.data) / ainp2File.sr; % in seconds
    
    %% loop through wires
    for iChan = 1:size(chans, 1)
        
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
            exByWC  = cat(1, exByWC, [iSub, chanIdx]); % bookkeeping
            continue;
        end
        
        % load decision whether to use clusters (based on inspection)
        c4aFile = strcat(chans(iChan).folder, filesep, chans(iChan).name, filesep, param.clusterName, '.mat');
        c4a     = load(c4aFile);
        
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
            
            %% firing rate (Viskontas et al., 2007)
            % procedure: total number of spikes divided by the total
            % duration of the recording session
            
            % mean firing rate
            meanFR  = size(thisCluster, 1) / totalDur;
            
            %% burstiness = interspike interval (ISI) ratio (Viskontas et al., 2007)
            % procedure: divide the total number of spikes that occur at
            % less than 10 ms of each other by the total number of spikes
            % that occur at greater than or equal to 10 ms of each other
            
            % inter-spike intervals
            ISI         = diff(thisCluster(:, 2)); % (ms)
            
            % estimate inter-spike interval ratio to quantify burstiness
            burstiness  = sum(ISI < 10) / sum(ISI >= 10);
            
            %% spike width (Ison et al., 2011)
            
            % original wave form
            waveForm    = mean(thisSpike, 1);
            waveTime    = (0:(numel(waveForm) - 1)) / times_datacut.par.sr * 1000; % in milliseconds
            
            % global minimum
            [~, minIdx] = min(waveForm); % global minimum
            tMin        = waveTime(minIdx);
            
            % maximum after global minimum
            tmpWaveForm             = waveForm;
            tmpWaveForm(1:minIdx)  	= nan;
            [~, maxIdx]         	= max(tmpWaveForm);
            tMax                  	= waveTime(maxIdx);
            
            % spike width
            spikeWidth  = tMax - tMin;
            
            %% results from this unit
            
            % basics
            unitRes           	= [];
            unitRes.idx        	= [iSub, chanIdx, iClus];
            unitRes.chanRegion 	= chanRegion;
            
            % average waveform
            unitRes.waveForm    = waveForm;
            unitRes.waveTime  	= waveTime;
            unitRes.tMin        = tMin;
            unitRes.tMax        = tMax;
            
            % metrics for classification
            unitRes.meanFR   	= meanFR; % Hz
            unitRes.burstiness 	= burstiness; % ratio
            unitRes.spikeWidth	= spikeWidth; % duration (ms)
            
            % concatenate across units
            allRes  = cat(1, allRes, unitRes);
            
            %% close figures
            close all;
        end
    end
end

%% save results
save(strcat(paths.save, 'results'), '-v7.3');

%% load results
r = load(strcat(paths.save, 'results.mat'));

% all indices
allIdx          = cell2mat({r.allRes.idx}');
allChanRegions  = {r.allRes.chanRegion}';

% report
fprintf('Number of cells: %d.\n', size(allIdx, 1));

%% all characteristics

% average firing rates, burstiness, and spike width
allMeanFR    	= cell2mat({r.allRes.meanFR}'); % non-normal distribution
allBurstiness   = cell2mat({r.allRes.burstiness}'); % non-normal distribution
allSpikeWidth   = cell2mat({r.allRes.spikeWidth}');

%% kmeans-clustering to differentiate between putative pyramidal cells and putative interneurons

% report
fprintf('\nSeparating the cells based on k-means-clustering.\n');

% kmeans clustering using firing rate, burstiness, and spike width
groupIdx        = kmeans([allMeanFR, allBurstiness, allSpikeWidth], 2);

% ensure that group #1 codes for the pyramidal neurons
if mean(allMeanFR(groupIdx == 1)) > mean(allMeanFR(groupIdx == 2))
    groupIdx    = (groupIdx == 1) + 1; % 1 = pyramidal cells; 2 = interneurons
end
groupColors     = {rgb('blue'), rgb('red')};

% report number of cells
fprintf('\nNumber of putative pyramidal cells: %d. Number of putative interneurons: %d.\n', sum(groupIdx == 1), sum(groupIdx == 2));

% report differences between the groups regarding firing rate
[~, p, ~, stats] = ttest2(allMeanFR(groupIdx == 1), allMeanFR(groupIdx == 2));
fprintf('Mean firing rate (Hz): %.3f for group #1; %.3f for group #2; t(%d) = %.3f, p = %.3f.\n', mean(allMeanFR(groupIdx == 1)), mean(allMeanFR(groupIdx == 2)), ...
    stats.df, stats.tstat, p);

% report differences between the groups regarding burstiness
[~, p, ~, stats] = ttest2(allBurstiness(groupIdx == 1), allBurstiness(groupIdx == 2));
fprintf('Burstiness (ISI ratio): %.3f for group #1; %.3f for group #2; t(%d) = %.3f, p = %.3f.\n', mean(allBurstiness(groupIdx == 1)), mean(allBurstiness(groupIdx == 2)), ...
    stats.df, stats.tstat, p);

% report differences between spike widths
[~, p, ~, stats] = ttest2(allSpikeWidth(groupIdx == 1), allSpikeWidth(groupIdx == 2));
fprintf('Spike width (ms): %.3f for group #1; %.3f for group #2; t(%d) = %.3f, p = %.3f.\n', mean(allSpikeWidth(groupIdx == 1)), mean(allSpikeWidth(groupIdx == 2)), ...
    stats.df, stats.tstat, p);

%% figure for separating units into pyramidal cells and interneurons based on different features

% create figure: 3D plot
f = figure('units', 'centimeters', 'position', [2, 2, 8.5, 8.5]);
% cell separation
axes('units', 'centimeters', 'position', [1.3, 2.4, 6, 6]);
hold on;
for iGroup = min(groupIdx):max(groupIdx)
    plot3(allSpikeWidth(groupIdx == iGroup), allBurstiness(groupIdx == iGroup), allMeanFR(groupIdx == iGroup), '.', 'Color', groupColors{iGroup}, 'MarkerSize', 5);
end
xl = xlabel('Spike width (ms)', 'units', 'normalized', 'position', [0.55, -0.05]);
yl = ylabel('Burstiness', 'units', 'normalized', 'position', [0.3, 0]);
zl = zlabel('Firing rate (Hz)');
view([135, 25]);
t1 = text(0.6, 0.5, {'Putative', 'pyramidal cells'}, 'FontUnits', 'centimeters', 'FontSize', 0.4, 'Color', groupColors{1}, 'units', 'normalized');
t2 = text(0.6, 0.7, {'Putative', 'interneurons'}, 'FontUnits', 'centimeters', 'FontSize', 0.4, 'Color', groupColors{2}, 'units', 'normalized');
grid minor;
set([gca, xl, yl, zl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
% explanation of spike width
axes('units', 'centimeters', 'position', [5, 0.2, 2, 1.3]);
hold on;
plot(r.allRes(1).waveTime, r.allRes(1).waveForm, 'Color', [0, 0, 0]);
tmpY = min(r.allRes(1).waveForm);
plot([r.allRes(1).tMin, r.allRes(1).tMin, r.allRes(1).tMin, r.allRes(1).tMax, r.allRes(1).tMax, r.allRes(1).tMax], ...
    [tmpY * 1.1, tmpY * 1.3, tmpY * 1.2, tmpY * 1.2, tmpY * 1.3, tmpY * 1.1], '-', 'Color', [0.5, 0.5, 0.5]);
axis off;
% save figure
LK_print(f, strcat(paths.save, 'CellFeatures_3D_20230325'), '-dpng', '-r600');

%% investigate whether place cells and/or object cells are particularly prevalent among pyramidal cells or interneurons

% report
fprintf('\nLoading cell classifications.\n');

% object cells and object results
objCells    = load(strcat(paths.objCells, 'cellClassification.mat'));
bObjCell    = strcmp(objCells.cellClassification.cellType, 'ObjectCell');
fprintf('Number of object cells: %d (out of %d).\n', sum(bObjCell), size(bObjCell, 1));

% place cells
placeCells  = load(strcat(paths.placeCells, 'cellClassification.mat'));
bPlaceCell  = strcmp(placeCells.cellClassification.cellType, 'PlaceCell');
fprintf('Number of place cells: %d (out of %d).\n', sum(bPlaceCell), size(bPlaceCell, 1));

% double-check cell indices
if ~all(all(allIdx == objCells.cellClassification{:, 3:5}, 2))
    error('Problem with cell indices between this analysis and the object-cell analysis.');
end
if ~all(all(allIdx == placeCells.cellClassification{:, 3:5}, 2))
    error('Problem with cell indices between this analysis and the place-cell analysis.');
end

% cells that are an object cell and/or a place cell
bOPCell = bObjCell | bPlaceCell;
fprintf('Number of cells that are either an object cell or a place cell: %d.\n', sum(bOPCell));

%% overlap between object/place cells and pyramidal cells vs. interneurons

% Fisher exact test
X = table([sum(bOPCell & groupIdx == 1); sum(~bOPCell & groupIdx == 1)], [sum(bOPCell & groupIdx == 2); sum(~bOPCell & groupIdx == 2)], ...
    'VariableNames', {'Pyr', 'Int'}, 'RowNames', {'OP', 'nonOP'});
[~, pFisher, statsFisher] = fishertest(X);
fprintf('Are object/place cells more prevalent among pyramidal cells or interneurons? Fisher''s p = %.3f.\n', pFisher);

% figure: percentage of responsive cells among pyramidal cells and
% interneurons
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.5, 1.4, 4, 4]);
hold on;
bar(1, 1, 'FaceColor', [1, 1, 1]);
bar(1, sum(bOPCell & groupIdx == 1) / sum(groupIdx == 1), 'FaceColor', groupColors{1});
bar(2, 1, 'FaceColor', [1, 1, 1]);
bar(2, sum(bOPCell & groupIdx == 2) / sum(groupIdx == 2), 'FaceColor', groupColors{2});
xl = xlabel('Putative cell type');
yl = ylabel('Fraction object/place cells');
tl = title(sprintf('\\itP\\rm = %.3f', pFisher));
set(gca, 'xlim', [0.4, 2.6], 'xtick', [1, 2], 'xticklabel', {'Pyr', 'Int'}, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
LK_print(f, strcat(paths.save, 'OverlapObjectPlaceCellsWithPyrInt_20230326'), '-dpng', '-r300');

%% compare different characteristics between responsive cells and non-responsive cells

% conditions
conditions  = {'meanFR', 'Mean FR (log Hz)'; 'burstiness', 'Burstiness (log)'; 'spikeWidth', 'Spike width (ms)'};
for iCond = 1:size(conditions, 1)

    % select data
    if strcmp(conditions{iCond, 1}, 'meanFR')
        data    = allMeanFR;
    elseif strcmp(conditions{iCond, 1}, 'burstiness')
        data    = allBurstiness;
    elseif strcmp(conditions{iCond, 1}, 'spikeWidth')
        data    = allSpikeWidth;
    end
    
    % separate data for different groups of cells
    d   = {data(bOPCell), data(~bOPCell)};
    m   = {mean(d{1}, 'omitnan'), mean(d{2}, 'omitnan')};
    ste = {LK_ste(d{1}), LK_ste(d{2})};
    
    % non-parametric statistics
    [p, ~, stats] = ranksum(d{1}, d{2});
    fprintf('\nIs "%s" higher in object/place cells than in other cells? z = %.3f, p = %.3f.\n', conditions{iCond, 1}, stats.zval, p);
    LK_ReportMeanAndSEM_20220322(['Object/place cells, ', conditions{iCond, 2}], d{1});
    LK_ReportMeanAndSEM_20220322(['Other cells, ', conditions{iCond, 2}], d{2});

    %% figure
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
    axes('units', 'centimeters', 'position', [1.5, 1.4, 4, 4]);
    hold on;
    for iM = 1:numel(m)
        plot(iM + rand(numel(d{iM}), 1) * 0.5 - 0.25, d{iM}, '.', 'Color', [0.8, 0.8, 0.8]);
        bar(iM, m{iM}, 'FaceColor', 'none');
        plot([iM, iM], [m{iM} - ste{iM}, m{iM} + ste{iM}], '-', 'Color', [0, 0, 0], 'LineWidth', 4);
    end
    xl = xlabel('cells');
    yl = ylabel(conditions{iCond, 2}, 'units', 'normalized', 'position', [-0.1, 0.5]);
    tl = title(sprintf('\\itP\\rm = %.3f', p));
    if p < 0.001
        tl = title('\itP\rm < 0.001');
    end
    if contains(conditions{iCond, 1}, {'meanFR', 'burstiness'})
        set(gca, 'yscale', 'log');
    end
    tmpAx = get(gca);
    set(gca, 'xlim', [0.4, 2.6], 'xtick', [1, 2], 'xticklabel', {'Object/Place', 'Other'}, 'xticklabelrotation', 0, ...
        'ylim', [min(tmpAx.YLim), max(tmpAx.YLim)], 'ytick', [min(tmpAx.YTick), max(tmpAx.YTick)], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    LK_print(f, strcat(paths.save, conditions{iCond, 1}, '_ObjectPlaceVSOtherCells_20230326'), '-dpng', '-r300');
end
