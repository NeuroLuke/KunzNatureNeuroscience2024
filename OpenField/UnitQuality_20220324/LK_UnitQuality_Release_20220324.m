%==========================================================================
% This script performs a unit quality assessment of the recorded units.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\WaveClus3\')); % Chaure et al., 2018

% variables
param               = [];
param.wcFileName    = 'times_datacut.mat';
param.c4aName       = 'LK_Cluster4Analysis_20210521.mat';
param.ccmName       = 'cluster_class_macrotime.mat';
param.peakIdx       = 20; % peak index of the waveform

% paths
paths         	= [];
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.micro     = 'E:\OpenField\Micro_20201215\';
paths.spike   	= 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.place     = 'E:\OpenField\PlaceCellAnalysis_20210511\20211116_Loc25x25_facThreshPF0.75\';
paths.save    	= 'E:\OpenField\UnitQuality_20220324\Evaluation_20230330\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% subjects
subjects        = load(strcat(paths.subjects, 'subjects.mat'));
subjects        = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% save settings
save(strcat(paths.save, 'settings'));

%% load results from place-cell analysis

% in order to only examine cells that were also included in the analyses
placeRes        = load(strcat(paths.place, 'results.mat'), 'allRes');

%% preallocations

% results for all units
allRes          = [];

% bookkeeping
exByWC          = []; % excluded by wave-clus
exByNoise       = []; % excluded because of being the noise cluster
exByVisInsp     = []; % excluded based on visual inspection
exByPlace       = []; % excluded from the place-cell analysis

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report and reset rng for reproducibility
    fprintf('\n\nSUBJECT: %s.\n', subjects{iSub});
    rng(999);
    if ~strcmp(subjects{iSub, 2}, 'Microwire')
        continue;
    end
    
    %% wires to investigate
    
    % available microwires
    wires   = LK_GetChans(fullfile(paths.spike, subjects{iSub}, 'chan*'));
    
    %% loop through wires
    for iWire = 1:size(wires, 1)
        
        % report
        fprintf('\tWire: %s.\n', wires(iWire).name);
        
        % original wire index
        wireIdx     = str2double(replace(wires(iWire).name, lettersPattern, ''));
        
        %% get spike times
        
        % load wave-clus output
        try
            times_datacut   = load(fullfile(wires(iWire).folder, wires(iWire).name, param.wcFileName));
        catch
            fprintf('\t- No wave-clus for this wire.\n');
            exByWC          = cat(1, exByWC, [iSub, wireIdx]); % bookkeeping
            continue;
        end
        
        % load decision whether to use the clusters (based on visual
        % inspection)
        c4a = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, param.c4aName));
        
        % load spike times in macrotime for sanity checking
        ccm = load(fullfile(wires(iWire).folder, wires(iWire).name, param.ccmName));
        
        % sanity check
        if size(times_datacut.cluster_class, 1) ~= size(ccm.cluster_class_macrotime, 1)
            error('Size of original "cluster_class" not congruent with new "cluster_class_macrotime".');
        end
        
        %% create bandpass-filtered raw signal
        
        % load raw data
        raw                 = load(strcat(paths.micro, subjects{iSub}, filesep, wires(iWire).name, '\datacut.mat'));
        
        % filter between 300 and 3,000 Hz using wave-clus
        par                 = [];
        par.sr              = raw.sr;
        par.detect_fmin     = times_datacut.par.detect_fmin;
        par.detect_fmax     = times_datacut.par.detect_fmax;
        par.detect_order    = times_datacut.par.detect_order;
        xf_detect           = spike_detection_filter(raw.data, par);
        
        % calculate STD of the filtered data
        STDnoise            = std(xf_detect);
        
        %% loop through clusters of this wire
        for iClus = 0:max(times_datacut.cluster_class(:, 1))
            
            % report
            fprintf('\t\tCluster: %d.\n', iClus);
            
            % continue if this cluster is the noise cluster
            if iClus == 0
                exByNoise   = cat(1, exByNoise, [iSub, wireIdx, iClus]); % bookkeeping
                continue;
            end
            
            % continue if you decided that this cluster is not sufficient
            if strcmp(c4a.Cluster4Analysis(cell2mat({c4a.Cluster4Analysis.cluster}) == iClus).decision, 'no')
                fprintf('\t\t- You decided not to analyse this cluster.\n');
                exByVisInsp = cat(1, exByVisInsp, [iSub, wireIdx, iClus]); % bookkeeping
                continue;
            end
            
            % continue if this cluster was not included in the place-cell
            % analysis
            if ~any(all([iSub, wireIdx, iClus] == cell2mat({placeRes.allRes.idx}'), 2))
                fprintf('\t\t- Was not included in the analysis, thus skipping.\n');
                exByPlace     = cat(1, exByPlace, [iSub, wireIdx, iClus]); % bookkeeping
                continue;
            end
            
            % data for this cluster
            thisCluster     = times_datacut.cluster_class(times_datacut.cluster_class(:, 1) == iClus, :); % cluster-number, microtime (msec)
            thisSpike       = times_datacut.spikes(times_datacut.cluster_class(:, 1) == iClus, :); % spike-idx X time-within-spike
            
            %% ISI refractoryness
            % = percentage of ISIs < 3ms
            
            ISI                 = diff(thisCluster(:, 2)); % inter-spike-intervals (ms)
            nISI                = size(ISI, 1);
            percISIlessThan3ms  = 100 * sum(ISI < 3) / nISI;
            
            %% mean firing rate
            
            % mean firing rate relative to the entire recording
            expDur      = (size(raw.data, 2) - 1) / raw.sr;
            meanFR  	= size(thisCluster, 1) / expDur; % (Hz)
            
            %% waveform peak SNR
            % = ratio between the peak amplitude of the mean waveform and
            % the STD of the noise
            
            peakAmpl   	= abs(mean(thisSpike(:, param.peakIdx)));
            peakSNR  	= peakAmpl / STDnoise;
                        
            %% collect information across units
            
            % this unit's results
            unitRes                     = [];
            unitRes.idx                 = [iSub, wireIdx, iClus];
            unitRes.percISIlessThan3ms  = percISIlessThan3ms;
            unitRes.meanFR              = meanFR;
            unitRes.peakSNR             = peakSNR;
            
            % collapse across units
            allRes  = cat(1, allRes, unitRes);
        end
    end
end

%% save important output
save(strcat(paths.save, 'results'), '-v7.3');

%% results

% load previously saved results
r   = load(strcat(paths.save, 'results.mat'));

% unit indices
allUnitIdx  = cell2mat({r.allRes.idx}');

% report
fprintf('\nResults.\n');
fprintf('Number of units: %d.\n', size(allUnitIdx, 1));

%% units per wire

% number of units per wire
uniqueAllWireIdx    = unique(allUnitIdx(:, 1:2), 'rows', 'stable');
numUnitsPerWire     = nan(size(uniqueAllWireIdx, 1), 1);
for iWire = 1:size(uniqueAllWireIdx, 1)
    numUnitsPerWire(iWire)  = sum(all(uniqueAllWireIdx(iWire, 1:2) == allUnitIdx(:, 1:2), 2));
end

% average number of units per wire
fprintf('Number of wires with at least one unit that was analyzed: %d.\n', size(uniqueAllWireIdx, 1));
LK_ReportMeanAndSEM_20220322('Number of units per wire', numUnitsPerWire);

% histogram showing the number of units per wire
cfg             = [];
cfg.figPosition = [2, 2, 6, 6];
cfg.axPosition  = [1.75, 1.75, 3.5, 3.5];
cfg.data        = numUnitsPerWire;
cfg.xlabel      = 'Units per wire';
cfg.ylabel      = 'Number of wires';
cfg.fontSize    = 0.4;
f               = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'UnitsPerWire_20220324'), '-dtiff', '-r300');

%% ISI refractoriness

% ISI refractoriness from all cells
allPercISILessThan3ms   = cell2mat({r.allRes.percISIlessThan3ms}');

% report
LK_ReportMeanAndSEM_20220322('Average percentage of ISIs <3 ms across units', allPercISILessThan3ms);
fprintf('Number of units with a percentage of ISIs <3 ms higher than 5%%: %d.\n', ...
    sum(allPercISILessThan3ms >= 5));

% histogram showing the ISI refractoriness per unit
cfg             = [];
cfg.figPosition = [2, 2, 6, 6];
cfg.axPosition  = [1.75, 1.75, 3.5, 3.5];
cfg.data        = allPercISILessThan3ms;
cfg.binEdges    = 0:0.5:6;
cfg.xlabel      = '% ISI <3 ms';
cfg.ylabel      = 'Number of units';
cfg.fontSize    = 0.4;
f               = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'ISIRefractoriness_20220324'), '-dtiff', '-r300');

%% firing rate

% mean firing rate from all cells
allMeanFR       = cell2mat({r.allRes.meanFR}');

% report
LK_ReportMeanAndSEM_20220322('Mean firing rate across all units (Hz)', allMeanFR);

% histogram showing the mean firing rate per unit
cfg             = [];
cfg.figPosition = [2, 2, 6, 6];
cfg.axPosition  = [1.75, 1.75, 3.5, 3.5];
cfg.data        = allMeanFR;
cfg.binEdges    = 0:1:25;
cfg.xlabel      = 'Mean FR (Hz)';
cfg.ylabel      = 'Number of units';
cfg.fontSize    = 0.4;
f               = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'MeanFR_20220324'), '-dtiff', '-r300');

%% peak SNR

% peak SNR from all cells
allPeakSNR      = cell2mat({r.allRes.peakSNR}');

% report
LK_ReportMeanAndSEM_20220322('Average peak SNR', allPeakSNR);

% histogram showing the peak SNR per unit
cfg             = [];
cfg.figPosition = [2, 2, 6, 6];
cfg.axPosition  = [1.75, 1.75, 3.5, 3.5];
cfg.data        = allPeakSNR;
cfg.binEdges    = 0:3:33;
cfg.xlabel      = 'Peak SNR';
cfg.ylabel      = 'Number of units';
cfg.fontSize    = 0.4;
f               = LK_HistogramFigure_20220321(cfg);
LK_print(f, strcat(paths.save, 'PeakSNR_20220324'), '-dtiff', '-r300');
