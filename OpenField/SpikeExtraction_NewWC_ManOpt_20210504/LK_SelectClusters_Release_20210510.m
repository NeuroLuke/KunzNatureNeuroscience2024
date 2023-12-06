%==========================================================================
% Double-check which units to use for the analyses.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;

% paths
paths           = [];
paths.data      = 'C:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';

% parameters
param           = [];
param.bRecheck  = false;
param.saveName  = 'LK_Cluster4Analysis_20210521';

% subjects
subjects    = {...
    'FR006a'; 'FR006b'; ...
    'FR007'; ...
    'FR008a'; 'FR008b'; ...
    'FR009a'; 'FR009b'; ...
    'FR010'; ...
    'FR011'; ...
    'FR013'; ...
    'FR014'; ...
    'FR017'; ...
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

% bookkeeping
allUnits    = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)

    % report
    fprintf('\n==== SUBJECT: %s.\n', subjects{iSub});
    
    % channels
    chans   = dir(strcat(paths.data, subjects{iSub}, '\chan*'));
    splits  = split({chans.name}', 'n');
    [~, I]  = sort(cellfun(@str2double, splits(:, 2)));
    chans   = chans(I);
        
    %% loop through channels
    for iChan = 1:size(chans, 1)
        
        % preallocate output
        Cluster4Analysis    = [];
        
        % output file
        outputFile  = strcat(chans(iChan).folder, '\', chans(iChan).name, '\', param.saveName, '.mat');
        if exist(outputFile, 'file')
            % load previous output if you already performed the analysis
            c4a     = load(outputFile);
        else
            c4a     = [];
        end
        
        % load wave-clus output for this channel
        wcFile  = strcat(chans(iChan).folder, '\', chans(iChan).name, '\times_datacut.mat');
        if exist(wcFile, 'file')
            td  = load(wcFile);
        else
            % save empty structure and continue if there is no waveclus
            % output
            save(outputFile, 'Cluster4Analysis');
            continue;
        end
        
        % different units recorded on this channel
        clusters    = unique(td.cluster_class(:, 1));
        
        %% loop through different clusters
        for iClus = min(clusters):max(clusters)
            
            % cluster index
            thisCluster4Analysis                = [];
            thisCluster4Analysis.cluster        = iClus;
            
            % exclude the rest cluster
            if iClus == 0
                thisCluster4Analysis.decision   = 'no';
                Cluster4Analysis                = cat(1, Cluster4Analysis, thisCluster4Analysis);
                continue;
            end
            
            %% create figure with information about this unit
            
            % create figure
            f = figure('units', 'centimeters', 'position', [1, 10, 49, 12]);
            myColors = colormap('lines');
            if param.bRecheck == false
                set(f, 'visible', 'off');
            end
            
            % waveforms
            ax1 = axes('units', 'centimeters', 'position', [2, 2, 10, 8]);
            hold on;
            thisSpikes  = td.spikes(td.cluster_class(:, 1) == iClus, :);
            thisTime    = 1000 .* ((1:size(thisSpikes, 2)) - 1) ./ td.par.sr;
            idx2plot    = unique(round(linspace(1, size(thisSpikes, 1), 1000))); % plot subselection
            numTotal    = sum(td.cluster_class(:, 1) == iClus);
            numForced   = sum(td.forced(td.cluster_class(:, 1) == iClus) == false);
            plot(thisTime, thisSpikes(idx2plot, :), 'Color', myColors(iClus, :), 'LineWidth', 0.1);
            plot(thisTime, mean(thisSpikes, 1), 'Color', [0, 0, 0], 'LineWidth', 2);
            set(gca, 'xlim', [min(thisTime), max(thisTime)]);
            xlabel('Time (ms)');
            ylabel('Voltage (\muV)');
            title(sprintf('Cluster %d: # %d (%d)', iClus, numTotal, numForced));
            
            % inter-spike-intervals
            axes('units', 'centimeters', 'position', [14, 2, 10, 8]);
            thisTimes   = td.cluster_class(td.cluster_class(:, 1) == iClus, 2);
            thisISI     = diff(thisTimes);
            numLess3    = sum(thisISI < 3);
            histogram(thisISI, 0:1:100, 'EdgeColor', [0, 0, 0], 'FaceColor', [0, 0, 0]);
            set(gca, 'xlim', [0, 100]);
            xlabel('ISI (ms)');
            ylabel('Number');
            title(sprintf('%d in < 3ms (%.1f%%)', numLess3, 100 * numLess3 / numTotal));
            
            % amplitude over time
            axes('units', 'centimeters', 'position', [26, 2, 10, 8]);
            thisAmp     = min(thisSpikes, [], 2); % global minimum
            plot(thisTimes ./ 60000, thisAmp, 'k.');
            set(gca, 'xlim', [min(thisTimes), max(thisTimes)] ./ 60000, 'ylim', ax1.YLim);
            xlabel('Time (min)');
            ylabel('Voltage (\muV)');
            title('Peak voltages over time');
            
            % all waveforms from this channel
            axes('units', 'centimeters', 'position', [38, 2, 10, 8]);
            hold on;
            set(gca, 'xlim', [min(thisTime), max(thisTime)]);
            for i = 1:max(td.cluster_class(:, 1))
                tmpSpikes   = td.spikes(td.cluster_class(:, 1) == i, :);
                idx2plot    = unique(round(linspace(1, size(tmpSpikes, 1), 200))); % plot subselection
                plot(thisTime, tmpSpikes(idx2plot, :), 'Color', [myColors(i, :), 0.5], 'LineWidth', 0.1);
                plot(thisTime, mean(tmpSpikes, 1), 'Color', [0, 0, 0], 'LineWidth', 2);
            end
            xlabel('Time (ms)');
            ylabel('Voltage (\muV)');
            
            %% decision
            
            % decide whether you want to analyze this cluster later on
            if ~isempty(c4a)
                
                % look up previous decision
                logIdx          = cell2mat({c4a.Cluster4Analysis.cluster}) == iClus;
                prevDecision    = c4a.Cluster4Analysis(logIdx).decision;
                
                % recheck or do not recheck the previous decision
                if param.bRecheck == false
                    
                    % use previous decision
                    myDecision 	= double(strcmp(prevDecision, 'yes'));
                    fprintf('%s, %s, cluster %d. Do you want to analyze this cluster? (1 = yes), previous decision was "%s".\n', ...
                        subjects{iSub}, chans(iChan).name, iClus, prevDecision)
                else
                    
                    % recheck previous decision
                    myDecision 	= input(sprintf('%s, %s, cluster %d. Do you want to analyze this cluster? (1 = yes), previous decision was "%s": ', ...
                        subjects{iSub}, chans(iChan).name, iClus, prevDecision));
                end
            else
                
                % make decision whether this cluster shall be analyzed
                myDecision      = input(sprintf('%s, %s, cluster %d. Do you want to analyze this cluster? (1 = yes): ', ...
                    subjects{iSub}, chans(iChan).name, iClus));
            end
            drawnow;
            
            % convert decision to string
            if myDecision == 1
                thisCluster4Analysis.decision   = 'yes'; % yes, use in future analyses
            else
                myDecision = 0;
                thisCluster4Analysis.decision   = 'no'; % no, do not use in future analyses
            end
            
            % concatenate informaton about different clusters
            Cluster4Analysis    = cat(1, Cluster4Analysis, thisCluster4Analysis);
            
            % bookkeeping
            allUnits            = cat(1, allUnits, [iSub, iChan, iClus, myDecision]);
            
            % close figures
            close all;
        end
        
        % save relevant output
        save(outputFile, 'Cluster4Analysis');
    end    
end

%% report

% number of cells
fprintf('\nCell count:\n');
fprintf('Number of cells before selection: %d.\n', size(allUnits, 1));
fprintf('Number of cells for analysis: %d.\n', sum(allUnits(:, 4)));
