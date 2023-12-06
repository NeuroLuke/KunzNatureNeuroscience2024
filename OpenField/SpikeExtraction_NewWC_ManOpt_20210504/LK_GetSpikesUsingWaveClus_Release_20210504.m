%==========================================================================
% This script performs spike-detection and -clustering using wave-clus.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('C:\WaveClus3\')); % Chaure et al., 2018

% paths
paths           = [];
paths.data      = 'E:\OpenField\Micro_20201215\';
paths.spike     = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\';
paths.pics      = strcat(paths.spike, 'AllPics\');
if ~exist(paths.pics, 'dir')
    mkdir(paths.pics);
end

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

% bookkeeping: wires for which wave clus cannot extract any spikes
allSubjNoWC     = cell(size(subjects, 1), 1);

%% loop through subjects
parfor iSub = 1:size(subjects, 1)
    
    % set rng for reproducibility
    rng(1, 'twister');
    
    % report
    fprintf('\nWorking on subject %s ...\n', subjects{iSub});
    cd(paths.spike);
    
    % bookkeeping
    noWC    = [];
    
    %% list wires
    
    % trigger channels
    trigs   = dir(strcat(paths.data, subjects{iSub}, filesep, 'ainp*'));
    
    % data channels
    wires   = dir(strcat(paths.data, subjects{iSub}, filesep, 'chan*'));
    tmp     = split(transpose({wires.name}), 'n');
    [B, I]  = sort(cellfun(@str2num, tmp(:, 2)));
    wires   = wires(I);
    
    % combine trigger- and wire-files
    files   = [trigs; wires];
    
    % report
    fprintf('\n------------ Performing wave_clus\n');
    
    %% loop through files and extract action potentials
    for iFile = 1:size(files, 1)
        
        % report
        fprintf('\tWorking on file "%s" ...\n', files(iFile).name);
        
        % storage destination
        thisSavePath    = fullfile(paths.spike, subjects{iSub}, files(iFile).name, filesep);
        
        % copy previously cut data to storage destination
        mkdir(thisSavePath);
        copyfile(fullfile(files(iFile).folder, files(iFile).name, filesep, 'datacut.mat'), thisSavePath);
        
        %% wave clus 3
        % needs input file containing "data" and "sr"
        
        % perform wave clus
        if ~isempty(regexp(files(iFile).name, 'chan', 'once'))
            
            % go to directory with data of this file
            cd(thisSavePath);
            
            % extract spikes and do clustering
            try
                % perform spike extraction
                file2use4wc             = {strcat(thisSavePath, 'datacut.mat')};
                myPar                   = [];
                myPar.detection         = 'neg'; % detection type
                myPar.randomseed        = 1; % random seed
                Get_spikes(file2use4wc, 'par', myPar);
                
                % perform spike clustering
                file2use4wc             = {strcat(thisSavePath, 'datacut_spikes.mat')};
                myPar                   = [];
                myPar.randomseed        = 1; % random seed
                myPar.template_sdnum    = 1.5; % for each unsorted spike, look up the closest template - but only if it is within some range
                myPar.min_clus          = 60; % minimum size of a cluster (default 20)
                myPar.max_clus        	= 10; % maximum number of clusters allowed (default 13)
                myPar.mintemp           = 0.05; % minimum temperature
                Do_clustering(file2use4wc, 'par', myPar);
                
                % copy output picture(s) to overview folder
                pics    = dir(fullfile(thisSavePath, '*.png'));
                for iPic = 1:length(pics)
                    copyfile(fullfile(pics(iPic).folder, pics(iPic).name), ...
                        strcat(paths.pics, subjects{iSub}, '_', files(iFile).name, '_', pics(iPic).name));
                end
            catch
                warning('Wave-clus could not be performed.');
                noWC    = cat(1, noWC, {subjects{iSub}, files(iFile).name});
            end
            
            %% save space
            
            % delete original input file
            delete(strcat(thisSavePath, 'datacut.mat'));
        end
    end
    
    % collect sanity-check information across subjects
    allSubjNoWC{iSub, 1}    = noWC;
end