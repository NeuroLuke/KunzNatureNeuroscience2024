%==========================================================================
% This script imports the iEEG-data.
%--------------------------------------------------------------------------
% This results in the fieldtrip structure "eeg_original" with
%   .label:     cell-array containing strings, Nchan x 1
%   .fsample:   sampling frequency in Hz (single number)
%   .trial:     cell-array containing a data matrix for each trial (1 x 
%               Ntrial), each data matrix is a Nchan x Nsamples matrix
%   .time:      cell-array containing a time axis for each trial (1 x
%               Ntrial), each time axis is a 1 x Nsamples vector
%==========================================================================

%% settings

% start
clear; clc; close all;

% paths
behPath     = 'E:\OpenField\Beh_20201215\';
dataPath    = 'E:\OpenField\DataComplete_20180321\';
savePath    = 'E:\OpenField\Macro_20201215\';

% subjects
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

%% loop through subjects and import eeg data
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nWorking on subject: %s ...\n', subjects{iSub, 1});
    
    % get files that contain eeg data
    files0      = dir(strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro\EEGData*.0*'));
    files1      = dir(strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro\EEGData*.1*'));
    realFiles   = [files0; files1];
    
    % preallocate with example file
    fid         = fopen(fullfile(realFiles(1).folder, realFiles(1).name));
    tmp         = fread(fid, inf, 'int32');
    fclose(fid);
    eegdata     = nan(size(realFiles, 1), size(tmp, 1));
    
    % loop through files and load them
    for iFile = 1:size(realFiles, 1)
        
        tmpFilename     = fopen(fullfile(realFiles(iFile).folder, realFiles(iFile).name));
        fprintf('Reading "%s" ...\n', fullfile(realFiles(iFile).folder, realFiles(iFile).name));
        
        % load raw data
        tmpEegdata     	= fread(tmpFilename, inf, 'int32');
        
        % subtract value of first data point so that it always starts at 0
        eegdata(iFile, :)   = bsxfun(@minus, tmpEegdata, tmpEegdata(1));
        
        % close file
        fclose(tmpFilename);
    end
    
    % get the multiplication factor from the ini-file
    iniFile         = dir(strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro\EEGData.ini'));
    if isempty(iniFile)
        multFactor  = 7.960022e-008 * 1000000; % uV per Bit
        fprintf('\nWARNING: Using the standard bit-to-uV conversion factor.\n');
    else
        iniFid      = fopen(fullfile(iniFile.folder, iniFile.name));
        iniData     = textscan(iniFid, '%s');
        bitInfo     = iniData{1}(contains(iniData{1}, 'Bit='));
        bitInfo     = split(bitInfo, '=');
        multFactor  = str2double(bitInfo{2}) * 1000000; % uV per Bit
    end
    
    % sanity check for the multiplication factor
    if multFactor ~= 7.960022e-008 * 1000000
        error('Problem with "multFactor". Check the "EEGData.ini" file of this subject.\n');
    end
    
    % multiply original eegdata with multiplication factor
    eegdata_uV = bsxfun(@times, eegdata, multFactor);
    
    % free memory
    clear eegdata;
    
    % import channel-info from the channel-info file
    tmp         = importdata(strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro', filesep, 'channel-info.txt'));
    elecLabel   = cell(size(tmp, 1), 1);
    for iLine = 1:size(tmp, 1)
        myString            = strsplit(regexprep(tmp{iLine, :}, '\t', ','), ',');
        elecLabel{iLine, 1} = myString{2};
    end
    
    % get the sampling frequency from the cnv-info file
    cnvFile     = dir(strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro\cnv-info.txt'));
    cnvFid      = fopen(fullfile(cnvFile.folder, cnvFile.name));
    cnvData     = textscan(cnvFid, '%s');
    sampleInfo  = cnvData{1}(~cellfun(@isempty, (regexp(cnvData{1}, '^sample_rate='))));
    sampleInfo  = split(sampleInfo, '=');    
    fsample     = str2double(sampleInfo{2});
    
    % sanity check for the sampling rate
    if fsample ~= 2000
        error('Sampling frequency is not 2 kHz. Check this subject.\n');
    end
    
    % insert data into fieldtrip structure
    eeg_original.label      = elecLabel;
    eeg_original.fsample    = fsample;
    eeg_original.trial      = {eegdata_uV}; % per trial: channel x time
    eeg_original.time       = {0 : 1/fsample : (size(eegdata_uV, 2) - 1) / fsample};
    
    %% trigger data
    
    % use adjusted ("reconstructed") trigger data for easier alignment
    retFile     = strcat(dataPath, subjects{iSub, 1}, filesep, 'Macro\reconstructed_eeg_trigger.mat');        
    ret         = load(retFile);
    trigIdx     = strcmp(eeg_original.label, 'Trigger');
    eeg_original.trial{1}(trigIdx, :)   = ret.reconstructed_eeg_trigger;
    
    %% double-check alignment between behavioral and neural triggers
    
    % get behavioral data
    tmp                 = load(strcat(behPath, subjects{iSub, 1}, filesep, 'triggers.mat'));
    behtrig_timepoints  = tmp.triggers(~isnan(tmp.triggers));
    
    % plot trigger data
    figure('units', 'normalized', 'position', [0.025, 0.025, 0.95, 0.8]);
    plot(eeg_original.time{1}, eeg_original.trial{1}(strcmp(eeg_original.label, 'Trigger'), :));
    xlabel('Time (sec)');
    ylabel('Voltage');
    title({strrep(['Trigger data for ', subjects{iSub, 1}], '_', '\_'), ...
        ['Number of behavioral triggers = ', num2str(size(behtrig_timepoints, 1))], ...
        ['Length of trigger-data = ', num2str(range(eeg_original.time{1})), ' sec'], ...
        ['Time between first and last cue-onset = ', num2str(range(behtrig_timepoints)), ' sec']});
    drawnow;
    
    % eeg-trigger timepoints (smooth to get rid of small bumps)
    tmpTrig          	= smoothdata(eeg_original.trial{1}(strcmp(eeg_original.label, 'Trigger'), :), 2, 'movmean', 10);
    tmpTrig           	= tmpTrig - min(tmpTrig);
    tmpTrig           	= tmpTrig ./ max(tmpTrig);
    eegtrig_timepoints  = transpose(eeg_original.time{1}(diff(tmpTrig > 0.5) == 1));
    
    % correlation between inter-trigger intervals
    [rho, pval]         = corr(diff(behtrig_timepoints), diff(eegtrig_timepoints), 'type', 'spearman');
    fprintf('%s: Correlation between eeg and behavioral triggers: rho = %.6f, p = %.6f.\n', subjects{iSub, 1}, rho, pval);
    
    %% save data
    
    % save
    mkdir(strcat(savePath, subjects{iSub, 1}));
    save(strcat(savePath, subjects{iSub, 1}, filesep, 'eeg_original'), 'eeg_original', '-v7.3');
    
    % free memory
    clear eegdata_uV eeg_original;
    
    % close figures
    close all;
end

%% report

fprintf('\nProcessing completed.\n');
