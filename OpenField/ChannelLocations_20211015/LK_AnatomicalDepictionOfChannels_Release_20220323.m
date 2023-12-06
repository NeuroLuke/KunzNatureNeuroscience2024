%==========================================================================
% This script:
% - computes an average MRI scan from all subject-specific MRIs.
% - creates a probability map of channels from a specific region.
% - creates Pylocator marker files to display the elctrode positions of
%   select electrodes in Pylocator. Channels of interest can be from, for
%   example: AMY, EC, HC, PHC, TP.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;

% paths
paths                   = [];
paths.spm               = 'E:\SPM\spm12\';
paths.implantations     = 'E:\Implantations_20201214\';
paths.subjects          = 'E:\OpenField\SubjectData_20210616\'; % subject information
paths.save              = 'E:\OpenField\ChannelLocations_20211015\20220323\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% add paths
addpath(genpath(paths.spm));

% settings
param                   = [];
param.region            = 'AMY'; % AMY, EC, HC, PHC, TP
param.choi              = {[param.region, 'leftBipolar']; [param.region, 'rightBipolar']}; % channels of interest
param.pylocatorSaveFile = [param.region, '_pylocator_bipolarChannels_20220323.txt'];
param.MRISaveFile       = 'meanPostMRI_20220323.nii';
param.channelSize       = 1; % spatial extent of the MNI coordinates (mm)
param.MRIChanSaveFile   = [param.region, '_chansOnMRIScan_percentBipolarChannels_extent', num2str(param.channelSize), '_20220323.nii'];

% subjects
subjects                = load(strcat(paths.subjects, 'subjects.mat'));
subjects                = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% preallocate
allMNIBipolar           = []; % MNI coordinates of all bipolar channels
allPylocator            = {}; % pylocator information for all channels

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    %% subject information
    
    % load subject information
    subjectFile = strcat(paths.subjects, subjects{iSub}, '\subjectdata_20210616.mat');
    subjectdata = load(subjectFile);
    subjectdata = subjectdata.subjectdata;
    
    % report
    fprintf('\nLoading subjectdata from "%s".\n', subjectFile);
    
    %% MNI coordinates of the channels
    
    % load pylocator file of this subject
    pylocatorFile   = dir(strcat(paths.implantations, subjectdata.implantation, '\Pylocator\i*'));
    fprintf('Name of the pylocator file: "%s".\n', pylocatorFile.name);
    pylocatorData   = importdata(fullfile(pylocatorFile.folder, pylocatorFile.name));
    chanNames       = pylocatorData.textdata;
    chanMNI         = pylocatorData.data(:, 1:3);
    
    % concatenate relevant channels from this subject
    chans           = [];
    for iChoi = 1:size(param.choi, 1)
        switch param.choi{iChoi}
            case 'AMYleftBipolar'
                chans   = cat(1, chans, subjectdata.macro.AMYleftBipolar);
            case 'AMYrightBipolar'
                chans   = cat(1, chans, subjectdata.macro.AMYrightBipolar);
            case 'ECleftBipolar'
                chans   = cat(1, chans, subjectdata.macro.ECleftBipolar);
            case 'ECrightBipolar'
                chans   = cat(1, chans, subjectdata.macro.ECrightBipolar);
            case 'HCleftBipolar'
                chans   = cat(1, chans, subjectdata.macro.HCleftBipolar);
            case 'HCrightBipolar'
                chans   = cat(1, chans, subjectdata.macro.HCrightBipolar);
            case 'PHCleftBipolar'
                chans   = cat(1, chans, subjectdata.macro.PHCleftBipolar);
            case 'PHCrightBipolar'
                chans   = cat(1, chans, subjectdata.macro.PHCrightBipolar);
            case 'TPleftBipolar'
                chans   = cat(1, chans, subjectdata.macro.TPleftBipolar);
            case 'TPrightBipolar'
                chans   = cat(1, chans, subjectdata.macro.TPrightBipolar);
        end
    end
    
    % loop through this subject's channels
    for iChan = 1:size(chans, 1)
    
        % channels that contribute to this bipolar channel
        contribChans    = split(chans{iChan}, '-');
        
        % MNI coordinate of both contributing channels
        MNI1            = chanMNI(strcmp(chanNames, contribChans{1}), :);
        MNI2            = chanMNI(strcmp(chanNames, contribChans{2}), :);
        
        % MNI coordinate of this bipolar channel
        MNIBipolar      = mean([MNI1; MNI2], 1); % average across the two contacts
        
        % collect across all bipolar channels
        allMNIBipolar   = cat(1, allMNIBipolar, MNIBipolar);
        
        % collect pylocator information across all bipolar channels
        allPylocator    = cat(1, allPylocator, ...
            {[subjects{iSub}, '_', chans{iChan}, ',', ...
            num2str(MNIBipolar(1), '%.6f'), ',', num2str(MNIBipolar(2), '%.6f'), ',', num2str(MNIBipolar(3), '%.6f'), ... % MNI coordinates of the channel
            ',2.0,1.0,1.0,1.0']}); % size and color of the channel
    end
    
    %% MRI scan
    
    % skip if this is not the first session of this subject (so that each
    % subject contributes only once to the average MRI scan)
    if iSub ~= 1 && strcmp(previousImplantation, subjectdata.implantation) == true
        fprintf('--- Skipping: %s.\n', subjectdata.name);
        continue;
    end
    
    % load post-MRI scan of this subject
    MRIFile     = dir(strcat(paths.implantations, subjectdata.implantation, '\Imaging\postMRI\*_inv.nii'));
    fprintf('Name of the MRI folder: "%s".\n', MRIFile.folder);
    V           = spm_vol(fullfile(MRIFile.folder, MRIFile.name));
    [Y, XYZ]    = spm_read_vols(V);
    
    % add this MRI image to all other MRI images
    if iSub == 1
        
        % summary MRI image
        sumMRI      = Y;
        
        % track number of subjects
        numSubjects = 1;
    else
        sumMRI      = sumMRI + Y;
        numSubjects = numSubjects + 1;
    end
    
    % track implantation date of this subject to compare with the next
    % implantation date
    previousImplantation    = subjectdata.implantation;    
end

%% write the bipolar channels into a Pylocator file

% report
fprintf('\nRESULTS:\n');
fprintf('Writing the pylocator information to a text file: "%s".\n', param.pylocatorSaveFile);
fprintf('The average y-coordinate is: %.3f.\n', mean(allMNIBipolar(:, 2), 'omitnan'));

% write pylocator file with MNI coordinates of the channels-of-interest
fileID = fopen(strcat(paths.save, param.pylocatorSaveFile), 'w');
for iP = 1:size(allPylocator, 1)
    fprintf(fileID, '%s\n', allPylocator{iP});
end
fclose(fileID);

%% process the grand-average MRI image

% divide the sum image by the number of subjects
fprintf('\nEstimating the average MRI scan of all subjects.\n');
fprintf('Number of subjects: %d.\n', numSubjects);
meanMRI     = sumMRI ./ numSubjects;

% write average image to file
meanV       = V;
meanV.fname = strcat(paths.save, param.MRISaveFile);
spm_write_vol(meanV, meanMRI);

%% probability map of all bipolar channels

% report
fprintf('\nProcessing all bipolar channels.\n');
fprintf('Number of bipolar channels: %d.\n', size(allMNIBipolar, 1));

% create image indicating the location of the bipolar channels
chanY       = zeros(size(meanMRI)); % data of all voxels
chanXYZ     = XYZ; % xyz coordinates of all voxels

% loop through channels
numChans    = 0; % number of channels
for iBipolarChan = 1:size(allMNIBipolar, 1)
    
    % identify all MNI coordinates within a range of X mm around the
    % channel's MNI coordinate
    D           = pdist2(allMNIBipolar(iBipolarChan, :), chanXYZ');
    thisY       = reshape(D <= param.channelSize, size(chanY, 1), size(chanY, 2), size(chanY, 3));
    
    % add this channel's area to all other channels' area
    chanY       = chanY + thisY;
    
    % keep track of the total number of channels
    numChans    = numChans + 1;
end

% average the bipolar channel areas across channels
meanY   = 100 .* chanY ./ numChans; % in percent

% write average image to file
meanV.fname = strcat(paths.save, param.MRIChanSaveFile);
spm_write_vol(meanV, meanY);
