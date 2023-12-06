%==========================================================================
% This script imports the microwire data.
%
% Lukas Kunz, 2020
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\NPMK-master\20201215\'));

% paths
beh_path    = 'E:\OpenField\Beh_20201215\';
data_path   = 'E:\OpenField\DataComplete_20180321\';
cut_path  	= 'E:\OpenField\Micro_20201215\';

% subjects
subjects = {...
    'FR006a', 111, 2678; ...
    'FR006b', 500, 3755; ...
    'FR007', 6337, 8200; ...
    'FR008a', 357, 3565; ...
    'FR008b', 122, 3437; ...
    'FR009a', 39, 3644; ...
    'FR009b', 63.02, 3112; ...
    'FR010', 10, 2500; ...
    'FR011', 246, 4674; ...
    'FR013', 178, 2033; ...
    'FR014', 377, 3170; ...
    'FR017', 800, 5500; ...
    'FR020', 621, 3633; ...
    'FR021', 882, 3907; ...
    'FR022', 71, 3888; ...
    'FR023a', 144, 3996; ...
    'FR023b', 202, 3370; ...
    'FR024', 281, 3362; ...
    'FR025a', 680, 3650; ...
    'FR025b', 686, 2850; ...
    'FR025c', 666, 3040; ...
    'FR026', 1148, 4370; ...
    'FR027', 2766, 5200; ...
    'FR028a', 1400, 10430; ...
    'FR028b', 970, 2800; ...
    'FR029', 640, 4700; ...
    'FR030', 758, 4000; ...
    };

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nWorking on subject %s ...\n', subjects{iSub, 1});
    
    % read header of file
    files 	= dir(strcat(data_path, subjects{iSub, 1}, filesep, 'Micro', filesep, '*.ns5'));
    hdr   	= openNSx(fullfile(files(1).folder, files(1).name), 't:1:2', 'sample', 'uV');
    wireidx = 1:length(hdr.ElectrodesInfo);
    chanidx = wireidx(cellfun(@isempty, regexp({hdr.ElectrodesInfo.Label}, 'ainp')));
    
    % sampling frequency for wave-clus
    sr      = hdr.MetaTags.SamplingFreq;
    
    %% check trigger data
    
    % get behavioral trigger data
    tmp         = load(strcat(beh_path, subjects{iSub, 1}, filesep, 'triggers.mat'));
    triggers    = tmp.triggers;
    
    % get microwire trigger data
    rmtFile     = strcat(data_path, subjects{iSub, 1}, filesep, 'Micro', filesep, 'reconstructed_microwire_trigger.mat');
    rmt         = load(rmtFile);
    trigdata    = rmt.reconstructed_microwire_trigger;
    
    % plot trigger data
    figure('units', 'normalized', 'position', [0.05, 0.025, 0.9, 0.8]);
    plot([1:20:length(trigdata) - 1] ./ sr, trigdata(1:20:end));
    xlabel('Time (sec)');
    ylabel('Voltage');
    title({strrep(['Trigger data for ', subjects{iSub, 1}], '_', '\_'), ...
        ['Number of triggers = ', num2str(size(triggers, 1))], ...
        ['Length of microwire data = ', num2str(length(trigdata) / sr), ' sec'], ...
        ['Time between first and last trigger onset = ', num2str(range(triggers)), ' sec']});
    drawnow;
    
    % determine borders for the data
    if ~isempty(subjects{iSub, 2})
        % use previously determined new start and new end
        newStart    = subjects{iSub, 2};
        newEnd      = subjects{iSub, 3};
        fprintf('\nNew start and end of data: %.3f to %.3f sec.\n\n', newStart, newEnd);
    else
        % decide whether you want to continue
        bProcess    = input('Do you want to process this subject? (0 = no) ');
        if bProcess == 0
            continue;
        end
        
        % define new start and end of data
        newStart    = input('New starting point of data [sec]: ');
        newEnd      = input('New ending point of data [sec]: ');
    end
    
    %% load entire data and cut
    
    % create folder and save header information
    mkdir(strcat(cut_path, subjects{iSub, 1}));
    save(strcat(cut_path, subjects{iSub, 1}, filesep, 'original_hdr_1stfile'), 'hdr');
    
    % loop through wires
    for iwire = 1:max(wireidx)
        
        % report
        fprintf('\tWorking on wire %d ...\n', iwire);
        
        % load data of this channel from all available files (max. 2)
        if size(files, 1) == 1
            
            % load data
            chan        = openNSx(fullfile(files.folder, files.name), ['c:', num2str(wireidx(iwire))], 'uV'); % import as microvolt
            chandata    = chan.Data;
                        
        elseif size(files, 1) == 2
            
            % load data
            chan1   = openNSx(fullfile(files(1).folder, files(1).name), ['c:', num2str(wireidx(iwire))], 'uV'); % import as microvolt
            chan2   = openNSx(fullfile(files(2).folder, files(2).name), ['c:', num2str(wireidx(iwire))], 'uV');
            
            % fill small gap between the two files if necessary
            time2add    = sum(chan2.MetaTags.DateTimeRaw(5:8) .* [3600, 60, 1, 0.001]) - (sum(chan1.MetaTags.DateTimeRaw(5:8) .* [3600, 60, 1, 0.001]) + chan1.MetaTags.DataPointsSec);
            samples2add = round(time2add * sr);
            
            % combine the files and the gap data
            chandata    = [chan1.Data, repmat(chan1.Data(end), 1, samples2add), chan2.Data];
        end
        
        % double-check length of the data
        if numel(chandata) ~= numel(trigdata)
            error('The channel data has a different length than the trigger data.\n');
        end
        
        % trigger data
        if ~isempty(regexp(chan.ElectrodesInfo.Label, 'ainp2', 'once'))
            chandata    = trigdata;
        end
        
        % cut data
        data    = chandata(newStart * sr : newEnd * sr);
    
        % save data and sampling rate
        mkdir(strcat(cut_path, subjects{iSub, 1}, filesep, hdr.ElectrodesInfo(iwire).Label));
        save(strcat(cut_path, subjects{iSub, 1}, filesep, hdr.ElectrodesInfo(iwire).Label, filesep, 'datacut'), 'data', 'sr', '-v7.3');
        
        % close open files
        fclose all;
        
        %% check alignment between behavioral and microwire triggers
        
        % double-check trigger
        if ~isempty(regexp(chan.ElectrodesInfo.Label, 'ainp2', 'once'))
            
            % behavioral triggers
            behtrig_timepoints  = triggers;
            
            % microwire trigger data
            microtrig_data     	= transpose(double(data));
            microtrig_time    	= transpose((0:length(microtrig_data) - 1) ./ sr); % start at zero (sec)
    
            % microwire trigger timepoints
            microtrig_data          = conv(microtrig_data, ones(10, 1) ./ 10, 'same'); % smooth to avoid tiny bumps
            microtrig_timepoints    = microtrig_time(diff(microtrig_data > 2000) == 1);
                        
            % report
            rhoAlignment    = corr(diff(microtrig_timepoints), diff(behtrig_timepoints), 'type', 'spearman');
            fprintf('Trigger information:\n');
            fprintf('The ITI-correlation between microtriggers and behtriggers is rho = %.6f.\n', rhoAlignment);
            fprintf('There are %d triggers in total.\n', length(microtrig_timepoints));
            if rhoAlignment < 0.999
                input('"rhoAlignment" is below 0.999. Do you want to continue? (Press any key to continue with the next subject.)');
            end
        end        
    end
    
    % close open figures
    close all;
end

%% end
fprintf('\nCompleted.\n');
