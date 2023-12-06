%==========================================================================
% This script extracts behavioral information from the logfiles.
%==========================================================================

% start
clear; clc; close all;

% paths
behPath     = 'E:\OpenField\Beh_20201215\';
dataPath    = 'E:\OpenField\DataComplete_20180321\';

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

%% loop through subjects
for iSub = 1:length(subjects)
    
    % report progress
    fprintf('\n============================================================\n');
    fprintf('SUBJECT: %s\n', subjects{iSub, 1});
    
    % subject-specific paths
    dataPathSpec    = strcat(dataPath, subjects{iSub, 1}, filesep, 'Beh');
    behPathSpec     = strcat(behPath, subjects{iSub, 1});
    if ~exist(behPathSpec, 'dir')
        mkdir(behPathSpec);
    end
    
    %% read in data
    
    % open logfile
    logfile       	= dir(strcat(dataPathSpec, '\*.log'));
    fid             = fopen(fullfile(logfile.folder, logfile.name));
    textdata        = [];
    
    % initialize line index
    iline = 1;
    
    % read logfile line-by-line
    tline = fgetl(fid);
    while ischar(tline)
        
        % combine all text lines into one variable
        textdata{iline, 1} = tline;
        
        % increase line index
        iline = iline + 1;
        
        % read next line
        tline = fgetl(fid);
    end
    
    % close logfile
    fclose(fid);
    
    % report progress
    fprintf('You have read %d text lines...\n', iline);
    
    %% save scriptlog data
    
    % extract scriptlog data
    data    = textdata(contains(textdata, 'ScriptLog'), :);
    
    % save scriptlog data
    save([behPathSpec, filesep, 'data'], 'data');
    
    %% extract information about behavioral triggers
    
    % cue onsets are the trigger timepoints
    cueData     = textdata(contains(textdata, 'CUE_start'), :);
    splits      = split(cueData);
    triggers    = cellfun(@str2double, splits(:, 3));

    % save trigger information
    save([behPathSpec, filesep, 'triggers'], 'triggers');    
end

% end report
fprintf('\nExtracted behavioral information for all subjects.\n');
