function out = LK_RippleDetection_ExtractData_20220319(mycfg, myeeg)
%
% LK_RippleDetection_ExtractData_20220319 extracts LFP data around ripples.
%
% Use as: out = LK_RippleDetection_ExtractData_20220319(mycfg, myeeg);
%
% Input #1 is a structure with multiple parameter settings.
%
% Input #2 is a fieldtrip-like structure with fields
%   label    	--> channel name;
%   fsample    	--> sampling frequency;
%   trial     	--> original data;
%   time       	--> timepoints of the data;
%   index     	--> sample indices of the data;
%   dataB      	--> bandpass-filtered data;
%   dataEB    	--> smoothed envelope of the bandpass-filtered data, with
%                   nans during artifacts.
%
% Output is a fieldtrip structure that contains the data around ripples.
%   bRipple     --> Boolean indicating which timepoints are ripples
%   eegRipples  --> raw LFP data around each ripple
%   ripples     --> structure with information about all ripples
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Extraction of data around ripples (channel: %s) ...\n', myeeg.label{:});

%% ripples

% logical vector indicating the presence of ripples over the entire time
startEndIdx = [cell2mat({mycfg.ripples.startIdx}'), cell2mat({mycfg.ripples.endIdx}')];
bRipple     = LK_idx2log(startEndIdx, size(myeeg.time{1}, 2));

% eeg data from ripples
if size(mycfg.ripples, 1) > 0
    cfg         = [];
    cfg.trl     = [cell2mat({mycfg.ripples.peakIdx}') + round(min(mycfg.twoi4Data) * myeeg.fsample), ...
        cell2mat({mycfg.ripples.peakIdx}') + round(max(mycfg.twoi4Data) * myeeg.fsample), ...
        ones(size(mycfg.ripples, 1), 1) * round(min(mycfg.twoi4Data) * myeeg.fsample)];
    eegRipples  = ft_redefinetrial(cfg, myeeg); % fieldtrip adds nans if trial indices go beyond the data
else
    eegRipples  = [];
end

%% create output

% output structure
out             = [];

% ripples
out.ripples     = mycfg.ripples; % information about each ripple
out.bRipple     = bRipple; % true if ripple is present, otherwise false
out.eegRipples  = eegRipples; % EEG data for each ripple

% report
fprintf('====================================================== Done: extraction of data.\n');
