function out = LK_RippleDetection_Preprocessing_20220319(mycfg, myeeg)
%
% LK_RippleDetection_Preprocessing_20220319 preprocesses the data for
% subsequent ripple detection.
%
% Preprocessing includes:
% - band-pass filtering;
% - envelope estimation;
% - envelope smoothing;
% - artifact removal from the envelope signal.
%
% Use as: out = LK_RippleDetection_Preprocessing_20220319(mycfg, myeeg);
%
% Input #1 is a structure with parameter settings including:
%   bpfilter        --> whether to use a bandpass filter;
%   bpfilttype      --> type of bandpass filter;
%   bpfreq          --> frequency range for the bandpass filter;
%   envSmoothTime   --> time period for smoothing the envelope (seconds).
%
% Input #2 is a fieldtrip structure with fields:
%   label           --> channel name;
%   fsample         --> sampling frequency (Hz);
%   trial           --> data;
%   time            --> time.
%
% Output is a structure with fields:
%   label           --> channel name;
%   fsample         --> sampling frequency;
%   trial           --> original data;
%   time            --> timepoints of the data;
%   index           --> sample indices of the data;
%   dataB           --> bandpass-filtered data;
%   dataEB          --> smoothed envelope of the bandpass-filtered data,
%                       where artifacts are removed.
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Preprocessing for ripple detection (channel: %s) ...\n', myeeg.label{:});

%% check input

% check whether artifact information is available
if ~isfield(mycfg, 'bArtifact')
    mycfg.bArtifact = false(size(myeeg.trial{1}));
    warning('No artifact information available.');
end

% check whether artifact information has the same size as the data
if size(mycfg.bArtifact, 2) ~= size(myeeg.trial{1}, 2)
    error('Artifact information has a different size than the data.');
end

%% band-pass filter

% report
fprintf('\n------------------------------------------------------ Band-pass filtering (%d to %d Hz) ...\n', ...
    min(mycfg.bpfreq), max(mycfg.bpfreq));

% bandpass-filter the data in the ripple band
cfg             = [];
cfg.bpfilter    = mycfg.bpfilter;
cfg.bpfilttype  = mycfg.bpfilttype;
cfg.bpfreq      = mycfg.bpfreq;
bpData          = ft_preprocessing(cfg, myeeg);

% band-pass data and its associated time and sample indices
dataB           = bpData.trial{1}; % bandpass
time            = bpData.time{1}; % timepoints
index           = 1:size(time, 2); % samples

%% envelope signal

% report
fprintf('\n------------------------------------------------------ Analytic envelope ...\n');

% compute envelope signal as instantaneous analytic amplitude
dataEB  = double(abs(hilbert(dataB - mean(dataB, 2, 'omitnan'))));

%% smooth the envelope signal

% report
fprintf('\n------------------------------------------------------ Envelope smoothing (%.1f ms) ...\n', ...
    1000 * mycfg.envSmoothTime);

% smooth the data
envSmoothSamples        = myeeg.fsample * mycfg.envSmoothTime; % smoothing window in samples
if mod(envSmoothSamples, 2) == 0
    envSmoothSamples    = envSmoothSamples + 1; % ensure that smoothing is centered on the mid time point
end
dataEB  = smoothdata(dataEB, 2, 'movmean', envSmoothSamples, 'omitnan'); % ignore nans

%% remove artifacts

% report
fprintf('\n------------------------------------------------------ Artifact removal (%.2f%% of the data) ...\n', ...
    100 * sum(mycfg.bArtifact) / numel(mycfg.bArtifact));

% set artifact samples to nan
dataEB(mycfg.bArtifact) = nan;

%% create output

% output structure
out             = [];

% data
out.label       = myeeg.label;
out.fsample     = myeeg.fsample;
out.trial       = myeeg.trial;
out.time        = {time};
out.index       = {index};
out.dataB       = {dataB};
out.dataEB      = {dataEB};
out.sampleinfo  = bpData.sampleinfo;

% report
fprintf('\n====================================================== Done: preprocessing for ripple detection.\n');
