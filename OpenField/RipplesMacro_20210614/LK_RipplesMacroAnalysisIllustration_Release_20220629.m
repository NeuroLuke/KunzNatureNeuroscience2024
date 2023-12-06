%==========================================================================
% This script illustrates the procedure for identifying ripples in the
% macro LFP data.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));
addpath(genpath('E:\WaveClus3\wave_clus-master'));

% fieldtrip
addpath('E:\fieldtrip\fieldtrip-20210614\');
ft_defaults;

% paths
paths           = [];
paths.macro     = 'E:\OpenField\MacroPreprocessing_20210627\'; % macro LFP
paths.ripple    = 'E:\OpenField\RipplesMacro_20210614\20220320_FBS_80to140Hz_tMF2\'; % ripples
paths.micro     = 'E:\OpenField\Micro_20201215\'; % micro LFP
paths.spikes    = 'E:\OpenField\SpikeExtraction_NewWC_ManOpt_20210504\'; % spiking
paths.save      = 'E:\OpenField\RipplesMacro_20210614\Illustration_20220607\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% settings
param                   = [];
param.subject           = 'Freiburg_20190130'; % example subject
param.channel           = 'HAR1-HAR2'; % example macro LFP channel
param.microChannel      = 'chan11'; % example microwire channel
param.threshMinFac      = 2; % minimum threshold for ripples
param.threshMaxFac      = Inf; % maximum threshold for ripples
param.threshPeakMinFac  = 3; % minimum peak threshold for ripples

% specific settings
param.plot.rippleIdx    = 477; % example ripple
param.plot.toi          = [-2, 2]; % time limits around the ripple
param.plot.colorbarLim  = [0, 7]; % colorbar limits for power spectrogram
param.plot.freqTicks    = [80, 140];
param.plot.microYLim    = [-120, 45]; % y-limits for the microwire data
param.plot.fontUnits    = 'centimeters';
param.plot.fontSize     = 0.4;
param.plot.saveName     = strcat(param.subject, '_', param.channel, '_', param.microChannel, '_IllustrationRippleDetection');

% save the settings
save(strcat(paths.save, param.subject, '_', param.channel, '_', param.microChannel, '_settings'));

%% iEEG data

% load data for ripples - contains raw LFP, bandpass-filtered LFP, and the
% envelope of the bandpass-filtered LFP (artifacts removed)
data4Ripples    = load(strcat(paths.ripple, param.subject, '\', param.channel, '\data4Ripples.mat'));
data4Ripples    = data4Ripples.data4Ripples;

% load ripples - contains ripple times and ripple samples
ripples         = load(strcat(paths.ripple, param.subject, '\', param.channel, '\ripples.mat'));
ripples         = ripples.ripples;

% load artifacts - contains artifact times and artifact samples
artifacts       = load(strcat(paths.ripple, param.subject, '\', param.channel, '\artifacts.mat'));
artifacts       = artifacts.artifacts;

%% microwire data

% load spiking data
t               = load(strcat(paths.spikes, param.subject, '\', param.microChannel, '\times_datacut.mat'));
spikesClass     = t.cluster_class(:, 1);
spikesTime      = t.cluster_class(:, 2) ./ 1000; % convert spike times to seconds

% load microwire LFP
datacut         = load(strcat(paths.micro, param.subject, '\', param.microChannel, '\datacut.mat'));
microTime       = ((1:numel(datacut.data)) - 1) ./ datacut.sr;
microLFP        = datacut.data;

% filter between 300 and 3,000 Hz using wave-clus
par                 = [];
par.sr              = datacut.sr;
par.detect_fmin     = t.par.detect_fmin;
par.detect_fmax     = t.par.detect_fmax;
par.detect_order    = t.par.detect_order;
microLFPFiltered    = spike_detection_filter(microLFP, par);

% calculate amplitude threshold for detecting spikes
spikesThresh        = -1 * t.par.stdmin * median(abs(microLFPFiltered)) / 0.6745; % MAD approximation of the STD

%% align microwire timeline to macro iEEG timeline

% data for macro triggers
macroTrig           = load(strcat(paths.macro, param.subject, '\Trigger.mat'));
macroTrigTime       = macroTrig.Trigger.time{1};
macroTrigData       = macroTrig.Trigger.trial{1};
macroTrigData       = smoothdata(macroTrigData, 2, 'movmean', 10); % smooth to get rid of tiny bumps

% macro-trigger timepoints
thisThresh          = max(macroTrigData) / 2; % threshold for detecting triggers
logIdx              = diff(macroTrigData > thisThresh) == 1;
macroTrigTimepoints = transpose(macroTrigTime(logIdx));

% data for micro triggers
microTrig           = load(strcat(paths.micro, param.subject, '\ainp2\datacut.mat')); % trigger channel
microTrigData       = microTrig.data;
microTrigData       = smoothdata(microTrigData, 2, 'movmean', 10); % smooth to get rid of tiny bumps
microTrigTime       = ((1:size(microTrigData, 2)) - 1) ./ microTrig.sr; % first sample has the timepoint 0

% micro triggers and inter-trigger-intervals
thisThresh          = max(microTrigData) / 2; % threshold for detecting triggers
logIdx              = diff(microTrigData > thisThresh) == 1;
microTrigTimepoints = transpose(microTrigTime(logIdx));

% convert microtime into macrotime
P                           = polyfit(microTrigTimepoints, macroTrigTimepoints, 1);
microTrigTimeInMacroTime    = P(1) .* microTrigTime + P(2);
spikesTimeInMacroTime       = P(1) .* spikesTime + P(2);

%% macro LFP processing

% apply low high-pass filter to get rid of slow drifts
cfg             = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.hpfiltord   = 5;
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 800;
rawLFP          = ft_preprocessing(cfg, data4Ripples);

% macro LFP
macroTime       = rawLFP.time{1};
macroLFP        = rawLFP.trial{1};

%% macro ripple band processing

% ripple band
macroRippleBand = data4Ripples.dataB{1};

% estimate threshold (settings are from the ripple analysis)
macroEnvelope   = data4Ripples.dataEB{1};
center          = mean(macroEnvelope, 2, 'omitnan'); % center of the envelope
spread          = std(macroEnvelope, [], 2, 'omitnan'); % spread of the envelope
threshMin       = center + spread .* param.threshMinFac; % minimum threshold
threshMax       = center + spread .* param.threshMaxFac; % maximum threshold
threshPeakMin   = center + spread .* param.threshPeakMinFac; % minimum threshold for peaks

%% macro power spectrogram processing

% time-frequency power spectrogram of the entire data
cfg             = [];
cfg.method      = 'wavelet';
cfg.width       = 7;
cfg.toi         = macroTime;
cfg.foi         = 30:2:190;
cfg.pad         = 'nextpow2';
cfg.output      = 'pow';
tfData          = ft_freqanalysis(cfg, data4Ripples);

% power at all time points
macroPow        = squeeze(tfData.powspctrm); % freq x time

% grand-average power (average across time)
macroMeanPow    = mean(macroPow, 2, 'omitnan'); % freq x 1

% normalize power by dividing by grand-average power
macroNormPow    = macroPow ./ repmat(macroMeanPow, 1, size(macroPow, 2));

%% figure

% create figure
dt                          = [];
dt.param                    = param;
dt.macroTime                = macroTime;
dt.macroLFP                 = macroLFP;
dt.macroRippleBand          = macroRippleBand;
dt.macroEnvelope            = macroEnvelope;
dt.threshMin                = threshMin;
dt.threshPeakMin            = threshPeakMin;
dt.ripples                  = ripples;
dt.artifacts                = artifacts;
dt.tfData                   = tfData;
dt.macroNormPow             = macroNormPow;
dt.microTrigTimeInMacroTime = microTrigTimeInMacroTime;
dt.microLFPFiltered         = microLFPFiltered;
dt.spikesThresh             = spikesThresh;                      
dt.spikesClass              = spikesClass;
dt.spikesTimeInMacroTime    = spikesTimeInMacroTime;
dt.time2Plot                = round(ripples.ripples(param.plot.rippleIdx).peakTime + param.plot.toi, 1); % time around example ripple
f = LK_PlotRipplesMacroAnalysisIllustration_20231030(dt);

% save figure
set(f, 'InvertHardcopy', 'off');
LK_print(f, strcat(paths.save, param.plot.saveName), '-dtiff', '-r600');
