%==========================================================================
% This script estimates the raw and relative power spectrum of ripples.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; clc; close all;

% paths
paths           = [];
paths.info      = 'E:\OpenField\SubjectData_20210616\';
paths.ripple    = 'E:\OpenField\RipplesMacro_20210614\20220320_FBS_80to140Hz_tMF2\';
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';
paths.functions = 'E:\OpenField\Functions\';
paths.save      = 'E:\OpenField\RipplesMacroPowerSpectrum_20220323\20230317\';
mkdir(paths.save);

% settings
param               = [];
param.myRNG         = 111;
param.tfr.method  	= 'wavelet';
param.tfr.foi       = [1:10, 12:2:30, 35:5:200];
param.tfr.toi       = 'all';
param.tfr.pad       = 'nextpow2';
param.tfr.width     = 7;
param.tfr.output    = 'pow';

% add paths
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
ft_warning('off');

% subjects
subjects    = load(strcat(paths.info, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

% preallocate
allRes      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report
    fprintf('\nSUBJECT: %s.\n\n', subjects{iSub});
    rng(param.myRNG); % reset rng for reproducibility
    
    % channels of interest on which you also identified the ripples
    choi = LK_dir(strcat(paths.ripple, subjects{iSub}, '\*'));
    
    %% loop through channels of interest
    for iChan = 1:numel(choi)
        
        % path for this channel
        chanPath        = strcat(choi(iChan).folder, '\', choi(iChan).name, '\');
        
        % load macro LFP data used for ripple detection
        d4R             = load(strcat(chanPath, 'data4Ripples.mat'));
        d4R             = d4R.data4Ripples;
        
        % load artifacts
        artifacts       = load(strcat(chanPath, 'artifacts.mat'));
        artifacts       = artifacts.artifacts;
        
        % load ripples
        ripples         = load(strcat(chanPath, 'ripples.mat'));
        ripples         = ripples.ripples;
        
        % load surrogate ripples
        ripplesSurro    = load(strcat(chanPath, 'ripplesSurro.mat'));
        ripplesSurro    = ripplesSurro.ripplesSurro;
                
        %% power across the entire time using wavelets
        
        % perform wavelet analysis using fieldtrip
        cfg             = param.tfr;
        TFR             = ft_freqanalysis(cfg, d4R);
        
        % extract power
        pow             = squeeze(TFR.powspctrm);
        
        % exclude artifacts
        logIdx          = artifacts.bArtifact;
        pow(:, logIdx)  = nan;
        
        % relative power (relative to grand average)
        meanPow         = mean(pow, 2, 'omitnan');
        normPow         = pow ./ repmat(meanPow, 1, size(pow, 2));
        
        %% power spectrum for ripples and surrogate ripples
        
        % raw power spectrum
        ripplesPow          = transpose(mean(pow(:, ripples.bRipple), 2, 'omitnan')); % 1 x freqs
        ripplesSurroPow     = transpose(mean(pow(:, ripplesSurro.bRipple), 2, 'omitnan'));
        
        % normalized power spectrum
        ripplesNormPow      = transpose(mean(normPow(:, ripples.bRipple), 2, 'omitnan')); % 1 x freqs
        ripplesSurroNormPow = transpose(mean(normPow(:, ripplesSurro.bRipple), 2, 'omitnan'));
        
        % sanity check
        if sum(ripples.bRipple) ~= sum(ripplesSurro.bRipple)
            error('Ripples and surrogate ripples do not have the same data length.');
        end
        
        %% collect results
        
        % results from this ripple channel
        thisRes                     = [];
        thisRes.idx                 = [iSub, iChan];
        thisRes.subject             = subjects{iSub};
        thisRes.channel             = choi(iChan).name;
        
        % mean power across the entire experiment
        thisRes.meanPow             = transpose(meanPow);
        
        % mean power during ripples
        thisRes.ripplesPow          = ripplesPow; % absolute power
        thisRes.ripplesSurroPow     = ripplesSurroPow;
        thisRes.ripplesNormPow      = ripplesNormPow; % relative power
        thisRes.ripplesSurroNormPow = ripplesSurroNormPow;
        
        % collect across ripple channels
        allRes                      = cat(1, allRes, thisRes);
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load results
r = load(strcat(paths.save, 'results.mat'));

%% figure comparing the raw power spectra

% unfold the power spectra
allRipplesPow       = cell2mat({r.allRes.ripplesPow}'); % channels x freqs
allRipplesSurroPow  = cell2mat({r.allRes.ripplesSurroPow}');

% different data groups
groups  = {'Surrogates', [0.7, 0.7, 0.7]; 'Ripples', [0, 0, 0]};
m       = {mean(allRipplesSurroPow, 1); mean(allRipplesPow, 1)};
ste     = {LK_ste(allRipplesSurroPow); LK_ste(allRipplesPow)};

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5.5, 6], 'Color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [1.65, 1.75, 3.5, 3.5]);
hold on;
for iM = 1:numel(m)
    patch([r.param.tfr.foi, fliplr(r.param.tfr.foi)], [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], groups{iM, 2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(r.param.tfr.foi, m{iM}, 'Color', groups{iM, 2});
end
t1 = text(0.09, 0.3, 'Empirical', 'units', 'normalized', 'Color', [0.1, 0.1, 0.1]);
t2 = text(0.09, 0.15, 'Surrogate', 'units', 'normalized', 'Color', [0.7, 0.7, 0.7]);
set(gca, ...
    'xtick', [1, 10, 100], 'ytick', [10^3, 10^6], ...
    'xscale', 'log', 'yscale', 'log', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Frequency (Hz)');
yl = ylabel('Power (\muV^{2})');
tl = title('Absolute power');
set([gca, t1, t2, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
% save figure
LK_print(f, strcat(paths.save, 'RawPowerRipplesVSSurrogateRipples_20211022'), '-dtiff', '-r300');

%% figure comparing the relative power spectra (relative to the entire experiment)

% unfold the power spectra
allRipplesNormPow       = cell2mat({r.allRes.ripplesNormPow}'); % channels x freqs
allRipplesSurroNormPow  = cell2mat({r.allRes.ripplesSurroNormPow}');

% different data groups
groups  = {'Surrogates', [0.7, 0.7, 0.7]; 'Ripples', [0, 0, 0]};
m       = {mean(allRipplesSurroNormPow, 1), mean(allRipplesNormPow, 1)};
ste     = {LK_ste(allRipplesSurroNormPow); LK_ste(allRipplesNormPow)};

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7, 6], 'Color', [1, 1, 1]);
% show the relative power spectra
ax1 = axes('units', 'centimeters', 'position', [1.15, 1.75, 3.5, 3.5]);
hold on;
for iM = 1:numel(m)
    patch([r.param.tfr.foi, fliplr(r.param.tfr.foi)], [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], groups{iM, 2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(r.param.tfr.foi, m{iM}, 'Color', groups{iM, 2});
end
set(gca, 'xscale', 'log', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Frequency (Hz)');
yl = ylabel('Relative power');
tl = title('Relative power');
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
% zoom in on the low frequencies
ax2 = axes('units', 'centimeters', 'position', [2, 3.35, 1.5, 1.5]);
hold on;
for iM = 1:numel(m)
    patch([r.param.tfr.foi, fliplr(r.param.tfr.foi)], [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], groups{iM, 2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(r.param.tfr.foi, m{iM}, 'Color', groups{iM, 2});
end
set(gca, 'xlim', [0, 30], 'xtick', [0, 15, 30], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set(gca, 'fontunits', 'centimeters', 'fontsize', 0.35);
% zoom in on the higher frequencies
ax3 = axes('units', 'centimeters', 'position', [5, 3.35, 1.5, 1.5]);
hold on;
for iM = 1:numel(m)
    patch([r.param.tfr.foi, fliplr(r.param.tfr.foi)], [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], groups{iM, 2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(r.param.tfr.foi, m{iM}, 'Color', groups{iM, 2});
end
set(gca, 'xlim', [80, 140], 'xtick', [80, 110, 140], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set(gca, 'fontunits', 'centimeters', 'fontsize', 0.35);

% save figure
LK_print(f, strcat(paths.save, 'RelativePowerRipplesVSSurrogateRipples_20211022'), '-dtiff', '-r300');

%% figure showing the power spectrum across the entire experiment (gives similar results as the surrogate ripples)

% unfold the power spectrum and estimate mean and standard error
allMeanPow  = cell2mat({r.allRes.meanPow}'); % channels x freqs
m           = mean(allMeanPow, 1);
ste         = LK_ste(allMeanPow);

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5.5, 6], 'Color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [1.65, 1.75, 3.5, 3.5]);
hold on;
patch([r.param.tfr.foi, fliplr(r.param.tfr.foi)], [m + ste, fliplr(m - ste)], [0.7, 0.7, 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(r.param.tfr.foi, m, 'Color', [0.7, 0.7, 0.7]);
set(gca, ...
    'xtick', [1, 10, 100], 'ytick', [10^3, 10^6], ...
    'xscale', 'log', 'yscale', 'log', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Frequency (Hz)');
yl = ylabel('Power (\muV^{2})');
tl = title('Absolute power');
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
% save figure
LK_print(f, strcat(paths.save, 'RawMeanPower_20230318'), '-dtiff', '-r300');
