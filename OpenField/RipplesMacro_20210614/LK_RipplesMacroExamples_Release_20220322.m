%==========================================================================
% This script plots example ripples.
%
% Lukas Kunz, 2021
%==========================================================================

% start
clear; close all; clc;
addpath(genpath('E:\OpenField\Functions\'));

% fieldtrip
addpath('E:\fieldtrip\fieldtrip-20210614\');
ft_defaults;

% paths
paths           = [];
paths.ripple    = 'E:\OpenField\RipplesMacro_20210614\20220320_FBS_80to140Hz_tMF2\';
paths.save      = 'E:\OpenField\RipplesMacro_20210614\RippleExamples_20220608\';
mkdir(paths.save);

% settings
param          	= [];
param.subject 	= 'Freiburg_20190130'; % example subject
param.channel   = 'HAR1-HAR2'; % example iEEG channel

%% ripple data

% load ripples from this subject and this channel
ripples         = load(strcat(paths.ripple, param.subject, '\', param.channel, '\ripples.mat'));
ripples         = ripples.ripples;

% loop through possible ripples
for iRipple = 1:size(ripples.ripples, 1)
    
    %% processing of the time domain data
    
    % apply low- and high-pass filters
    cfg             = [];
    cfg.trials      = iRipple;
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 0.5;
    cfg.hpfiltord   = 5;
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 800;
    erpData         = ft_preprocessing(cfg, ripples.eegRipples);
    
    %% filtering of the time domain data in the ripple band
    
    % apply band-pass filter
    cfg             = [];
    cfg.trials      = iRipple;
    cfg.bpfilter    = 'yes';
    cfg.bpfilttype  = 'but';
    cfg.bpfreq      = [80, 140];
    bpData          = ft_preprocessing(cfg, ripples.eegRipples);
    
    %% processing of the time-frequency power data
    
    % apply wavelet power extraction
    cfg             = [];
    cfg.trials      = iRipple;
    cfg.method      = 'wavelet';
    cfg.width       = 7;
    cfg.toi         = 'all';
    cfg.foi         = 30:2:190;
    cfg.pad         = 'nextpow2';
    cfg.output      = 'pow';
    powData         = ft_freqanalysis(cfg, ripples.eegRipples);
    
    % normalize the power data by dividing by the average power
    pow             = squeeze(powData.powspctrm);
    meanPow         = mean(pow, 2, 'omitnan');
    normPow         = pow ./ repmat(meanPow, 1, size(pow, 2));
    
    %% figure for time-domain data and for time-frequency-domain data
    
    % create figure
    dt          = [];
    dt.ripple   = ripples.ripples(iRipple);
    dt.erpData  = erpData;
    dt.bpData   = bpData;
    dt.powData  = powData;
    dt.normPow  = normPow;
    [f, g] = LK_PlotRippleExample_20231030(dt);
    
    % save figure
    set(f, 'InvertHardcopy', 'off');
    LK_print(f, strcat(paths.save, param.subject, '_', param.channel, '_Ripple', num2str(iRipple),'_ERP_20220608'), '-dtiff', '-r300');

    % save figure
    set(g, 'InvertHardcopy', 'off');
    LK_print(g, strcat(paths.save, param.subject, '_', param.channel, '_Ripple', num2str(iRipple),'_Pow_20220608'), '-dtiff', '-r300');

    % save figure data
    if iRipple == 490
        save(strcat(paths.save, 'Fig_2g'), 'dt');
    end

    % close open figures
    close all;
end
