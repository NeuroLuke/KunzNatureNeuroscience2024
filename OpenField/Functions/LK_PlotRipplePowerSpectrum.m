function f = LK_PlotRipplePowerSpectrum(mycfg)
%
% LK_PlotRipplePowerSpectrum plots the average power spectrogram across
% multiple ripples.
%
% The input is a structure with multiple fields.
%
% The ouptut is a figure handle.
%
% Lukas Kunz, 2021

%% ripple TF power

% extract average power across ripples
cfg                 = [];
cfg.trials          = 'all'; % trials to use
cfg.method          = 'wavelet';
cfg.width           = 7;
cfg.output          = 'pow';
cfg.foi             = 2:4:200;
cfg.toi             = 'all';
cfg.keeptrials      = 'no';
cfg.pad             = 'nextpow2';
powRipple           = ft_freqanalysis(cfg, mycfg.ripples.eegRipples);

% baseline correction
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.baseline        = mycfg.param.ripple.twoi4BaselineTF;
cfg.baselinetype    = 'relchange';
powRipple           = ft_freqbaseline(cfg, powRipple);

%% figure for power spectrum

% create figure
f = figure('units', 'centimeters', 'position', mycfg.figEx, 'visible', mycfg.visible);
axes('units', 'centimeters', 'position', [1.5, 1.4, 3.5, 3.5]);
imagesc(powRipple.time, powRipple.freq, squeeze(powRipple.powspctrm));
cb = colorbar('location', 'eastoutside', 'units', 'centimeters', 'position', [5.2, 1.9, 0.25, 2.1]);
cb.Title.String = 'RC';
cb.FontSize = 11;
xl = xlabel('Time (s)');
yl = ylabel('Frequency (Hz)');
tl = title(mycfg.title);
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
set(gca, ...
    'xlim', mycfg.param.ripple.twoi4PlotERP, ...
    'ydir', 'normal', ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
