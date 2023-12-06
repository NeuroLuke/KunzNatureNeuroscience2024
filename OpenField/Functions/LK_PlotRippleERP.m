function f = LK_PlotRippleERP(mycfg)
%
% LK_PlotRippleERP plots the ERP across multiple ripples.
%
% The input is a structure with multiple fields.
%
% The ouptut is a figure handle.
%
% Lukas Kunz, 2021

%% baseline correction
cfg             = [];
cfg.baseline    = mycfg.param.ripple.twoi4BaselineERP;
erpRipple       = ft_timelockbaseline(cfg, mycfg.ripples.eegRipples);

%% ERP
cfg             = [];
erpRipple       = ft_timelockanalysis(cfg, erpRipple);

%% figure for ERP
f = figure('units', 'centimeters', 'position', mycfg.figEx, 'visible', mycfg.visible);
axes('units', 'centimeters', 'position', [1.7, 1.4, 3.5, 3.5]);
plot(erpRipple.time, erpRipple.avg, 'k-');
xl = xlabel('Time (s)');
yl = ylabel('Voltage (\muV)');
tl = title(mycfg.title);
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
set(gca, ...
    'xlim', mycfg.param.ripple.twoi4PlotERP, ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);