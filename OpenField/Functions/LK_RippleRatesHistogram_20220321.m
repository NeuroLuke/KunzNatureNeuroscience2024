function f = LK_RippleRatesHistogram_20220321(cfg)
%
% LK_RippleRatesHistogram_20220321 creates a histogram figure showing the
% ripple rates across channels.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.75, 1.4, 3.5, 3.5]);
histogram(cfg.chanRippleRate, cfg.binEdges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
xl = xlabel('Ripple rate (Hz)');
yl = ylabel('Count');
tl = title(cfg.title);
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
set(gca, 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);