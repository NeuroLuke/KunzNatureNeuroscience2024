function f = LK_PlotRipplePhaseLockingAcrossChannels_20231030(dt)
%
% LK_PlotRipplePhaseLockingAcrossChannels_20231030 plots the phase locking
% of ripples to delta-frequency activity.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'Color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [1.75, 1.4, 3.5, 3.5]);
hold on;
histogram(dt.allMeanRipplePhaseSurro(:), dt.xEdges, 'normalization', 'probability', ...
    'FaceColor', [1, 1, 1], 'EdgeColor', [0.9, 0.9, 0.9], 'LineWidth', 1);
histogram(dt.allMeanRipplePhase(:), dt.xEdges, 'normalization', 'probability', ...
    'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
t1 = text(0.05, 0.925, 'Empirical', 'units', 'normalized', 'Color', [0.1, 0.1, 0.1]);
t2 = text(0.05, 0.775, 'Surrogate', 'units', 'normalized', 'Color', [0.7, 0.7, 0.7]);
xl = xlabel('Phase (rad)');
yl = ylabel('Probability');
tl = title({'Delta-phase locking', ['(', LK_IndicateThousands(size(dt.allMeanRipplePhase, 1)), ' channels)']});
set(gca, ...
    'xlim', [min(dt.xEdges), max(dt.xEdges)], 'xtick', [-pi, 0, pi], 'xticklabel', {'-\pi', '0', '\pi'}, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off', 'Layer', 'top');
set([gca, t1, t2, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');