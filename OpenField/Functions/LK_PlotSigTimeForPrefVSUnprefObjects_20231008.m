function f = LK_PlotSigTimeForPrefVSUnprefObjects_20231008(dt)
%
% LK_PlotSigTimeForPrefVSUnprefObjects_20231008 plots the significant times
% for preferred versus unpreferred objects, across object cells.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.75, 1.75, 3.5, 3.5]);
hold on;
% mean and SEM
patch([dt.x, fliplr(dt.x)], [dt.m - dt.ste, fliplr(dt.m + dt.ste)], rgb('orange'), 'EdgeColor', 'none', 'facealpha', 0.3);
plot(dt.x, dt.m, 'Color', rgb('orange'));
% cue period
xline(min(dt.param.CBPTTimeEdges), ':', 'Color', [0, 0, 0]);
xline(max(dt.param.CBPTTimeEdges), ':', 'Color', [0, 0, 0]);
xl = xlabel('Time (s)');
yl = ylabel('Fraction of object cells');
set(gca, ...
    'xlim', [min(dt.param.trialTimeEdges), max(dt.param.trialTimeEdges)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off', 'layer', 'top');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);