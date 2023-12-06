function f = LK_PlotTemporalStabilityOfCueObjCells_20231008(dt)
%
% LK_PlotTemporalStabilityOfCueObjCells_20231008 plots the temporal
% stability of the object cells.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% histogram of temporal stabilities
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.75, 1.75, 3.5, 3.5]);
hold on;
histogram(dt.temporalStab, -1:0.1:1, ...
    'FaceColor', [0.9, 0.9, 0.9], 'FaceAlpha', 1, 'EdgeColor', [0, 0, 0]);
xline(mean(dt.temporalStab, 'omitnan'), 'Color', [1, 0, 0], 'LineWidth', 2);
xl = xlabel('Temporal stability (\itr\rm)');
yl = ylabel('Count');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
set(gca, ...
    'xlim', [-1, 1], 'xtick', -1:0.5:1, ...
    'xticklabelrotation', 0, 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'layer', 'top');