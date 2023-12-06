function f = LK_PlotNumObjInsidePF_20231013(dt)
%
% LK_PlotNumObjInsidePF_20231013 plots a figure for showing the number of
% objects within place fields.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5, 6]);
xEdges = -0.5:1:5.5;
yLim = [0, 0.5];

% empirical histogram
axes('units', 'centimeters', 'position', [1.4, 3.75, 3, 1.5]);
histogram(dt.data1, xEdges, 'normalization', 'probability', 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
xline(mean(dt.data1, 'omitnan'), '-', 'Color', [1, 0, 0], 'LineWidth', 2);
t1 = text(0.55, 0.5, 'Empirical', 'units', 'normalized', 'Color', [0.1, 0.1, 0.1]);
tl = title(sprintf('\\itP\\rm = %.3f', dt.pObj));
set([gca, t1, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
set(gca, 'xtick', [], 'ylim', yLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'layer', 'top', 'box', 'off');

% surrogate histogram
axes('units', 'centimeters', 'position', [1.4, 1.75, 3, 1.5]);
histogram(dt.data2, xEdges, 'normalization', 'probability', 'FaceColor', [1, 1, 1], 'EdgeColor', [0.7, 0.7, 0.7]);
t2 = text(0.55, 0.5, 'Surrogate', 'units', 'normalized', 'Color', [0.7, 0.7, 0.7]);
xl = xlabel('# objects in place field');
yl = ylabel('Probability', 'units', 'centimeters', 'position', [-0.9, 1.8], 'HorizontalAlignment', 'center');
set([gca, t2, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'ylim', yLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'layer', 'top', 'box', 'off');