function f = LK_PlotTrialPeriodDurations_20231005(dt)
%
% LK_PlotTrialPeriodDurations_20231005 plots the durations of the trial
% periods.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023.

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5, 5]);

% histogram for durations
axes('units', 'centimeters', 'position', [1.85, 2, 2.6, 2.6]);
hold on;
histogram(dt.dataAll, min(dt.myXLim):2:max(dt.myXLim), 'normalization', 'probability', ...
    'edgecolor', dt.condition, 'facecolor', dt.condition, 'facealpha', 0.2, 'linewidth', 1);
xline(mean(dt.dataAll, 'omitnan'), '-', 'Color', [0, 0, 0], 'LineWidth', 2);
xl = xlabel('Duration (s)');
yl = ylabel('% trials', 'units', 'normalized', 'position', [-0.4, 0.5]);
set(gca, 'xlim', dt.myXLim, 'ylim', dt.myYLim, 'YTick', dt.myYLim', 'YTickLabel', dt.myYLim .* 100, ...
    'tickdir', 'out', 'ticklength', [0.04, 0.04], 'Layer', 'Top');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.6);