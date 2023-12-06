function f = LK_PlotNeuronalActivityDuringRipples_20231005(dt)
%
% LK_PlotNeuronalActivityDuringRipples_20231005 plots neuronal firing rates
% during hippocampal ripples.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1]);

% subpanel for ripple
ax1 = axes('units', 'centimeters', 'position', [1.75, 4, 3.75, 1]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
plot(dt.tRipple, dt.mRipple, 'Color', [0, 0, 0]);
text(0.575, 0.7, 'HC ripple', 'units', 'normalized');
tl = title({dt.regionGroup, ['(\color[rgb]{0, 0, 0.8}', num2str(dt.nSame), '\color[rgb]{0, 0, 0}/\color[rgb]{0.5, 0.5, 0.5}', num2str(dt.nDiff), '\color[rgb]{0, 0, 0} combinations)']}, ...
    'units', 'normalized', 'position', [0.5, 1.05]);
set(gca, 'xlim', dt.xLim, 'xticklabel', [], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
set([gca, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
axis off;

% subpanel for single-neuron activity
ax2 = axes('units', 'centimeters', 'position', [1.75, 0.9, 3.75, 3]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
yline(0, '--', 'Color', [0.3, 0.3, 0.3]);
% different hemispheres
patch([dt.t, fliplr(dt.t)], [dt.mDiff + dt.steDiff, fliplr(dt.mDiff - dt.steDiff)], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt.t, dt.mDiff, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% same hemispheres
patch([dt.t, fliplr(dt.t)], [dt.mSame + dt.steSame, fliplr(dt.mSame - dt.steSame)], [0, 0, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt.t, dt.mSame, 'Color', [0, 0, 0.8], 'LineWidth', 1);
% significance
LK_SigLine(dt.t, [max(dt.yLim); max(dt.yLim) - range(dt.yLim) * 0.025], dt.logIdxSigClus);
% enhance axes
xl = xlabel('Time (s)', 'units', 'normalized', 'position', [0.5, -0.09]);
yl = ylabel(dt.yLabel, 'units', 'normalized', 'position', [-0.05, 0.5]);
set(gca, ...
    'xlim', dt.xLim, 'xtick', [min(dt.xLim), 0, max(dt.xLim)], 'xticklabel', {num2str(min(dt.xLim)), '', num2str(max(dt.xLim))}, ...
    'ylim', dt.yLim, 'ytick', dt.yLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'fontsize', 0.4);
drawnow;

% link single-neuron activity and ripples
linkaxes([ax1, ax2], 'x');
