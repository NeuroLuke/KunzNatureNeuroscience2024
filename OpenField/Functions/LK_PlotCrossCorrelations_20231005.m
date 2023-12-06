function f = LK_PlotCrossCorrelations_20231005(dt)
%
% LK_PlotCrossCorrelations_20231005 plots cross-correlations between
% hippocampal ripples and other ripples.
%
% Lukas Kunz, 2023

% figure showing cross-correlations
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1]);

% HC ripples
axes('units', 'centimeters', 'position', [2.1, 4, 3.25, 1]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
plot(dt.meanRippleTime, dt.meanRipple, 'Color', [0, 0, 0]);
text(0.575, 0.7, 'HC ripple', 'units', 'normalized');
tl = title({dt.LFPROI, ['(\color[rgb]{0, 0, 0.8}', num2str(dt.nSame), '\color[rgb]{0, 0, 0}/\color[rgb]{0.5, 0.5, 0.5}', num2str(dt.nDiff), '\color[rgb]{0, 0, 0} combinations)']}, ...
    'units', 'normalized', 'position', [0.5, 1.05]);
set(gca, 'xlim', dt.twoi, 'xticklabel', [], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
set([gca, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
axis off;

% cross-correlations
axes('units', 'centimeters', 'position', [2.1, 0.9, 3.25, 3]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
yline(0, '--', 'Color', [0.3, 0.3, 0.3]);
% different hemispheres
patch([dt.x, fliplr(dt.x)], [dt.mDiff + dt.steDiff, fliplr(dt.mDiff - dt.steDiff)], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt.x, dt.mDiff, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% same hemispheres
patch([dt.x, fliplr(dt.x)], [dt.mSame + dt.steSame, fliplr(dt.mSame - dt.steSame)], [0, 0, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt.x, dt.mSame, 'Color', [0, 0, 0.8], 'LineWidth', 1);
% enhance axes
xl = xlabel('Time lag (s)', 'units', 'normalized', 'position', [0.5, -0.09]);
yl = ylabel({'X-Correlation', '(unbiased,', 'z-scored)'}, 'units', 'normalized', 'position', [-0.05, 0.5]);
set(gca, 'xlim', dt.twoi, 'xtick', [min(dt.twoi), 0, max(dt.twoi)], 'xticklabel', {num2str(min(dt.twoi)), '', num2str(max(dt.twoi))}, ...
    'ylim', [min(dt.yLim), max(dt.yLim)], 'ytick', [min(dt.yLim), max(dt.yLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
tmpAx = get(gca);
LK_SigLine(dt.x, [max(tmpAx.YLim); max(tmpAx.YLim) - 0.025 * range(tmpAx.YLim)], dt.logIdxSigClus);
drawnow;