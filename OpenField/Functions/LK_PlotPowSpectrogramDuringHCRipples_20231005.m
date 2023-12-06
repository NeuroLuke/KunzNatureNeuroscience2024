function f = LK_PlotPowSpectrogramDuringHCRipples_20231005(dt)
%
% LK_PlotPowSpectrogramDuringHCRipples_20231005 plots power spectrograms
% during hippocampal ripples.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1]);

% depict HC ripples
ax1 = axes('units', 'centimeters', 'position', [2.3, 4, 3.25, 1]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
plot(dt.t, dt.m, 'Color', [0, 0, 0]);
tx = text(0.575, 0.7, 'HC ripple', 'units', 'normalized');
tl = title({dt.LFPROI, ['(', num2str(dt.n), ' combinations)']}, 'units', 'normalized', 'position', [0.5, 1.05]);
set(gca, 'xlim', dt.twoi, 'xticklabel', [], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
set([gca, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
axis off;

% LFP power
ax2 = axes('units', 'centimeters', 'position', [2.3, 0.9, 3.25, 3]);
hold on;
xline(0, '--', 'Color', [0.3, 0.3, 0.3]);
imagesc(dt.timeOI, 1:numel(dt.freq), dt.mPow);
% significances
contour(dt.timeOI, 1:numel(dt.freq), dt.posCBPT.logIdxSigClus, 1, 'k-');
contour(dt.timeOI, 1:numel(dt.freq), dt.negCBPT.logIdxSigClus, 1, 'w-');
xl = xlabel('Time (s)', 'units', 'normalized', 'position', [0.5, -0.09]);
yl = ylabel('Frequency (Hz)', 'units', 'normalized', 'position', [-0.3, 0.5]);
cb = colorbar('southoutside');
cb.Units = 'centimeters';
cb.Position = [0.55, 4.3, 1.2, 0.35];
cb.Limits = [min(dt.cbLim), max(dt.cbLim)];
cb.Ticks = [min(dt.cbLim), max(dt.cbLim)];
cb.Ruler.TickLabelRotation = 0;
xlabel(cb, 'z', 'position', [0, 0]);
caxis([min(dt.cbLim), max(dt.cbLim)]);
set(gca, ...
    'xlim', dt.twoi, 'xtick', [min(dt.twoi), 0, max(dt.twoi)], 'xticklabel', {num2str(min(dt.twoi)), '', num2str(max(dt.twoi))}, ...
    'ylim', [1, numel(dt.freq)], 'ytick', [1, 25, numel(dt.freq)], 'yticklabel', round([dt.freq(1), dt.freq(25), dt.freq(numel(dt.freq))]), ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'on', 'layer', 'top');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'fontsize', 0.4);

% link axes
linkaxes([ax1, ax2], 'x');
