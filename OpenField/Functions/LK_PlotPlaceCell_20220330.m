function f = LK_PlotPlaceCell_20220330(dt)
%
% LK_PlotPlaceCell_20220330 plots firing rate as a function of location.
%
% Input is a structure with multiple fields:
%   visible     --> whether the figure shall be shown
%   place       --> place settings
%   FR          --> firing-rate map
%   PF          --> place field
%   path        --> navigation path
%   thisSpike   --> spike waveforms
%   sr          --> sampling rate
%   nspk        --> number of spikes used for the analysis
%   figTitle    --> figure title
%   idx         --> index (subject, channel, cluster)
%   tPlace      --> t-value from statistical evaluation
%   tPlaceSurro --> surrogate t-values from statistical evaluation
%
% Output is a handle to the figure, f.
%
% Lukas Kunz, 2021

%% general settings

% font size
myFontSize  = 0.4;

% minimum and maximum firing rate
minFR       = min(dt.FR(:));
maxFR    	= max(dt.FR(:));

% set the maximum FR to 1 if it is 0
if maxFR == 0
    minFR   = 0;
    maxFR   = 1;
end

% arena radius
arenaR      = max(dt.place.xEdges);

%% resize the place field for better plotting

% resize the place field
augFac        	= 20; % augmentation factor
PF            	= dt.PF .* 0 + minFR;
PF(dt.PF == 1)  = maxFR;
resPF        	= LK_resizem(PF, augFac); % resize place field

% resize the corresponding x- and y-values for better plotting
resXEdges     	= linspace(min(dt.place.xEdges), max(dt.place.xEdges), size(resPF, 2) + 1);
resXCenters   	= movmean(resXEdges, 2, 'endpoints', 'discard');
resYEdges     	= linspace(min(dt.place.yEdges), max(dt.place.yEdges), size(resPF, 1) + 1);
resYCenters   	= movmean(resYEdges, 2, 'endpoints', 'discard');

%% create figure

% basics
f = figure('units', 'centimeters', 'position', [2, 5, 8, 6], 'Color', [1, 1, 1]);
if isfield(dt, 'visible') && strcmp(dt.visible, 'off')
    set(f, 'visible', 'off');
end
axis off;

% firing-rate map and place field
ax1 = axes('units', 'centimeters', 'position', [2.5, 0.65, 3.75, 3.75]);
hold on;
% firing-rate map
imagesc(dt.place.xCenters, fliplr(dt.place.yCenters), dt.FR, ...
    'AlphaData', ~isnan(dt.FR)); % low y-values are at the bottom
% indicate range of firing rates
cb = colorbar;
cb.Units = 'centimeters';
cb.Position = [6.3, 0.5, 0.2, 1];
cb.YTick = [minFR, maxFR];
cb.YTickLabel = {num2str(minFR, '%.1f'), num2str(maxFR, '%.1f')};
cb.FontSize = 11;
cb.Ruler.TickLabelRotation = 0;
ylabel(cb, 'Hz', 'Rotation', 0, 'VerticalAlignment', 'middle', 'units', 'normalized', 'position', [3, 0.5]);
% colormap
colormap(ax1, 'jet');
caxis([minFR, maxFR]);
% place field
contour(resXCenters, fliplr(resYCenters), resPF, 1, ...
    'Color', [0, 0, 0], 'linewidth', 0.5); % flipping the y-values is necessary so that the correct y-values are associated with their corresponding matrix values
hold off;
% flip y-axis
set(gca, ...
    'ydir', 'reverse', ... % negative y-values are pointing north
    'xlim', arenaR * [-1.05, 1.05], 'ylim', arenaR * [-1.05, 1.05]);
axis equal off;

% boundary
ax2 = axes('units', 'centimeters', 'position', ax1.Position);
hold on;
% plot circle
alpha = -pi:0.001:pi;
patch(arenaR .* [cos(alpha), 1.1 .* fliplr(cos(alpha))], ...
    arenaR .* [sin(alpha), 1.1 .* fliplr(sin(alpha))], [1, 1, 1], 'EdgeColor', 'none');
plot(arenaR .* cos(alpha), arenaR .* sin(alpha), '.', ...
    'Color', [0, 0, 0], 'MarkerSize', 4);
% plot xy-value text
t1 = text(1.15 * arenaR, 0, '(5000/0)', 'Rotation', -90);
t2 = text(0, 1.15 * arenaR, '(0/5000)');
t3 = text(-1.15 * arenaR, 0, '(-5000/0)', 'Rotation', 90);
t4 = text(0, -1.15 * arenaR, '(0/-5000)');
set([t1, t2, t3, t4], ...
    'FontUnits', 'centimeters', 'FontSize', 0.35, 'HorizontalAlignment', 'center');
hold off;
set(gca, ...
    'ydir', 'reverse', ... % negative y-values are pointing north
    'xlim', arenaR * [-1.05, 1.05], 'ylim', arenaR * [-1.05, 1.05]);
axis off;
title(dt.figTitle, ...
    'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.125, 0]);

% link axes of the firing-rate map and the boundary
linkaxes([ax1, ax2]);

% navigation path
axes('units', 'centimeters', 'position', [0.2, 0.25, 2, 2]);
hold on;
plot(dt.path(:, 1), dt.path(:, 2), '-', 'Color', [0.5, 0.5, 0.5]);
if isfield(dt, 'bPathOfInterest')
    plot(dt.path(dt.bPathOfInterest, 1), dt.path(dt.bPathOfInterest, 2), '.', 'Color', [0, 0, 0], 'MarkerSize', 2);
end
plot(cos(alpha) .* arenaR, sin(alpha) .* arenaR, '.', 'Color', [0, 0, 0], 'MarkerSize', 2);
set(gca, ...
    'xlim', arenaR * [-1.025, 1.025], 'ylim', arenaR * [-1.025, 1.025], ...
    'ydir', 'reverse'); % flip y-axis
if numel(dt.idx) == 3
    idxText = [num2str(dt.idx(1)), ' / ', num2str(dt.idx(2)), ' / ', num2str(dt.idx(3))];
elseif numel(dt.idx) == 4
    idxText = [num2str(dt.idx(1)), ' / ', num2str(dt.idx(2)), ' / ', num2str(dt.idx(3)), ' / ', num2str(dt.idx(4))];
end
text(0.5, 1.1, idxText, ...
    'units', 'normalized', 'fontunits', 'centimeters', 'fontsize', 0.2, 'HorizontalAlignment', 'center');
axis equal off;

% statistics
axes('units', 'centimeters', 'position', [6.3, 4, 1.25, 1]);
hold on;
histogram(dt.tPlaceSurro, 10, 'normalization', 'probability', ...
    'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
if ~isnan(dt.tPlace)
    xline(dt.tPlace, 'r-', 'LineWidth', 2);
end
xl = xlabel('\itt', 'Units', 'normalized', 'position', [0.5, -0.175]);
tmpAx = get(gca);
set(gca, ...
    'xlim', [0, max(tmpAx.XLim)], 'xtick', [0, max(tmpAx.XLim)], 'xticklabel', {num2str(0), num2str(round(max(tmpAx.XLim), 1))}, 'xticklabelrotation', 0, ...
    'ytick', [], 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl], 'FontUnits', 'centimeters', 'FontSize', myFontSize);

% spike density plot of the waveforms
axes('units', 'centimeters', 'position', [1, 4, 1.5, 1.5]);
hold on;
outPlot = LK_SpikeDensityPlot(dt.thisSpike, dt.sr);
set(gca, ...
    'xlim', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), 'xtick', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), ...
    'ylim', [outPlot.yMin, outPlot.yMax], 'ytick', [outPlot.yMin, outPlot.yMax], ...
    'ticklength', [0, 0]);
xl = xlabel('ms', 'units', 'normalized', 'position', [0.5, -0.1, 0]);
yl = ylabel('\muV', 'units', 'normalized', 'position', [-0.12, 0.5, 0]);
tx = text(0.5, 1.15, ['n=', num2str(dt.nspk)], ...
    'units', 'normalized', 'HorizontalAlignment', 'center');
set([gca, xl, yl, tx], 'FontUnits', 'centimeters', 'FontSize', myFontSize);

