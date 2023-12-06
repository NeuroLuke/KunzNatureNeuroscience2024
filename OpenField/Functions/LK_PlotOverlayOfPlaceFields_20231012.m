function f = LK_PlotOverlayOfPlaceFields_20231012(dt)
%
% LK_PlotOverlayOfPlaceFields_20231012 plots the overlay of all place
% fields.
%
% Lukas Kunz, 2023

% arena radius for plotting the boundary
arenaR  = 5000;

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'Color', [1, 1, 1]);

% distribution of all place fields
ax1 = axes('units', 'centimeters', 'position', [0.8, 1.25, 3.75, 3.75]);
hold on;
% place-field distribution
imagesc(dt.place.xCenters, fliplr(dt.place.yCenters), dt.sumPF, ...
    'AlphaData', ~isnan(dt.sumPF)); % low y-values are at the bottom
% colorbar
caxis([min(dt.sumPF(:)), LK_ceil(max(dt.sumPF(:)), 0)]);
cb = colorbar('location', 'southoutside', 'YTick', [min(dt.sumPF(:)), LK_ceil(max(dt.sumPF(:)), 0)], ...
    'Units', 'centimeters', 'Position', [4.25, 4.5, 1.25, 0.25], 'FontSize', 11);
title(cb, '%', 'units', 'normalized', 'position', [0.5, 1.5]);
cb.Ruler.TickLabelRotation = 0;
% colormap
colormap(ax1, parula);
hold off;
% flip y-axis
set(gca, ...
    'ydir', 'reverse', ... % low y-values are pointing north due to the way the arena is programmed in the task
    'xlim', arenaR * [-1.05, 1.05], 'ylim', arenaR * [-1.05, 1.05]);
axis equal off;

% boundary
ax2 = axes('units', 'centimeters', 'position', ax1.Position);
hold on;
% plot circle
alpha = -pi:0.001:pi;
patch(arenaR .* [cos(alpha), 1.1 .* fliplr(cos(alpha))], ...
    arenaR .* [sin(alpha), 1.1 .* fliplr(sin(alpha))], [1, 1, 1], 'EdgeColor', 'none');
plot(arenaR .* cos(alpha), arenaR .* sin(alpha), 'Color', [0, 0, 0], 'LineWidth', 2);
% coordinates
t1 = text(1.175 * arenaR, 0, '(5000/0)', 'Rotation', -90);
t2 = text(0, 1.175 * arenaR, '(0/5000)');
t3 = text(-1.175 * arenaR, 0, '(-5000/0)', 'Rotation', 90);
t4 = text(0, -1.175 * arenaR, '(0/-5000)');
set([t1, t2, t3, t4], ...
    'FontUnits', 'centimeters', 'FontSize', 0.4, 'HorizontalAlignment', 'center');
hold off;
set(gca, ...
    'ydir', 'reverse', ... % negative y-values are pointing north
    'xlim', arenaR * [-1.05, 1.05], 'ylim', arenaR * [-1.05, 1.05]);
axis off;

% link axes of the firing-rate map and the boundary
linkaxes([ax1, ax2]);

% adjust figure layout
set(f, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
