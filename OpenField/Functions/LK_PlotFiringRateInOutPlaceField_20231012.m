function f = LK_PlotFiringRateInOutPlaceField_20231012(dt)
%
% LK_PlotFiringRateInOutPlaceField_20231012 plots the firing rate inside
% vs. outside the place fields (using a boxplot).
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

%% groups

% mean and standard error
m       = {mean(dt.data{1}, 1), mean(dt.data{2}, 1)};
ste     = {LK_ste(dt.data{1}), LK_ste(dt.data{2})};

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5, 6.5]);

% bar plot for means
axes('units', 'centimeters', 'position', [1.2, 1.4, 3, 2.5]);
hold on;
for iG = 1:size(dt.groups, 1)
    % mean +/- standard error
    bar(iG, m{iG}, 'FaceColor', dt.groups{iG, 2}, 'FaceAlpha', 0.5);
    plot([iG, iG], [m{iG} - ste{iG}, m{iG} + ste{iG}], '-', 'Color', [0, 0, 0], 'LineWidth', 2);
end
xl = xlabel('place field');
yl = ylabel('Firing rate (Hz)');
set(gca, ...
    'xlim', [0.4, size(dt.groups, 1) + 0.6], 'xtick', 1:size(dt.groups, 1), ...
    'xticklabel', dt.groups(:, 1), 'xticklabelrotation', 0, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);

% histogram of differences
axes('units', 'centimeters', 'position', [2.2, 4.75, 2.4, 1.5]);
xEdges = 0:0.2:LK_ceil(max(dt.data{1} - dt.data{2}), 1);
histogram(dt.data{1} - dt.data{2}, xEdges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
xline(mean(dt.data{1} - dt.data{2}), 'color', [1, 0, 0], 'LineWidth', 2);
xl = xlabel('\DeltaFR (Hz)', 'units', 'normalized', 'position', [0.45, -0.15]);
yl = ylabel({'# place', 'cells'}, 'units', 'normalized', 'position', [-0.1, 0.5]);
tx = text(0.45, 0.65, {'Inside-', 'Outside'}, 'units', 'normalized', 'Color', [0.5, 0.5, 0.5]);
tmpAx = get(gca);
set(gca, ...
    'xlim', [min(xEdges), max(xEdges)], 'xtick', [min(xEdges), max(xEdges)], ...
    'ylim', [min(tmpAx.YLim), max(tmpAx.YLim)], 'ytick', [min(tmpAx.YLim), max(tmpAx.YLim)], 'tickdir', 'out', ...
    'box', 'off')
set([gca, xl, yl, tx], 'FontUnits', 'centimeters', 'FontSize', 0.4);
