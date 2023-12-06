%==========================================================================
% This script performs simulations to explain the conversion between drop
% errors and memory-performance values.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));

% paths
paths           = [];
paths.save      = 'E:\OpenField\MemoryPerformance_20230324\';

% settings
param           = [];
param.arenaR    = 5000;
param.arenaCtr  = [0, 0];
param.corrLocs  = [0, 0; 2500, 0; 5000, 0]; % example correct locations
param.respLocs  = param.corrLocs + [-3000, 0]; % example response locations
param.cLim      = [0, 2 * param.arenaR];

% random locations in the environment to normalize drop errors
cfg             = [];
cfg.maxR        = param.arenaR;
cfg.minR        = 0;
cfg.N           = 1000001;
cfg.centerX     = param.arenaCtr(1);
cfg.centerY     = param.arenaCtr(2);
randLocs        = LK_RandomPointsInCircle(cfg);

%% loop through example correct locations
for iL = 1:size(param.corrLocs, 1)
    
    % drop errors
    dropError       = pdist2(param.respLocs(iL, :), param.corrLocs(iL, :));
    dropErrorSurro  = transpose(pdist2(param.corrLocs(iL, :), randLocs));
    memPerf         = sum(dropError < dropErrorSurro) / numel(dropErrorSurro);

    %% create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 12, 6], 'color', [1, 1, 1]);
    
    % plot the surrogate drop error across the arena
    axes('units', 'centimeters', 'position', [0.5, 1.75, 4, 4]);
    hold on;
    scatter(randLocs(:, 1), randLocs(:, 2), 5, dropErrorSurro, 'filled');
    plot(cosd(0:0.01:360) .* param.arenaR, sind(0:0.01:360) .* param.arenaR, '.', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10);
    plot(param.corrLocs(iL, 1), param.corrLocs(iL, 2), 'x', 'Color', [0, 0, 0], 'MarkerSize', 10, 'LineWidth', 2);
    plot(param.respLocs(iL, 1), param.respLocs(iL, 2), 'x', 'Color', [1, 0, 0], 'MarkerSize', 10, 'LineWidth', 2);
    plot([param.corrLocs(iL, 1), param.respLocs(iL, 1)], [param.corrLocs(iL, 2), param.respLocs(iL, 2)], '--', 'Color', [1, 0, 0], 'LineWidth', 1);
    cm = colormap(flipud(winter));
    cb = colorbar('SouthOutside', 'units', 'centimeters', 'position', [0.75, 0.75, 3.5, 0.6], 'FontSize', 11);
    ylabel(cb, 'Drop error (vu)', 'units', 'normalized', 'position', [0.5, 1], 'fontunits', 'centimeters', 'fontsize', 0.4, 'color', [1, 1, 1]);
    caxis(param.cLim);
    axis off;
    
    % plot a histogram of all surrogate drop errors and their relation to
    % the empirical drop error
    axes('units', 'centimeters', 'position', [7, 1.4, 4, 4]);
    [N, edges] = histcounts(dropErrorSurro, linspace(min(param.cLim), max(param.cLim), size(cm, 1) + 1));
    binWidth = edges(2) - edges(1);
    centers = movmean(edges, 2, 'endpoints', 'discard');
    hold on;
    for iN = 1:numel(N)
        bar(centers(iN), N(iN), 'FaceColor', cm(iN, :), 'EdgeColor', 'none', 'BarWidth', binWidth);
    end
    xline(dropError, '--', 'LineWidth', 1, 'Color', [1, 0, 0]);
    xl = xlabel('Drop error (vu)');
    yl = ylabel('Frequency');
    tl = title(sprintf('Memory performance = %.3f', memPerf));
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    set(gca, 'ytick', [], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'Layer', 'Top');
    set(gcf, 'InvertHardcopy', 'off');
    LK_print(f, strcat(paths.save, 'DropError2MemPerf_X', num2str(param.corrLocs(iL, 1)), 'Y', num2str(param.corrLocs(iL, 2))), '-dpng', '-r300');
end