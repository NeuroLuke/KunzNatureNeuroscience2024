function [f, percCellsPerReg] = LK_PlotCellsPerRegion(cfg)
%
% LK_PlotCellsPerRegion creates a barplot to present the percentage of
% cells per region.
%
% Input is a structure with fields:
%   uniqueRegions   --> unique regions
%   allRegions      --> region of each cell
%   bCell           --> boolean whether cell belongs to specific cell type
%   ylabel          --> label of y-axis
%   
% Output:
%   f               --> figure handle
%   percPerReg      --> percentage of specific cell type per region
%
% Lukas Kunz, 2022

%% general

% font size
fontSize    = 0.4;
alphaLevel  = 0.05; % alpha level
chancePerc  = 100 * alphaLevel; % chance percentage

% report
fprintf('\nPlotting the percentage of cells per region.\n');

%% figure

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.4, 1.4, 4.3, 4.2]);
hold on;

% loop through regions
percCellsPerReg = nan(numel(cfg.uniqueRegions), 1); % percentage of specific cell type per region
for iReg = 1:numel(cfg.uniqueRegions)
    
    % neurons from this region
    bThisReg                    = strcmp(cfg.allRegions, cfg.uniqueRegions(iReg));
    
    % calculate percentage and p-value for this region
    percCellsPerReg(iReg, 1)    = 100 * sum(cfg.bCell & bThisReg) / sum(bThisReg);
    pBinom                      = myBinomTest(sum(cfg.bCell & bThisReg), sum(bThisReg), alphaLevel);

    % Bonferroni correction of the p-value
    if isfield(cfg, 'BonferroniCorrection') && cfg.BonferroniCorrection == true
        pBinom = pBinom * numel(cfg.uniqueRegions);
        fprintf('Region: %s. Bonferroni-corrected p-value: %.3f.\n', cfg.uniqueRegions{iReg}, pBinom);
    else
        fprintf('Region: %s. Uncorrected p-value: %.3f.\n', cfg.uniqueRegions{iReg}, pBinom);
    end
    
    % plot bar
    b1 = bar(iReg, percCellsPerReg(iReg, 1));
    set(b1, 'FaceColor', [0.9, 0.9, 0.9]);
    % number of cells in this region
    text(b1.XData, 0.25, num2str(sum(bThisReg)), ...
        'FontUnits', 'centimeters', 'fontsize', 0.3, 'Color', [0, 0, 0], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    % indicate significance
    if percCellsPerReg(iReg, 1) > chancePerc
        if pBinom < 0.001
            t = text(b1.XData, b1.YData, '***');
        elseif pBinom < 0.01
            t = text(b1.XData, b1.YData, '**');
        elseif pBinom < 0.05
            t = text(b1.XData, b1.YData, '*');
        end
        % enhance text
        if pBinom < 0.05
            set(t, ...
                'fontunits', 'centimeters', 'fontsize', fontSize, ...
                'verticalalignment', 'baseline', 'horizontalalignment', 'center');
        end
    end
end
hold off;

% indicate chance level
yline(chancePerc, ':', 'Color', rgb('red'), 'LineWidth', 1);

% enhance axes
set(gca, ...
    'xlim', [0.2, numel(cfg.uniqueRegions) + 0.8], ...
    'xtick', 1:numel(cfg.uniqueRegions), 'xticklabel', cfg.uniqueRegions, 'tickdir', 'out', 'XTickLabelRotation', 0);
xl = xlabel('Region');
yl = ylabel(cfg.ylabel);
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', fontSize);
