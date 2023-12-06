%==========================================================================
% This script performs simulations to examine the relationship between the
% firing rates of the neurons and their coactivity.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('E:\OpenField\Functions\'));
rng(1); % for reproducibility

% paths
paths       = [];
paths.save  = 'E:\OpenField\NeuronalCoactivityDuringRipples_20211025\CoactivityIllustration\';

% settings
numRipples  = 100;
FRs         = transpose(0:numRipples);
numTests   	= 1000; % number of tests per firing rate (e.g., 1000)

%% coactivity

% preallocate
allFR1      = nan(size(FRs, 1), size(FRs, 1), numTests);
allFR2      = nan(size(FRs, 1), size(FRs, 1), numTests);
allZ        = nan(size(FRs, 1), size(FRs, 1), numTests);

% loop through possible firing rates
for iFR1 = 1:size(FRs, 1)
    for iFR2 = 1:size(FRs, 1)
        
        % firing rates in this round
        FR1 = FRs(iFR1, 1);
        FR2 = FRs(iFR2, 1);
        fprintf('Firing rates: %d and %d.\n', FR1, FR2);
        
        % loop through tests
        for iTest = 1:numTests
            
            % activity cell 1
            cell1       = zeros(numRipples, 1);
            i1          = datasample(1:size(cell1, 1), FR1, 'replace', false); % ripples in which the cell is active
            cell1(i1)   = 1;
            
            % activity cell 2
            cell2       = zeros(numRipples, 1);
            i2          = datasample(1:size(cell2, 1), FR2, 'replace', false); % ripples in which the cell is active
            cell2(i2)   = 1;
            
            % coactivity z-score
            z = LK_CoactivityZScore_20220615(cell1, cell2);
            
            % collect all coactivity z-scores
            allFR1(iFR1, iFR2, iTest)   = FR1;
            allFR2(iFR1, iFR2, iTest)   = FR2;
            allZ(iFR1, iFR2, iTest)     = z;
        end
    end
end

%% figure: z-score as a function of the firing rates of cell 1 and 2

% x-values, y-values, and height values
x = allFR1;
y = allFR2;
h = allZ;

% plot the coactivity values for different firing-rate groups
groups  = {'FR1'; 'FR2'; 'FR1andFR2'};
for iG = 1:size(groups, 1)
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 8.5, 8]);
    axes('units', 'centimeters', 'position', [1.5, 1.5, 6, 6]);
    scatter3(x(:), y(:), h(:), 3, h(:), 'filled');
    if strcmp(groups{iG}, 'FR1')
        view([0, 0]);
    elseif strcmp(groups{iG}, 'FR2')
        view([90, 0]);
    elseif strcmp(groups{iG}, 'FR1andFR2')
        view([45, 20]);
    end
    grid minor;
    xl = xlabel('# ripples cell B active');
    yl = ylabel('# ripples cell A active');
    zl = zlabel('Coactivity (z)');
    if strcmp(groups{iG}, 'FR1andFR2')
        set(xl, 'units', 'normalized', 'position', [0.1, -0.1]);
        set(yl, 'units', 'normalized', 'position', [0.825, -0.1]);
        set(zl, 'units', 'normalized', 'position', [-0.125, 0.5]);
    end
    colormap parula;
    set(gca, 'xlim', [min(FRs), max(FRs)], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
    set([gca, xl, yl, zl], 'fontunits', 'centimeters', 'fontsize', 0.4);
    % save figure
    LK_print(f, strcat(paths.save, groups{iG, 1}, 'vsZ_20230330'), '-dpng', '-r300');
end

%% t-test across tests

% t-test against 0 for each combination of firing rates
[~, p, ~, stats]   = ttest(allZ, 0, 'dim', 3);
fprintf('How many of the t-tests are significant? %d (%.3f%%)\n', sum(p(:) < 0.05), 100 * sum(p(:) < 0.05) / numel(p));

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 13, 6]);

% show all t-values as a function of the firing rates
axes('units', 'centimeters', 'position', [1.5, 1.5, 4, 4]);
hold on;
imagesc(FRs, FRs, stats.tstat);
xl = xlabel('# ripples cell B active');
yl = ylabel('# ripples cell A active');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'xlim', [min(FRs), max(FRs)], 'ylim', [min(FRs), max(FRs)]);
myCLim = [-4, 4];
caxis(myCLim);
cb = colorbar('eastoutside', 'units', 'centimeters', 'position', [5.75, 2.5, 0.25, 2]);
cb.Title.String = '\itt';
cb.Ticks = myCLim;
cb.TickLabels = {sprintf('%.1f', myCLim(1)), sprintf('%.1f', myCLim(2))};

% show histogram of the t-values
axes('units', 'centimeters', 'position', [8.5, 1.5, 4, 4]);
histogram(stats.tstat(:), linspace(myCLim(1), myCLim(2), 21), 'FaceColor', 'none', 'EdgeColor', [0, 0 0]);
xl = xlabel('\itt');
yl = ylabel('Count');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'xlim', myCLim, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');

% save figure
LK_print(f, strcat(paths.save, 'FR2tStatAcrossTests_20230330'), '-dpng', '-r300');
