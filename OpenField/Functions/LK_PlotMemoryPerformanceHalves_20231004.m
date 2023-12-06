function f = LK_PlotMemoryPerformanceHalves_20231004(dt)
%
% LK_PlotMemoryPerformanceHalves_20231004 plots the subjects' memory
% performance as a function of trial half.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 5, 6]);
axes('units', 'centimeters', 'position', [1.5, 1.4, 3, 4]);
hold on;
for iSub = 1:size(dt.subjects, 1)
    if dt.bMaskSession(iSub) == 1
        p = plot([1, 2], dt.allMemPerfEarlyLate(iSub, :), '-', 'Color', rgb('gray'));
        if strcmp(dt.subjects{iSub, 2}, 'Microwire')
            set(p, 'Color', [0, 0, 0]); % indicate microwire subjects in black
        end
    end
end
plot([1, 2], mean(dt.allMemPerfEarlyLate(dt.bMaskSession, :)), '-', 'Color', rgb('blue'), 'LineWidth', 4); % average across sessions
if dt.pHalves < 0.001
    tl = title('\itP\rm < 0.001');
else
    tl = title('');
end
hold off;
% enhance axes
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'Early', 'Late'}, 'ylim', [0.4, 1.05], 'ytick', [0.4, 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('trials', 'units', 'normalized', 'position', [0.5, -0.21]);
yl = ylabel('Memory performance', 'units', 'normalized', 'position', [-0.275, 0.5, 0]);
set([gca, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
set(gcf, 'PaperPositionMode', 'auto');
