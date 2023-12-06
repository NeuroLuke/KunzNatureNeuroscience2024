function f = LK_PlotMemoryPerformanceNormTime_20231004(dt)
%
% LK_PlotMemoryPerformanceNormTime_20231004 plots the subjects' memory
% performance as a function of normalized time.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6]);
axes('units', 'centimeters', 'position', [1.625, 1.4, 4, 4]);
hold on;
patch([dt.param.timeCenters, fliplr(dt.param.timeCenters)], [dt.m + dt.ste, fliplr(dt.m - dt.ste)], [0.5, 0.5, 0.5], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt.param.timeCenters, dt.m, 'Color', rgb('black'), 'LineWidth', 1);
set(gca, 'ylim', [0.65, 0.85], 'tickdir', 'out');
xl = xlabel('Normalized time');
yl = ylabel('Memory performance', 'units', 'normalized', 'position', [-0.275, 0.5, 0]);
if dt.pvalMean < 0.001
    tl = title('\itP\rm < 0.001');
else
    tl = title('');
end
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
