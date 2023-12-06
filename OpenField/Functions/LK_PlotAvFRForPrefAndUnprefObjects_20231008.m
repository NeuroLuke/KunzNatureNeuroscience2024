function f = LK_PlotAvFRForPrefAndUnprefObjects_20231008(dt)
%
% LK_PlotAvFRForPrefAndUnprefObjects_20231008 plots the average firing rate
% for preferred and unpreferred objects.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.75, 1.75, 3.5, 3.5]);
hold on;
for iG = 1:size(dt.groups, 1)
    % standard error
    patch([dt.x, fliplr(dt.x)], [dt.m{iG} - dt.ste{iG}, fliplr(dt.m{iG} + dt.ste{iG})], dt.groups{iG, 2}, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    % mean
    plot(dt.x, dt.m{iG}, 'Color', dt.groups{iG, 2});
end
% enhance axes
LK_yline(0, '--', [0.1, 0.1, 0.1], 'bottom');
% cue period
LK_xline(min(dt.param.CBPTTimeEdges), ':', [0, 0, 0], 'bottom');
LK_xline(max(dt.param.CBPTTimeEdges), ':', [0, 0, 0], 'bottom');
xl = xlabel('Time (s)');
yl = ylabel('Relative FR (Hz)', 'units', 'normalized', 'position', [-0.2, 0.5, 0]);
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
set(gca, ...
    'xlim', [min(dt.param.trialTimeEdges), max(dt.param.trialTimeEdges)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);