function f = LK_PlotCoactivity_20231015(dt)
%
% LK_PlotCoactivity_20231015 plots a coactivity map.
%
% Input is a structure with multiple fields.
% Output is a figure handle.
% 
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7.4, 6]);

% plot data for this data group
axes('units', 'centimeters', 'position', [2, 1.8, 4, 4]);
hold on; % low y-values are at the bottom
imagesc(dt.t, dt.t, dt.m); % mean
contour(dt.t, dt.t, dt.sig, 1, '-', 'Color', [1, 1, 1], 'LineWidth', 1); % significant area
colormap parula;
xline(0, '--', 'Color', [0, 0, 0]); % time 0
yline(0, '--', 'Color', [0, 0, 0]); % time 0
set(gca, 'xtick', [-0.2, 0, 0.2], 'xticklabel', {'-0.2', '0', '0.2'}, 'ytick', [-0.2, 0, 0.2], 'yticklabel', {'-0.2', '0', '0.2'}, 'box', 'on');
axis equal;
% indicate p-value
if dt.pValue < 0.001
    t1 = text(0.03, 0.925, '\itP\rm<0.001', 'units', 'normalized', 'Color', [1, 1, 1]);
else
    t1 = text(0.03, 0.925, ['\itP\rm=', num2str(dt.pValue, '%.3f')], 'units', 'normalized', 'Color', [1, 1, 1]);
end
% indicate type of group contrast
tmpAx = get(gca);
plot([min(tmpAx.XLim), max(tmpAx.XLim), max(tmpAx.XLim), min(tmpAx.XLim), min(tmpAx.XLim)], ...
    [min(tmpAx.YLim), min(tmpAx.YLim), max(tmpAx.YLim), max(tmpAx.YLim), min(tmpAx.YLim)], '-', 'linewidth', 1.5, 'Color', dt.color);
% colorbar limits
if contains(dt.condition, 'ret', 'IgnoreCase', true)
    myCLim = dt.coac.clim.ret; % c-limits for retrieval
else
    myCLim = dt.coac.clim.enc; % c-limits for encoding
end
caxis(myCLim);
% axis labeling
xl = xlabel({'Time during ripple (s)', 'Place cell'}); % 2nd dim
yl = ylabel({'Object cell', 'Time during ripple (s)'}); % 1st dim
cb = colorbar('eastoutside', 'units', 'centimeters', 'position', [6.15, 3.05, 0.35, 1.5]);
cb.Ticks = myCLim;
cb.TickLabels = {num2str(min(myCLim), '%.1f'), num2str(max(myCLim), '%.1f')};
cb.LineWidth = 1;
title(cb, dt.coac.type);
set([gca, xl, yl, t1], 'fontunits', 'centimeters', 'fontsize', 0.4);
drawnow;