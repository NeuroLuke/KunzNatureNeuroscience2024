function f = LK_PlotCoactivityAdditionalInformation_20231015(dt)
%
% LK_PlotCoactivityAdditionalInformation_20231015 plots additional
% information related to the coactivity maps.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7.5, 6]);

% plot data
axes('units', 'centimeters', 'position', [2, 1.8, 4, 4]);
hold on;
imagesc(dt.t, dt.t, dt.data);
% colorbar
colormap parula;
thisCLim = [min(dt.data(:)), max(dt.data(:))];
if all(isnan(thisCLim))
    thisCLim = [-1, 1];
elseif min(thisCLim) == max(thisCLim)
    thisCLim = [min(thisCLim) - 1, max(thisCLim) + 1];
end
caxis(thisCLim);
cb = colorbar('eastoutside', 'units', 'centimeters', 'position', [6.15, 3.05, 0.35, 1.5]);
cb.Ticks = thisCLim;
cb.TickLabels = {num2str(min(thisCLim), '%.0f'), num2str(max(thisCLim), '%.0f')};
title(cb, 'n');
% axis
xl = xlabel({'Time during ripple (s)', 'Place cell'});
yl = ylabel({'Object cell', 'Time during ripple (s)'});
xline(0, '--', 'Color', [0, 0, 0]);
yline(0, '--', 'Color', [0, 0, 0]);
set(gca, 'xtick', [-0.2, 0, 0.2], 'xticklabel', {'-0.2', '0', '0.2'}, 'ytick', [-0.2, 0, 0.2], 'yticklabel', {'-0.2', '0', '0.2'}, 'box', 'on', 'Layer', 'Top');
axis equal;
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
drawnow;
