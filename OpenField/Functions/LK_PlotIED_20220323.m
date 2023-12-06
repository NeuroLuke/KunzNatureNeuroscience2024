function f = LK_PlotIED_20220323(cfg)
%
% LK_PlotIED_20220323 plots interictal epileptic discharges (IEDs).
%
% Input is a structure with multiple fields including:
%   visible         --> whether to show the figure (e.g., 'off');
%   artifactTime    --> time period of the artifact;
%   time            --> time of the data;
%   data            --> data;
%   idx             --> indexing (subject, channel, artifact).
%
% Output is a figure handle.
%
% Lukas Kunz, 2022.

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'visible', cfg.visible);

% create axes
axes('units', 'centimeters', 'position', [1.5, 1.4, 4, 4]);
set(gca, 'xlim', [cfg.artifactTime(1) - 1, cfg.artifactTime(1) + 3] - cfg.artifactTime(1));
hold on;
% plot data
plot(cfg.time - cfg.artifactTime(1), cfg.data, '-', 'Color', [0, 0, 0]);
% plot artifact
logIdx = cfg.time >= cfg.artifactTime(1) & cfg.time <= cfg.artifactTime(2);
plot(cfg.time(logIdx) - cfg.artifactTime(1), cfg.data(logIdx), 'r-', 'MarkerSize', 0.5);
% enhance axes
tmpAx = get(gca);
xl = xlabel('Time (s)');
yl = ylabel('\muV', 'units', 'normalized', 'position', [-0.075, 0.5]);
tl = title(sprintf('(%d / %d / %d)', cfg.idx));
set(gca, 'ylim', [floor(min(tmpAx.YLim)), ceil(max(tmpAx.YLim))], 'ytick', [floor(min(tmpAx.YLim)), ceil(max(tmpAx.YLim))], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
% adjust font size
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
set(tl, 'fontunits', 'centimeters', 'fontsize', 0.2, 'fontweight', 'normal')
drawnow;
