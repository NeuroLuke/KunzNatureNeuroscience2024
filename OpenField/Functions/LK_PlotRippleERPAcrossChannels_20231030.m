function f = LK_PlotRippleERPAcrossChannels_20231030(dt)
%
% LK_PlotRippleERPAcrossChannels_20231030 plots the ripple ERP across
% multiple channels.
%
% The input is a structure with multiple fields.
%
% The ouptut is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', dt.groups{dt.iGroup, 2}, 'visible', 'on');
axes('units', 'centimeters', 'position', dt.groups{dt.iGroup, 3});
hold on;

% mean and SEM
if ~contains(dt.groups{dt.iGroup, 1}, 'MeanOnly')
    patch([dt.time, fliplr(dt.time)], [dt.m - dt.ste, fliplr(dt.m + dt.ste)], [0.7, 0.7, 0.7], 'EdgeColor', 'none'); % SEM
end
plot(dt.time, dt.m, 'Color', [0, 0, 0]); % mean

% enhance axes
set(gca, ...
    'xlim', dt.groups{dt.iGroup, 4}, 'ylim', dt.groups{dt.iGroup, 5}, ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'xticklabelrotation', 0);
if ~isempty(dt.groups{dt.iGroup, 6})
    % indicate zoom-in window
    plot([min(dt.groups{dt.iGroup, 6}), max(dt.groups{dt.iGroup, 6}), max(dt.groups{dt.iGroup, 6}), min(dt.groups{dt.iGroup, 6}), min(dt.groups{dt.iGroup, 6})], ...
        [min(dt.groups{dt.iGroup, 5}), min(dt.groups{dt.iGroup, 5}), max(dt.groups{dt.iGroup, 5}), max(dt.groups{dt.iGroup, 5}), min(dt.groups{dt.iGroup, 5})], 'k--');
    % labels
    xl = xlabel('Time (s)');
    yl = ylabel('Voltage (\muV)');
    tl = title({'Grand-average ripple', ['(', num2str(dt.numChannels), ' channels)']});
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', dt.groups{dt.iGroup, 7}, 'fontweight', 'normal');
else
    % adjust ticks
    set(gca, 'xtick', dt.groups{dt.iGroup, 4}, 'ytick', dt.groups{dt.iGroup, 5});
    set(gca, 'fontunits', 'centimeters', 'fontsize', dt.groups{dt.iGroup, 7}, 'fontweight', 'normal');
end
if contains(dt.groups{dt.iGroup, 1}, 'MeanOnly')
    axis off;
end
