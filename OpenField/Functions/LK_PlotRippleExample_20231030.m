function [f, g] = LK_PlotRippleExample_20231030(dt)
%
% LK_PlotRippleExample_20231030 plots the ERP and the power spectrogram of
% an example ripple.
%
% The input is a structure with multiple fields
%
% The output are two figure handles.
%
% Lukas Kunz, 2023

%% ERP

% create figure
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6], 'Color', [1, 1, 1], 'visible', 'on');

% time-domain data
axes('units', 'centimeters', 'position', [1.6, 3.2, 4, 2.5]);
hold on;
plot(dt.erpData.time{1}, dt.erpData.trial{1}, '-', 'Color', [0, 0, 0]);
logIdx = dt.erpData.time{1} >= (dt.ripple.startTime - dt.ripple.peakTime) & ...
    dt.erpData.time{1} <= (dt.ripple.endTime - dt.ripple.peakTime);
plot(dt.erpData.time{1}(logIdx), dt.erpData.trial{1}(logIdx), '-', 'Color', rgb('darkgreen'));
set(gca, 'xlim', [-0.2, 0.2], 'xtick', -0.2:0.1:0.2, 'xticklabel', '', 'tickdir', 'out');
set(gca, 'fontunits', 'centimeters', 'fontsize', 0.4);

% time-domain filtered data
axes('units', 'centimeters', 'position', [1.6, 1.3, 4, 1.25]);
hold on;
plot(dt.bpData.time{1}, dt.bpData.trial{1}, '-', 'Color', [0, 0, 0]);
plot(dt.bpData.time{1}(logIdx), dt.bpData.trial{1}(logIdx), '-', 'Color', rgb('darkgreen'));
set(gca, 'xlim', [-0.2, 0.2], 'xtick', -0.2:0.1:0.2, 'xticklabelrotation', 0, 'tickdir', 'out');
xl = xlabel('Time (s)');
yl = ylabel('Voltage (\muV)', 'units', 'normalized', 'position', [-0.275, 1.5]);
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);

%% power spectrogram

% create figure
g = figure('units', 'centimeters', 'position', [13, 5, 6, 6], 'Color', [1, 1, 1], 'visible', 'off');

% time-frequency-domain power data
axes('units', 'centimeters', 'position', [1.5, 1.4, 3.5, 3.5]);
hold on;
imagesc(dt.powData.time, dt.powData.freq, dt.normPow);
plot([min(dt.powData.time(logIdx)), max(dt.powData.time(logIdx)), max(dt.powData.time(logIdx)), min(dt.powData.time(logIdx)), min(dt.powData.time(logIdx))], ...
    [min(dt.powData.freq), min(dt.powData.freq), max(dt.powData.freq), max(dt.powData.freq), min(dt.powData.freq)], '--', 'Color', [1, 1, 1], 'LineWidth', 1);
caxis([0, 5]);
cb = colorbar;
cb.Ticks = [0, 5];
cb.Units = 'normalized';
cb.Position = [0.875, 0.35, 0.05, 0.3];
cb.Title.String = 'RP'; % relative power
set(gca, ...
    'xlim', [-0.2, 0.2], 'xtick', -0.2:0.1:0.2, 'xticklabelrotation', 0, ...
    'ylim', [min(dt.powData.freq), max(dt.powData.freq)], 'ytick', [80, 140], ...
    'tickdir', 'out');
xl = xlabel('Time (s)');
yl = ylabel('Frequency (Hz)');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
drawnow;
