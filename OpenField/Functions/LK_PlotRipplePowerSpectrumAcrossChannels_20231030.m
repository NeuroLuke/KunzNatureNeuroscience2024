function f = LK_PlotRipplePowerSpectrumAcrossChannels_20231030(dt)
%
% LK_PlotRipplePowerSpectrumAcrossChannels_20231030 plots the power
% spectrogram of ripples across channels.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', dt.groups{dt.iGroup, 2}, 'visible', 'on');
axes('units', 'centimeters', 'position', dt.groups{dt.iGroup, 3});

% plot power spectrogram
imagesc(dt.data.time, dt.param.ripple.freq4PlotTF, dt.data.avg);

% enhance axes
set(gca, 'xlim', dt.groups{dt.iGroup, 4}, 'ydir', 'normal');
xl = xlabel('Time (s)');
yl = ylabel('Frequency (Hz)');
tl = title({'Ripple power', ['(', num2str(dt.numChannels), ' channels)']});
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
cb = colorbar;
colormap parula;
cb.Units = 'centimeters';
cb.Position = [5.25, 2.1, 0.3, 1.8];
cb.Ticks = dt.groups{dt.iGroup, 5};
cb.Title.String = 'RC';
caxis(dt.groups{dt.iGroup, 5});
