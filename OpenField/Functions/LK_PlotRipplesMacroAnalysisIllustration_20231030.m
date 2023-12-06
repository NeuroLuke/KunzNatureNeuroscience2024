function f = LK_PlotRipplesMacroAnalysisIllustration_20231030(dt)
%
% LK_PlotRipplesMacroAnalysisIllustration_20231030 plots an illustration of
% the ripple detection in macro LFP (iEEG) data.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% spike colors
spikeColors = [0.5, 0.5, 0.5; distinguishable_colors(max(dt.spikesClass(:, 1)))]; % cluster 0 is displayed in gray

% overall figure
f = figure('units', 'centimeters', 'position', [2, 2, 24, 16.5], 'Color', [1, 1, 1]);

%--- macro LFP

% plot raw macro LFP and indicate ripples and artifacts
ax1 = axes('units', 'centimeters', 'position', [3, 13.5, 20, 2.75]);
hold on;
% LFP
bPlot = dt.macroTime >= min(dt.time2Plot) & dt.macroTime <= max(dt.time2Plot);
plot(dt.macroTime(bPlot), dt.macroLFP(bPlot), '-', 'Color', [0.5, 0.5, 0.5]);
% indicate ripples
for iRipple = 1:size(dt.ripples.ripples, 1)
    if dt.ripples.ripples(iRipple).peakTime >= min(dt.time2Plot) && dt.ripples.ripples(iRipple).peakTime <= max(dt.time2Plot)
        startEndIdx = dt.ripples.ripples(iRipple).startIdx:dt.ripples.ripples(iRipple).endIdx;
        plot(dt.macroTime(startEndIdx), dt.macroLFP(startEndIdx), '-', 'Color', rgb('darkgreen')); % ripples
    end
end
% indicate artifacts
for iArt = 1:size(dt.artifacts.artifactOnsEndsSample, 1)
    if (dt.artifacts.artifactOnsEndsTime(iArt, 1) >= min(dt.time2Plot) && dt.artifacts.artifactOnsEndsTime(iArt, 1) <= max(dt.time2Plot)) || ...
            (dt.artifacts.artifactOnsEndsTime(iArt, 2) >= min(dt.time2Plot) && dt.artifacts.artifactOnsEndsTime(iArt, 2) <= max(dt.time2Plot))
        startEndIdx = dt.artifacts.artifactOnsEndsSample(iArt, 1):dt.artifacts.artifactOnsEndsSample(iArt, 2);
        plot(dt.macroTime(startEndIdx), dt.macroLFP(startEndIdx), '-', 'color', rgb('red')); % artifacts
    end
end
tx = text(0.15, 0.4, '0.5-800 Hz', 'units', 'centimeters', 'BackgroundColor', [1, 1, 1, 0.5]); % 4th number in vector is opacity
set(gca, ...
    'xlim', dt.time2Plot, 'xtick', min(dt.time2Plot):max(dt.time2Plot), 'xticklabel', '', ...
    'tickdir', 'out', 'ticklength', [0.005, 0.005]);
yl = ylabel({'Macroelec-', 'trode LFP', '(\muV)'}, ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, tx, yl], 'fontunits', dt.param.plot.fontUnits, 'fontsize', dt.param.plot.fontSize);

%--- macro ripple-band data

% plot ripple band data and indicate ripples and artifacts
ax2 = axes('units', 'centimeters', 'position', [3, 10.5, 20, 2.5]);
hold on;
% ripple-band LFP
plot(dt.macroTime, dt.macroRippleBand, '-', 'Color', [0.5, 0.5, 0.5]);
% indicate ripples
for iRipple = 1:size(dt.ripples.ripples, 1)
    startEndIdx = dt.ripples.ripples(iRipple).startIdx:dt.ripples.ripples(iRipple).endIdx;
    plot(dt.macroTime(startEndIdx), dt.macroRippleBand(startEndIdx), '-', 'Color', rgb('darkgreen')); % ripples
end
% indicate artifacts
for iArt = 1:size(dt.artifacts.artifactOnsEndsSample, 1)
    startEndIdx = dt.artifacts.artifactOnsEndsSample(iArt, 1):dt.artifacts.artifactOnsEndsSample(iArt, 2);
    plot(dt.macroTime(startEndIdx), dt.macroRippleBand(startEndIdx), '-', 'color', rgb('red')); % artifacts
end
tx = text(0.15, 0.4, '80-140 Hz', 'units', 'centimeters', 'BackgroundColor', [1, 1, 1, 0.5]); % 4th number in vector is opacity
set(gca, ...
    'xlim', dt.time2Plot, 'xtick', min(dt.time2Plot):max(dt.time2Plot), 'xticklabel', '', ...
    'tickdir', 'out', 'ticklength', [0.005, 0.005]);
yl = ylabel({'Ripple band', 'LFP', '(\muV)'}, ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, tx, yl], 'fontunits', dt.param.plot.fontUnits, 'fontsize', dt.param.plot.fontSize);

%--- ripple-band envelope with detection threshold

% envelope data with the different thresholds
ax3 = axes('units', 'centimeters', 'position', [3, 7.5, 20, 2.5]);
hold on;
% envelope of ripple-band LFP
plot(dt.macroTime, dt.macroEnvelope, '-', 'Color', [0.5, 0.5, 0.5]);
% indicate ripples
for iRipple = 1:size(dt.ripples.ripples, 1)
    startEndIdx = dt.ripples.ripples(iRipple).startIdx:dt.ripples.ripples(iRipple).endIdx;
    plot(dt.macroTime(startEndIdx), dt.macroEnvelope(startEndIdx), '-', 'Color', rgb('darkgreen')); % ripples
end
% thresholds
plot(dt.macroTime, repmat(dt.threshMin, 1, numel(dt.macroTime)), '--', 'Color', rgb('red')); % minimum threshold
plot(dt.macroTime, repmat(dt.threshPeakMin, 1, numel(dt.macroTime)), '--', 'Color', rgb('orange')); % minimum peak threshold
tx = text(0.15, 0.4, '80-140 Hz', 'units', 'centimeters', 'BackgroundColor', [1, 1, 1, 0.5]); % 4th number in vector is opacity
set(gca, ...
    'xlim', dt.time2Plot, 'xtick', min(dt.time2Plot):max(dt.time2Plot), 'xticklabel', '', ...
    'tickdir', 'out', 'ticklength', [0.005, 0.005]);
yl = ylabel({'Ripple band', 'envelope', '(\muV)'}, ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, tx, yl], 'fontunits', dt.param.plot.fontUnits, 'fontsize', dt.param.plot.fontSize);

%--- power spectrum of the ripple

% plot normalized (i.e., relative) power
ax4 = axes('units', 'centimeters', 'position', [3, 4.5, 20, 2.5]);
imagesc(dt.macroTime, dt.tfData.freq, dt.macroNormPow);
caxis(dt.param.plot.colorbarLim);
cb = colorbar;
cb.Ticks = dt.param.plot.colorbarLim;
cb.Units = 'centimeters';
cb.Position = [22, 4.9, 0.25, 1.7];
cb.Label.String = 'RP';
cb.Label.Position = [-2.5, mean(dt.param.plot.colorbarLim), 0];
cb.Color = [1, 1, 1];
set(gca, ...
    'xlim', dt.time2Plot, 'xtick', min(dt.time2Plot):max(dt.time2Plot), 'xticklabel', '', ...
    'ytick', dt.param.plot.freqTicks, ...
    'tickdir', 'out', 'ticklength', [0.005, 0.005], ...
    'ydir', 'normal', 'box', 'off');
yl = ylabel({'Power', 'spectrogram', 'f (Hz)'}, ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, yl], 'fontunits', dt.param.plot.fontUnits, 'fontsize', dt.param.plot.fontSize);

% plot averaged normalized power across the ripple period
startEndIdx = dt.ripples.ripples(dt.param.plot.rippleIdx).startIdx:dt.ripples.ripples(dt.param.plot.rippleIdx).endIdx;
meanNormPow = mean(dt.macroNormPow(:, startEndIdx), 2, 'omitnan');
axes('units', 'centimeters', 'position', [14.5, 5.15, 1.6, 1.6], 'Color', 'none');
plot(dt.tfData.freq, meanNormPow, '-', 'Color', [1, 1, 1], 'LineWidth', 1);
yl = ylabel('RP', 'units', 'normalized', 'position', [-0.1, 0.5]);
set(gca, ...
    'xlim', [min(dt.tfData.freq), max(dt.tfData.freq)], 'xtick', dt.param.plot.freqTicks, ...
    'ylim', dt.param.plot.colorbarLim, 'ytick', dt.param.plot.colorbarLim, ...
    'LineWidth', 1, 'XColor', [1, 1, 1], 'YColor', [1, 1, 1], 'Color', 'none', 'tickdir', 'out', 'box', 'off');
set([gca, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);

%--- microwire data
ax6 = axes('units', 'centimeters', 'position', [3, 1, 20, 2.5]);
hold on;
% highlight each spike
for iSpike = 1:size(dt.spikesTimeInMacroTime, 1)
    plot([dt.spikesTimeInMacroTime(iSpike), dt.spikesTimeInMacroTime(iSpike)], ...
        [dt.param.plot.microYLim(1), dt.param.plot.microYLim(2)], ...
        'Color', [spikeColors(dt.spikesClass(iSpike, 1) + 1, :), 0.1], 'LineWidth', 4); % 4th input is opacity
end
% filtered microwire data and spike threshold
plot(dt.microTrigTimeInMacroTime, dt.microLFPFiltered, '-', 'Color', [0, 0, 0]);
plot(dt.microTrigTimeInMacroTime, repmat(dt.spikesThresh, 1, numel(dt.microTrigTimeInMacroTime)), '-', 'Color', rgb('red'));
tx = text(0.15, 0.4, '300-3000 Hz', 'units', 'centimeters', 'BackgroundColor', [1, 1, 1, 0.5]); % 4th number in vector is opacity
xl = xlabel('Time (s)', ...
    'units', 'normalized', 'position', [0.5, -0.1]);
yl = ylabel({'Microelec-', 'trode LFP', '(\muV)'}, ...
    'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, tx, xl, yl], 'fontunits', dt.param.plot.fontUnits, 'fontsize', dt.param.plot.fontSize);
set(gca, ...
    'xlim', dt.time2Plot, 'xtick', min(dt.time2Plot):max(dt.time2Plot), ...
    'xticklabel', {min(dt.time2Plot) - min(dt.time2Plot), '', '', '', max(dt.time2Plot) - min(dt.time2Plot)}, ...
    'ylim', [dt.param.plot.microYLim(1), dt.param.plot.microYLim(2)], ...
    'tickdir', 'out', 'ticklength', [0.005, 0.005], 'layer', 'top');

%--- link axes
linkaxes([ax1, ax2, ax3, ax4, ax6], 'x');