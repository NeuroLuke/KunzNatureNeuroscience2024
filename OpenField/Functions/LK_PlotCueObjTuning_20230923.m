function out = LK_PlotCueObjTuning_20230923(dt)
%
% LK_PlotCueObjTuning_20230923 plots various information about object
% tuning during the cue period.
%
% Input is a structure with multiple fields for plotting.
%
% Output is a structure with fields:
%   f   --> figure handle
%   m   --> average firing rates for the unpreferred and the preferred
%               objects
%   ste --> standard error of the mean for the firing rates of the
%               unpreferred and the preferred objects
%
% Lukas Kunz, 2022

%% figure

% font size
myFontSize  = 0.4;

% create figure
f = figure('units', 'centimeters', 'position', [2, 5, 16, 6], 'Color', [1, 1, 1], 'visible', dt.visible);

% average firing rate per object
axes('units', 'centimeters', 'position', [4, 1.4, 3.5, 3.5]);
hold on;
for iObj = 1:size(dt.prefObjTest.objFR, 1)
    % minimum and maximum
    thisMin = dt.prefObjTest.objFRMin(iObj, 1);
    thisMax = dt.prefObjTest.objFRMax(iObj, 1);
    plot(iObj + [-0.2, 0.2, 0, 0, -0.2, 0.2], [thisMin, thisMin, thisMin, thisMax, thisMax, thisMax], '-', 'Color', [0, 0, 0]);
    % 25th and 75th percentile
    thisP25 = dt.prefObjTest.objFRp25(iObj, 1);
    thisP75 = dt.prefObjTest.objFRp75(iObj, 1);
    thisPatch = patch(iObj + [-0.3, 0.3, 0.3, -0.3, -0.3], [thisP25, thisP25, thisP75, thisP75, thisP25], [1, 1, 1], 'EdgeColor', [0, 0, 0], 'LineWidth', 1);
    % preferred object
    if iObj == dt.prefObjTest.prefObjIdx
        set(thisPatch, 'FaceColor', rgb('orange'));
    end
    % median
    thisMed = dt.prefObjTest.objFRMed(iObj, 1);
    plot(iObj + [-0.3, 0.3], [thisMed, thisMed], '-', 'LineWidth', 1, 'Color', [0, 0, 0]);
end
% enhance the axes
minY = 0;
maxY = ceil(max(dt.prefObjTest.objFRMax) + 0.1);
set(gca, ...
    'xlim', [0, numel(dt.prefObjTest.objFR) + 1], 'xtick', 1:numel(dt.prefObjTest.objFR), ...
    'ylim', [minY, maxY], 'ytick', [minY, maxY], ...
    'fontunits', 'centimeters', 'fontsize', myFontSize, ...
    'xticklabelrotation', 0, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
xl = xlabel('Object');
yl = ylabel('FR (Hz)', 'units', 'normalized', 'position', [-0.04, 0.5, 0]);
if dt.objTP < 0.001
    figTitle = {'Object tuning', ['\itP\rm < 0.001 (', dt.chanRegion, ')']};
else
    figTitle = {'Object tuning', ['\itP\rm = ', num2str(dt.objTP, '%.3f'), ' (', dt.chanRegion, ')']};
end
tl = title(figTitle, 'Units', 'normalized', 'Position', [0.5, 1, 0]);
set([gca, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal');

% averaged time-resolved firing rates
axes('units', 'centimeters', 'position', [8.65, 1.4, 3.5, 3.5]);
hold on;
% time, mean, and standard error
t = dt.param.trialTimeCenters;
m = {mean(dt.trFR(dt.trialInfo.Object ~= dt.prefObjTest.prefObjName, :)); ...
    mean(dt.trFR(dt.trialInfo.Object == dt.prefObjTest.prefObjName, :))}; % unpreferred, preferred
ste = {LK_ste(dt.trFR(dt.trialInfo.Object ~= dt.prefObjTest.prefObjName, :)); ...
    LK_ste(dt.trFR(dt.trialInfo.Object == dt.prefObjTest.prefObjName, :))}; % unpreferred, preferred
patch([t, fliplr(t)], [m{1} - ste{1}, fliplr(m{1} + ste{1})], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % unpreferred
patch([t, fliplr(t)], [m{2} - ste{2}, fliplr(m{2} + ste{2})], rgb('orange'), 'EdgeColor', 'none', 'FaceAlpha', 0.5); % preferred
plot(t, m{1}, 'Color', [0.5, 0.5, 0.5]);
plot(t, m{2}, 'Color', rgb('orange'));
% indicate cue period
LK_xline(min(dt.param.CBPTTimeEdges), ':', [0, 0, 0], 'bottom');
LK_xline(max(dt.param.CBPTTimeEdges), ':', [0, 0, 0], 'bottom');
% enhance axes
LK_yline(0, '--', [0.1, 0.1, 0.1], 'bottom');
xl = xlabel('Time (s)');
yl = ylabel({'Relative FR', '(Hz)'}, 'units', 'normalized', 'position', [-0.04, 0.5, 0]);
tmpAx = get(gca);
set(gca, ...
    'xtick', dt.param.CBPTTimeEdges, ...
    'ylim', [LK_floor(min(tmpAx.YLim), 1), LK_ceil(max(tmpAx.YLim), 1)], 'ytick', [LK_floor(min(tmpAx.YLim), 1), LK_ceil(max(tmpAx.YLim), 1)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02]);
% show significance
if ~isempty(dt.objCBPTSigClusTime)
    sigTime = max(dt.objCBPTSigClusTime) >= dt.param.trialTimeCenters & min(dt.objCBPTSigClusTime) <= dt.param.trialTimeCenters;
    LK_SigLine(dt.param.trialTimeCenters, [min(tmpAx.YLim); min(tmpAx.YLim) + range(tmpAx.YLim) * 0.025], sigTime);
end
% title
figTitle = {'Time-resolved object tuning', ['\itP\rm = ', num2str(dt.objCBPTSigP, '%.3f')]};
tl = title(figTitle, 'Units', 'normalized', 'Position', [0.5, 1, 0]);
set([gca, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', myFontSize, 'FontWeight', 'normal');

% raster plot for time-resolved firing rates
ax3 = axes('units', 'centimeters', 'position', [12.3, 1.4, 3.5, 3.5]);
hold on;
% raster data
rasterData = [dt.trNumSpikes(dt.trialInfo.Object == dt.prefObjTest.prefObjName, :) > 0; ...
    (dt.trNumSpikes(dt.trialInfo.Object ~= dt.prefObjTest.prefObjName, :) > 0) .* 0.5]; % preferred, unpreferred
imagesc(dt.param.trialTimeCenters, 1:size(rasterData, 2), rasterData);
colormap(ax3, [1, 1, 1; [0.5, 0.5, 0.5]; rgb('orange')]);
% indicate cue period
LK_xline(min(dt.param.CBPTTimeEdges), '-', [0, 0, 0], 'top');
LK_xline(max(dt.param.CBPTTimeEdges), '-', [0, 0, 0], 'top');
% enhance axes
xl = xlabel('Time (s)');
yl = ylabel('Trials', 'units', 'normalized', 'position', [0.14, 0.5, 0]);
set(gca, ...
    'xtick', dt.param.CBPTTimeEdges, 'ytick', [], 'ydir', 'reverse', ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off', 'layer', 'top');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', myFontSize);
title([num2str(dt.idx(1)), ' / ', num2str(dt.idx(2)), ' / ', num2str(dt.idx(3))], ...
    'fontunits', 'centimeters', 'fontsize', 0.2, 'fontweight', 'normal');

% object locations
axes('units', 'centimeters', 'position', [0.9, 0.6, 1.5, 1.5]);
hold on;
for iObj = min(dt.trialInfo.Object):max(dt.trialInfo.Object)
    % this object location
    thisObjXY   = unique([dt.trialInfo.xCorrect(dt.trialInfo.Object == iObj), dt.trialInfo.yCorrect(dt.trialInfo.Object == iObj)], 'rows', 'stable');
    thisObjXY   = thisObjXY(~isnan(thisObjXY(:, 1)), :);
    myp         = plot(thisObjXY(1, 1), thisObjXY(1, 2), 'ko', 'MarkerSize', 4);
    % preferred object
    if iObj == dt.prefObjTest.prefObjName
        set(myp, 'MarkerFaceColor', rgb('orange'));
    end
end
% boundary
alpha = 0:0.001:(2*pi);
plot(cos(alpha) .* 5000, sin(alpha) .* 5000, 'k.', 'markerSize', 3);
t1 = text(0, 6500, 'x');
t2 = text(-6500, 0, 'y');
t3 = text(0, -7400, {'Object locations'});
set([t1, t2, t3], 'fontunits', 'centimeters', 'fontsize', myFontSize, 'HorizontalAlignment', 'center');
hold off;
set(gca, 'ydir', 'reverse', 'xlim', [-5000, 5000], 'ylim', [-5000, 5000]);
axis off square equal;

% waveform
ax5 = axes('units', 'centimeters', 'position', [1, 4, 1.5, 1.5]);
outPlot = LK_SpikeDensityPlot(dt.thisSpike, dt.sr);
colormap(ax5, 'parula');
set(gca, ...
    'xlim', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), 'xtick', round([min(outPlot.spikeTime), max(outPlot.spikeTime)]), ...
    'ylim', [outPlot.yMin, outPlot.yMax], 'ytick', [outPlot.yMin, outPlot.yMax], 'ticklength', [0, 0]);
xl = xlabel('ms', 'units', 'normalized', 'position', [0.5, -0.1, 0]);
yl = ylabel('\muV', 'units', 'normalized', 'position', [-0.12, 0.5, 0]);
% indicate number of spikes used for the analysis
tx = text(0.5, 1.15, ['n=', num2str(dt.nspk)], 'units', 'normalized', 'HorizontalAlignment', 'center');
set([gca, xl, yl, tx], 'FontUnits', 'centimeters', 'FontSize', myFontSize);

%% output

% create output
out         = [];
out.f       = f;
out.m       = m;
out.ste     = ste;
dt.label   = {'unpreferred'; 'preferred'};
