function f = LK_PlotRipplePSTHPerPhaseStart_20231030(dt)
%
% LK_PlotRipplePSTHPerPhaseStart_20231030 plots the ripple PSTH and
% time-resolved average ripple rate per trial phase. Good and bad trials
% are separated. Data is locked to the start of the trial phases.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 10]);

% ripple rate in each time bin, averaged within each channel, good vs. bad
axes('units', 'centimeters', 'position', [1.25, 0.9, 4.5, 2.6]);
hold on;
% mean and standard error across channels
m = {mean(dt.meanPSTHBad, 1, 'omitnan'), mean(dt.meanPSTHGood, 1, 'omitnan')}; % average across channels
ste = {LK_ste(dt.meanPSTHBad), LK_ste(dt.meanPSTHGood)}; % standard error across channels
for iM = 1:numel(m)
    pa = patch([dt.beh.PSTHStartTime, fliplr(dt.beh.PSTHStartTime)], ...
        [m{iM} + ste{iM}, fliplr(m{iM} - ste{iM})], dt.beh.trialPhasesColor(dt.iPhase, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    p = plot(dt.beh.PSTHStartTime, m{iM}, 'Color', dt.beh.trialPhasesColor(dt.iPhase, :));
    if iM == 1
        set(pa, 'FaceColor', [0.5, 0.5, 0.5]);
        set(p, 'Color', [0.5, 0.5, 0.5]);
    end
end
% plot significance
sigTime = zeros(1, numel(dt.beh.PSTHStartTime));
sigTime(ismember(dt.beh.PSTHStartTime, dt.outFT.time))  = dt.outFT.mask;
LK_SigLine(dt.beh.PSTHStartTime, [max(dt.thisParam.yLim); max(dt.thisParam.yLim) - 0.05 * range(dt.thisParam.yLim)], sigTime == 1);
% enhance axes
set(gca, ...
    'xlim', [min(dt.thisParam.xLim), max(dt.thisParam.xLim)], 'xtick', [min(dt.thisParam.xLim), max(dt.thisParam.xLim)], ...
    'ylim', dt.thisParam.yLim, 'ytick', dt.thisParam.yLim, ...
    'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
LK_xline(0, '--', [0, 0, 0], 'bottom');
xl = xlabel('Time (s)', 'units', 'normalized', 'position', [0.5, -0.14]);
yl = ylabel({'Ripple rate', '(Hz)'}, 'units', 'normalized', 'position', [-0.05, 0.5]);
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);

% ripple PSTH
axes('units', 'centimeters', 'position', [1.25, 4, 4.5, 5.4]);
hold on;
for iTrial = 1:size(dt.allPSTH, 1)
    myp = plot(dt.beh.PSTHStartTime(dt.allPSTH(iTrial, :) > 0), repmat(iTrial, 1, sum(dt.allPSTH(iTrial, :) > 0)), '.', ...
        'Color', dt.beh.trialPhasesColor(dt.iPhase, :), 'MarkerSize', 2);
    if dt.allBGoodMemory(iTrial) == 0
        set(myp, 'Color', [0.5, 0.5, 0.5]); % bad trials in gray
    end
end
LK_xline(0, '--', [0, 0, 0], 'bottom');
% enhance axes
yl = ylabel('Trial', 'units', 'normalized', 'position', [-0.05, 0.5]);
tl = title(dt.beh.trialPhases{dt.iPhase});
set(gca, ...
    'xlim', dt.thisParam.xLim, 'xtick', [], 'ylim', [1, size(dt.allPSTH, 1)], 'ytick', [1, size(dt.allPSTH, 1)], ...
    'tickdir', 'out', 'ticklength', [0.01, 0.01], 'box', 'off', 'layer', 'top');
set([gca, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
drawnow;
