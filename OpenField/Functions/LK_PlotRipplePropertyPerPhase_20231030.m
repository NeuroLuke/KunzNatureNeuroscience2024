function f = LK_PlotRipplePropertyPerPhase_20231030(dt)
%
% LK_PlotRipplePropertyPerPhase_20231030 plots a particular ripple property
% by trial phase.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 5, 6]);
axes('units', 'centimeters', 'position', [2.4, 1.4, 2, 3.5]);
hold on;

% loop through the different trial phases
for iPhase = 1:numel(dt.beh.trialPhases)
    
    % report mean and standard error
    thisM = mean(dt.myData(:, iPhase), 1, 'omitnan');
    thisSTE = LK_ste(dt.myData(:, iPhase));
    fprintf('\t\tPhase: %s, mean = %.3f, SEM = %.3f.\n', dt.beh.trialPhases{iPhase}, thisM, thisSTE);
    
    % plot outliers
    bOutlier = isoutlier(dt.myData(:, iPhase), 'quartiles'); % elements more than 1.5 IQRs above the upper quartile or below the lower quartile
    if any(bOutlier)
        plot(dt.myData(bOutlier, iPhase), iPhase, '.', 'Color', dt.beh.trialPhasesColor(iPhase, :));
    end
    
    % plot minimum and maximum
    thisMin = min(dt.myData(~bOutlier, iPhase)); % remove outliers from calculating the minimum and maximum
    thisMax = max(dt.myData(~bOutlier, iPhase));
    plot([thisMin, thisMin, thisMin, thisMax, thisMax, thisMax], iPhase + [-0.2, 0.2, 0, 0, -0.2, 0.2], '-', 'Color', dt.beh.trialPhasesColor(iPhase, :));
    
    % 25th and 75th percentile
    thisCI = prctile(dt.myData(:, iPhase), [25, 75]);
    patch([thisCI(1), thisCI(1), thisCI(2), thisCI(2), thisCI(1)], iPhase + [-0.4, 0.4, 0.4, -0.4, -0.4], dt.beh.trialPhasesColor(iPhase, :), 'FaceAlpha', 0.5, ...
        'EdgeColor', dt.beh.trialPhasesColor(iPhase, :), 'LineWidth', 1);
    
    % median
    thisMedian = median(dt.myData(:, iPhase), 'omitnan');
    plot([thisMedian, thisMedian], iPhase + [-0.4, 0.4], '-', 'LineWidth', 1, 'Color', dt.beh.trialPhasesColor(iPhase, :));
end

% enhance axes
set(gca, ...
    'xlim', [LK_floor(min(dt.myData(:)) * 0.98, dt.myPrec), LK_ceil(max(dt.myData(:)) * 1.02, dt.myPrec)], ...
    'xtick', [LK_floor(min(dt.myData(:)) * 0.98, dt.myPrec), LK_ceil(max(dt.myData(:)) * 1.02, dt.myPrec)], ...
    'ylim', [0.4, numel(dt.beh.trialPhases) + 0.6], 'ytick', 1:numel(dt.beh.trialPhases), 'yticklabel', dt.beh.trialPhases, 'ydir', 'reverse', ...
    'xticklabelrotation', 0, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
if ~ismember(dt.group, {'chanRippleRatePerPhase', 'chanFracArtifactsPerPhase', 'chanFracIEDsPerPhase'})
    set(gca, 'yticklabel', '');
end
xl = xlabel(dt.myXLabel);
tl = title(dt.myTitle, 'units', 'normalized', 'position', [0.5, 1.025]);
set([gca, xl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
