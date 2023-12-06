function f = LK_PlotRippleCorrelationsWithBeh_20231030(dt)
%
% LK_PlotRippleCorrelationsWithBeh_20231030 plots the correlations between
% ripple rates and different behavioral properties including memory
% performance and trial index.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 7.5, 6]);
axes('units', 'centimeters', 'position', [2.5, 1.4, 4.5, 3.5]);
hold on;

% loop through the different trial phases
for iPhase = 1:size(dt.myData, 2)
    
    % report mean and standard error
    thisM = mean(dt.myData(:, iPhase), 1, 'omitnan');
    thisSTE = LK_ste(dt.myData(:, iPhase));
    fprintf('\t\tPhase: %s, mean = %.3f, STE = %.3f.\n', dt.beh.trialPhases{iPhase}, thisM, thisSTE);

    % plot outliers
    bOutlier = isoutlier(dt.myData(:, iPhase), 'quartiles'); % elements more than 1.5 IQRs above the upper quartile or below the lower quartile
    if any(bOutlier)
        plot(dt.myData(bOutlier, iPhase), iPhase, '.', 'Color', dt.beh.trialPhasesColor(iPhase, :));
    end

    % plot minimum and maximum
    thisMin = min(dt.myData(~bOutlier, iPhase));
    thisMax = max(dt.myData(~bOutlier, iPhase));
    plot([thisMin, thisMin, thisMin, thisMax, thisMax, thisMax], iPhase + [-0.2, 0.2, 0, 0, -0.2, 0.2], '-', 'Color', dt.beh.trialPhasesColor(iPhase, :));

    % 25th and 75th percentile
    thisCI = prctile(dt.myData(:, iPhase), [25, 75]);
    patch([thisCI(1), thisCI(1), thisCI(2), thisCI(2), thisCI(1)], iPhase + [-0.4, 0.4, 0.4, -0.4, -0.4], dt.beh.trialPhasesColor(iPhase, :), 'FaceAlpha', 0.5, ...
        'EdgeColor', dt.beh.trialPhasesColor(iPhase, :), 'LineWidth', 1);
    
    % median
    thisMedian = median(dt.myData(:, iPhase), 'omitnan');
    plot([thisMedian, thisMedian], iPhase + [-0.4, 0.4], '-', 'LineWidth', 1, 'Color', dt.beh.trialPhasesColor(iPhase, :));

    % significance info
    if dt.p_corr2Beh(iPhase) < 0.001
        sigText = '***';
    elseif dt.p_corr2Beh(iPhase) < 0.01
        sigText = '**';
    elseif dt.p_corr2Beh(iPhase) < 0.05
        sigText = '*';
    else
        sigText = '';
    end
    text(max(dt.dataGroup{2}), iPhase, sigText, 'Color', [0, 0, 0], ...
        'FontUnits', 'centimeters', 'FontSize', 0.5, 'HorizontalAlignment', 'center', 'Rotation', 90);
end

% enhance axes
set(gca, ...
    'xlim', dt.dataGroup{2}, 'xtick', [min(dt.dataGroup{2}), 0, max(dt.dataGroup{2})], ...
    'ylim', [0.4, numel(dt.beh.trialPhases) + 0.6], 'ytick', 1:numel(dt.beh.trialPhases), ...
    'yticklabel', dt.beh.trialPhases, 'ydir', 'reverse', ...
    'xticklabelrotation', 0, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
if ~(strcmp(dt.corrGroup, 'chanCorr2Memory'))
    set(gca, 'yticklabel', '');
end
xl = xlabel('Pearson''s \itr');
tl = title(dt.dataGroup{3}, 'units', 'normalized', 'position', [0.5, 1.025]);
set([gca, xl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
