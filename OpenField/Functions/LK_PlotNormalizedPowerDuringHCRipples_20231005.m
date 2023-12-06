function f = LK_PlotNormalizedPowerDuringHCRipples_20231005(dt)
%
% LK_PlotNormalizedPowerDuringHCRipples_20231005 plots normalized power
% during hippocampal ripples.
%
% Input is a structure with multiple fields.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% loop through LFP-channel regions
x = dt.foi; % frequencies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'Color', [1, 1, 1]);
axes('units', 'centimeters', 'position', [1.5, 0.9, 4.25, 3], 'box', 'off', 'xlim', [min(x), max(x)]);
hold on;
LK_yline(0, '--', rgb('gray'), 'bottom');
p = nan(size(dt.LFPROIs, 1), 1);
for iRegion = 1:size(dt.LFPROIs, 1)

    % select channels from this region
    bMask   = ismember(dt.allLFPChanRegion, dt.LFPROIs{iRegion, 2});

    % assemble data for this figure
    m      	= mean(dt.thisSpectra(bMask, :), 1, 'omitnan'); % average across channels
    ste   	= LK_ste(dt.thisSpectra(bMask, :)); % standard error across channels

    % plot mean and standard error
    patch([x, fliplr(x)], [m - ste, fliplr(m + ste)], dt.LFPROIs{iRegion, 3}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    p(iRegion) = plot(x, m, 'Color', dt.LFPROIs{iRegion, 3}, 'LineWidth', 1);
end
set(gca, 'xscale', 'log', 'tickdir', 'out');
xl = xlabel('Frequency (Hz)', 'units', 'normalized', 'position', [0.44, -0.11]);
yl = ylabel('Normalized power');
tl = title({'Non-HC regions', ['(', num2str(sum(contains(dt.allLFPChanRegion, {'AMY', 'EC', 'PHC', 'TP'}))), ' combinations)']}, ...
    'units', 'normalized', 'position', [0.5, 1.38]);
lg = legend(p, dt.LFPROIs(:, 1), 'box', 'off', 'location', 'eastoutside', 'fontsize', 10, 'units', 'centimeters', 'position', [1.5, 2.85, 2, 2]);
lg.ItemTokenSize = [7, 3];
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
