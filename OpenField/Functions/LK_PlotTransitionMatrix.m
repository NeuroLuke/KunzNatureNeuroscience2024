function f = LK_PlotTransitionMatrix(tra)
%
% LK_PlotTransitionMatrix plots the transition matrix for converting
% original channels into bipolar channels.
%
% The output is a figure handle.
%
% Lukas Kunz, 2021

% create figure
f = figure('units', 'centimeters', 'position', [-5, -5, size(tra.labelold, 1), size(tra.labelnew, 1)], 'visible', 'off');
imagesc(tra.tra);
set(gca, ...
    'xtick', 1:length(tra.labelold), 'xticklabel', strrep(tra.labelold, '_', '\_'), 'xticklabelrotation', 90, ...
    'ytick', 1:length(tra.labelnew), 'yticklabel', strrep(tra.labelnew, '_', '\_'), ...
    'FontUnits', 'centimeters', 'FontSize', 0.25);
xl = xlabel('Original label');
yl = ylabel('New label');
tl = title('Re-referencing scheme for bipolar montage');
set([xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.5, 'FontWeight', 'normal');