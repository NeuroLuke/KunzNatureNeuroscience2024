function f = LK_HistogramFigure_20220321(dt)
%
% LK_HistogramFigure_20220321 creates a histogram figure.
%
% Input is a structure with multiple fields:
%   figPosition
%   axPosition
%   data
%   binEdges (optional)
%   normalization (optional)
%   xline (optional)
%   xlabel
%   ylabel
%   title (optional)
%   fontSize
%
% Output is a figure handle, f.
%
% Lukas Kunz, 2022

% create figure
f = figure('units', 'centimeters', 'position', dt.figPosition);
axes('units', 'centimeters', 'position', dt.axPosition);
hold on;
% histogram
if isfield(dt, 'binEdges')
    h = histogram(dt.data, dt.binEdges);
else
    h = histogram(dt.data);
end
% mean/median
if isfield(dt, 'xline')
    xline(dt.xline, 'r', 'LineWidth', 2);
end
% enhance the histogram
set(h, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
if isfield(dt, 'normalization')
    set(h, 'Normalization', dt.normalization);
end
% enhance the axes
xl = xlabel(dt.xlabel);
yl = ylabel(dt.ylabel);
if isfield(dt, 'title')
    tl = title(dt.title);
else
    tl = '';
end
set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', dt.fontSize, 'fontweight', 'normal');
set(gca, 'box', 'off', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
% bring axes to front
set(gca, 'Layer', 'Top');
