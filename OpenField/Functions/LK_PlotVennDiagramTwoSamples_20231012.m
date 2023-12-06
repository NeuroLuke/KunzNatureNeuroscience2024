function f = LK_PlotVennDiagramTwoSamples_20231012(dt)
%
% LK_PlotVennDiagramTwoSamples_20231012 plots a Venn diagram with two
% samples.
%
% Input is a structure, dt, with fields A (total areas) and I
% (intersection). Orange and blue are used to color the total areas.
%
% Output is a figure handle.
%
% Lukas Kunz, 2023

% create figure
f = figure('units', 'centimeters', 'position', [2, 5, 12, 8]);
venn(dt.A, dt.I, 'FaceColor', {rgb('orange'), rgb('blue')}, 'LineWidth', 2);
axis equal tight off;