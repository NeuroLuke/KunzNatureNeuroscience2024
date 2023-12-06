function [] = LK_print(f, filename, formattype, resolution)
%
% LK_print prints a figure and is the analog to Matlab's print function. It
% sets the PaperPositionMode to auto.
%
% Use as: LK_print(f, filename, formattype, resolution);
%
% Lukas Kunz, 2021

set(f, 'PaperPositionMode', 'auto');
print(f, filename, formattype, resolution);