function d = LK_dir(myPath)
%
% LK_dir functions as MATLAB's "dir" but removes '.' and '..' directories
% from the output.
%
% Use as: d = LK_dir(myPath);
%
% See also dir
%
% Lukas Kunz, 2021

% directory information
d   = dir(myPath);

% remove '.' and '..' paths
bRemove     = strcmp({d.name}, '.') | strcmp({d.name}, '..');
d           = d(~bRemove);