function ste = LK_ste(mat)
%
% LK_ste computes the standard error for each column of the input matrix
% "mat".
%
% The calculation is performed across the first dimension of the input
% matrix ("mat"): ste = std(mat, [], 1, 'omitnan') ./ sqrt(sum(~isnan(mat), 1));
%
% Use as: ste = LK_ste(mat);
%
% See also: std, nanstd
%
% Lukas Kunz, 2021

ste = std(mat, [], 1, 'omitnan') ./ sqrt(sum(~isnan(mat), 1));