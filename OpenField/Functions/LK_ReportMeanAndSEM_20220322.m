function LK_ReportMeanAndSEM_20220322(msg, data)
%
% LK_ReportMeanAndSEM_20220322 reports the mean and standard error of the
% mean of an input vector.
%
% Use as:
% LK_ReportMeanAndSEM_20220322(msg, data);
%
% Input:
%   msg     --> message;
%   data    --> data as n x 1 vector.
%
% Lukas Kunz, 2022

% compute mean and standard error along the first dimension
n   = sum(~isnan(data), 1);
m   = mean(data, 1, 'omitnan');
ste = std(data, [], 1, 'omitnan') / sqrt(n);

% report
fprintf(sprintf('%s: mean = %.3f, SEM = %.3f (n = %d).\n', msg, m, ste, n));