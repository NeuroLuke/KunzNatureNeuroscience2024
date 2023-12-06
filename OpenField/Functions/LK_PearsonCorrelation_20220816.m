function r = LK_PearsonCorrelation_20220816(v1, v2)
%
% LK_PearsonCorrelation_20220816 estimates a Pearson correlation between
% the two vectors v1 and v2. This implementation is faster than Matlab's
% corr function.
%
% Use as: r = LK_PearsonCorrelation_20220816(v1, v2);
%
% Input:
%   v1  --> n x 1 vector;
%   v2  --> n x 1 vector.
%
% Output:
%   r   --> Pearson's r value.
%
% See also: corr
% 
% Lukas Kunz, 2022

% sanity checks
[r1, c1]    = size(v1);
[r2, c2]    = size(v2);
if r1 ~= r2 || c1 ~= c2
    error('The input vectors v1 and v2 have different sizes.');
end
if c1 ~= 1
    error('The input vectors v1 and v2 are not one-dimensional.');
end

% remove nans
bValid  = ~isnan(v1) & ~isnan(v2);
v1      = v1(bValid);
v2      = v2(bValid);

% estimate the different components
N       = size(v1, 1); % sample size
v1Mean  = mean(v1); % mean of the first vector
v1Std   = std(v1); % standard deviation of the first vector
v2Mean  = mean(v2); % mean of the second vector
v2Std   = std(v2); % standard deviation of the second vector

% estimate the different subterms of the equation
NTerm   = 1 / (N - 1);
v1Term  = (v1 - v1Mean) / v1Std;
v2Term  = (v2 - v2Mean) / v2Std;

% estimate Pearson correlation
r       = NTerm * sum(v1Term .* v2Term);