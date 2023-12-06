function z = LK_CoactivityZScore_20220615(varargin)
%
% LK_CoactivityZScore_20220615 estimates a coactivity z-score between pairs
% of cells during ripples.
% 
% Use as: z = LK_CoactivityZScore_20220615(activityCellA, activityCellB)
%     or: z = LK_CoactivityZScore_20220615(N, nA, nB, nAB)
% 
% Input version #1:
%   activityCellA = number-of-ripples x 1 vector containing whether the
%                   first cell was active during the ripples
%                   1 = active; 0 = inactive
%   activityCellB = number-of-ripples x 1 vector containing whether the
%                   second cell was active during the ripples
%                   1 = active; 0 = inactive
%
% Input version #2:
%   N   = number of ripples
%   nA  = number of ripples in which cell A was active
%   nB  = number of ripples in which cell B was active
%   nAB = number of ripples in which cell A and cell B were both active
%
% Output:
%   z-score
%
% Notes:
% - If >= 1 cell did not spike during any ripple, the result is a NaN.
% - If both cells spiked during all ripples, the result is a NaN.
%
% References:
%   Cheng & Frank, Neuron, 2008
%   Singer & Frank, Neuron, 2009
%   Sosa et al., Neuron, 2020
%
% Lukas Kunz, 2022

% check the number of input arguments
narginchk(2, 4); % minimum and maximum number of inputs

% process input
if nargin == 2 % if cellular activity is given as two input vectors

    % cellular activity of cell A and B
    activityCellA   = varargin{1};
    activityCellB   = varargin{2};
    
    % estimate the different quantities
    N   = size(activityCellA, 1); % total number of ripples
    nA  = sum(activityCellA); % number of ripples in which cell A is active
    nB  = sum(activityCellB); % number of ripples in which cell B is active
    nAB = sum(activityCellA & activityCellB); % number of ripples in which both cell A and cell B are active

elseif nargin == 4 % if cellular activity is given as counts

    % extract the different quantities
    N   = varargin{1};
    nA  = varargin{2};
    nB  = varargin{3};
    nAB = varargin{4};
end

% sanity checks
if nAB > nA || nAB > nB
    error('nAB is larger than nA or nB.');
end
if (nA + nB) > (N + nAB)
    error('nA + nB is larger than N + nAB.')
end

% specify numerator (difference between observed and expected coactivities)
numerator       = nAB - (nA * nB) / N;

% specify denumerator (for normalization by the standard deviation)
denominatorTop  = nA * nB * (N - nA) * (N - nB);
denominatorLow  = (N .^ 2) * (N - 1);
denominator     = sqrt(denominatorTop / denominatorLow); % take the sqrt to obtain the standard deviation (instead of variance)

% calculate z-score
z = numerator / denominator;
