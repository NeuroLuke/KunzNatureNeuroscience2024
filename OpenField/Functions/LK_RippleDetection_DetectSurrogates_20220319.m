function out = LK_RippleDetection_DetectSurrogates_20220319(mycfg, myeeg)
%
% LK_RippleDetection_DetectSurrogates_20220319 finds surrogate ripples.
%
% Use as: out = LK_RippleDetection_DetectSurrogates_20220319(mycfg, myeeg);
%
% Input #1 is a structure with parameter settings, including:
%   bArtifact   --> timepoints that are part of artifacts;
%   ripples     --> information about "real" ripples.
%
% Input #2 is a fieldtrip structure with fields
%   label       --> channel name;
%   fsample     --> sampling frequency;
%   trial       --> original data;
%   time        --> timepoints of the data;
%   index       --> sample indices of the data;
%   dataB       --> bandpass-filtered data;
%   dataEB      --> smoothed envelope of the bandpass-filtered data, with
%                   artifacts removed.
%
% Output is an m x 2 matrix ("idxSurrogate") with start and end indices of
% surrogate ripples.
%
% References:
%   Ngo et al., eLife, 2020
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Detection of surrogate ripples (channel: %s) ...\n', myeeg.label{:});

%% samples available for extracting surrogate ripples

% check data size
if size(mycfg.bArtifact, 2) ~= size(mycfg.ripples.bRipple, 2)
    error('The boolean artifact vector has a different size as the boolean ripple vector.');
end

% check data plausibility
if any(mycfg.bArtifact & mycfg.ripples.bRipple)
    error('Artifacts and ripples are overlapping.');
end

% identify time epochs that can be used for extracting surrogate ripples
bPossible      = ~any([mycfg.bArtifact; mycfg.ripples.bRipple], 1);
idxPossible    = LK_log2idx(bPossible);

% report
fprintf('Initial number of data segments from which to draw the surrogate ripples: %d.\n', size(idxPossible, 1));

%% surrogate ripples

% loop through "real" ripples and extract corresponding surrogate ripples
idxSurrogate    = nan(size(mycfg.ripples.ripples, 1), 2); % start index, end index
for iR = 1:size(mycfg.ripples.ripples, 1)
    
    %% define start index of this surrogate ripple
    
    % start and duration of this empirical ripple
    thisStartIdx        = mycfg.ripples.ripples(iR).startIdx; % (samples)
    thisDurationIdx     = mycfg.ripples.ripples(iR).endIdx - mycfg.ripples.ripples(iR).startIdx; % (samples)
    
    % shift the end samples of the possible periods to earlier samples
    % based on the duration of this empirical ripple ==> surrogate ripples
    % always end within the available time period
    thisIdxPossible     = [idxPossible(:, 1), idxPossible(:, 2) - thisDurationIdx];
    
    % remove negative periods
    thisIdxPossible   	= thisIdxPossible(thisIdxPossible(:, 1) < thisIdxPossible(:, 2), :);
    
    % convert n x 2 matrix into 1 x m vector
    thisIdxPossible     = LK_idx2vec(thisIdxPossible);
    
    % only consider indices that are relatively close in time to the "real"
    % ripples
    bClose              = abs(thisIdxPossible - thisStartIdx) <= (mycfg.surro.timeWindow * myeeg.fsample);
    
    % sample randomly from the possible indices that are close enough to
    % obtain the start sample of this surrogate ripple
    surroStartIdx       = datasample(thisIdxPossible(bClose), 1); % sample one index
    
    % estimate the end sample for this surrogate ripple based on the
    % duration of the corresponding real ripple
    surroEndIdx         = surroStartIdx + thisDurationIdx;
    
    % store all start and end samples of surrogate ripples
    idxSurrogate(iR, :) = [surroStartIdx, surroEndIdx];
    
    %% update the samples available for extracting surrogate ripples
    
    % the time window for this surrogate ripple is no longer available
    thisIdxSurro            = surroStartIdx:surroEndIdx;
    bPossible(thisIdxSurro) = false;
    idxPossible             = LK_log2idx(bPossible);
end

%% output

% create output
out                 = [];
out.idxSurrogate    = idxSurrogate; % start and end indices of the surrogate candidate ripples

% report
fprintf('\nNumber of surrogate ripples: %d.\n', size(idxSurrogate, 1));
fprintf('====================================================== Done: detection of surrogate ripples.\n');
