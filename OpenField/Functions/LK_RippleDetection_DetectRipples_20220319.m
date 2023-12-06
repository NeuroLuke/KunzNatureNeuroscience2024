function out = LK_RippleDetection_DetectRipples_20220319(mycfg, myeeg)
%
% LK_RippleDetection_DetectRipples_20220319 detects ripples.
% Step 1: candidate ripples are detected.
% Step 2: candidates are selected as ripples based on several criteria:
%   - peak amplitude;
%   - duration;
%   - number of cycles;
%   - power spectrum.
%
% Candidate ripples that do not meet the power-spectrum criterion are
% classified as "false positives".
%
% Use as: out = LK_RippleDetection_DetectRipples_20220319(mycfg, myeeg);
%
% Input #1 is a structure with multiple parameter settings.
%
% Input #2 is a fieldtrip-like structure with fields:
%   label       --> channel name;
%   fsample  	--> sampling frequency;
%   trial     	--> original data;
%   time     	--> original time of the data;
%   index     	--> sample indices of the data;
%   dataB     	--> bandpass-filtered data;
%   dataEB    	--> smoothed envelope of the bandpass-filtered data, with
%                   NaNs during artifacts.
%
% Output is an m x 2 matrix with the onset and end indices of the ripples
% ("idxRipples") and an m x 2 matrix with the onset and end indices of
% false positives ("idxFalsePositives").
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Detection of ripples (channel: %s) ...\n', myeeg.label{:});

%% center and spread of the envelope ==> thresholds

% report
fprintf('\n------------------------------------------------------ Center, spread, and thresholds ...\n');

% average and standard deviation based on the smoothed envelope
center  = mean(myeeg.dataEB{1}, 2, 'omitnan');
spread  = std(myeeg.dataEB{1}, [], 2, 'omitnan');

% thresholds
threshMin       = center + spread .* mycfg.threshMinFac; % minimum threshold
threshMax       = center + spread .* mycfg.threshMaxFac; % maximum threshold
threshPeakMin   = center + spread .* mycfg.threshPeakMinFac; % minimum peak threshold

%% candidate ripples

% report
fprintf('\n------------------------------------------------------ Candidate ripples ...\n');

% find threshold crossings in the smoothed envelope
bCandidate      = myeeg.dataEB{1} >= threshMin & myeeg.dataEB{1} <= threshMax;

% start and end index of the candidate ripples
idxCandidates   = round(LK_log2idx(bCandidate));

% report
fprintf('Initial number of candidate ripples: %d.\n', size(idxCandidates, 1));

%% selection criterion: peak amplitude

% report
fprintf('\n------------------------------------------------------ Selection criterion: peak amplitude ...\n');

% only keep candidate ripples that cross the minimum threshold for the peak
bValidPeak  = false(size(idxCandidates, 1), 1);
for iC = 1:size(idxCandidates, 1)
    thisIdx = idxCandidates(iC, 1):idxCandidates(iC, 2);
    if max(myeeg.dataEB{1}(thisIdx)) >= threshPeakMin
        bValidPeak(iC, 1)   = true;
    end
end

% update information about candidates
idxCandidates   = idxCandidates(bValidPeak, :);

% report
fprintf('Number of candidate ripples: %d.\n', size(idxCandidates, 1));

%% selection criterion: duration

% report
fprintf('\n------------------------------------------------------ Selection criterion: duration ...\n');

% only keep candidate ripples that meet the necessary duration
bValidDur   = false(size(idxCandidates, 1), 1);
for iC = 1:size(idxCandidates, 1)
    thisDur = myeeg.time{1}(idxCandidates(iC, 2)) - myeeg.time{1}(idxCandidates(iC, 1)); % duration of this candidate
    if thisDur >= mycfg.durationMin && thisDur <= mycfg.durationMax
        bValidDur(iC, 1)    = true; % if the duration is in the desired range
    end
end

% restrict to valid candidates
idxCandidates   = idxCandidates(bValidDur, :);

% report
fprintf('Number of candidate ripples: %d.\n', size(idxCandidates, 1));

%% selection criterion: number of cycles in the bandpass signal

% report
fprintf('\n------------------------------------------------------ Selection criterion: number of cycles ...\n');

% remove candidates not fulfilling the number-of-cycles requirement
bValidCycle = false(size(idxCandidates, 1), 1);
for iC = 1:size(idxCandidates, 1)
    
    % data for this candidate
    thisIdx     = idxCandidates(iC, 1):idxCandidates(iC, 2);
    thisDataB	= myeeg.dataB{1}(thisIdx); % bandpass signal
    
    % local peaks and troughs of the band-pass signal
    if numel(thisDataB) > 2
        
        % peaks
        [~, peaks]  	= findpeaks(thisDataB);

        % troughs
        [~, troughs]    = findpeaks(-1 .* thisDataB);
        
        % test whether the number of peaks and troughs is sufficient
        if numel(peaks) >= mycfg.numCyclesMin && numel(troughs) >= mycfg.numCyclesMin
            bValidCycle(iC, 1)  = true;
        end
    end    
end

% restrict to valid candidates
idxCandidates   = idxCandidates(bValidCycle, :);

% report
fprintf('Number of candidate ripples: %d.\n', size(idxCandidates, 1));

%% selection criterion: power spectrum

% report
fprintf('\n------------------------------------------------------ Selection criterion: power spectrum ...\n');

% time-frequency power spectrum of the entire data
cfg         = [];
cfg.method  = mycfg.falsePos.method;
cfg.width   = mycfg.falsePos.width;
if strcmp(mycfg.falsePos.timeRes, 'all')
    cfg.toi = mycfg.falsePos.timeRes;
else
    cfg.toi = min(myeeg.time{1}):mycfg.falsePos.timeRes:max(myeeg.time{1});
end
cfg.foi     = mycfg.falsePos.foi;
cfg.pad     = mycfg.falsePos.pad;
cfg.output  = mycfg.falsePos.output;
tfData      = ft_freqanalysis(cfg, myeeg);

% power at all time points
pow       	= squeeze(tfData.powspctrm); % freq x time

% grand-average power (average across time)
meanPow   	= mean(pow, 2, 'omitnan'); % freq x 1

% normalize power by dividing by grand-average power
normPow   	= pow ./ repmat(meanPow, 1, size(pow, 2));

% loop through candidates and check validity of power spectrum
bValidPow   = false(size(idxCandidates, 1), 1);
for iC = 1:size(idxCandidates, 1)
    
    % sample indices and time of this candidate
    thisIdx         = idxCandidates(iC, 1):idxCandidates(iC, 2);
    thisTime        = myeeg.time{1}(thisIdx);
    
    % normalized power of this candidate
    logIdx          = tfData.time >= min(thisTime) & tfData.time <= max(thisTime);
    thisNormPow     = normPow(:, logIdx); % normalized power for this candidate
    thisMeanNormPow = mean(thisNormPow, 2, 'omitnan'); % average across time
    
    % peaks in the power spectrum of this candidate
    [peaksPow, peaksFreq]   = findpeaks(thisMeanNormPow, mycfg.falsePos.foi); % peaksFreq is in Hz
    
    % index of maximum peak in the power spectrum
    [~, maxPeakIdx]         = max(peaksPow);
    
    % test whether maximum peak is within the desired ripple band
    if ~isempty(maxPeakIdx)
        if peaksFreq(maxPeakIdx) >= min(mycfg.bpfreq) && peaksFreq(maxPeakIdx) <= max(mycfg.bpfreq)
            bValidPow(iC) = true;
        end
    end
end

%% indices of ripples and false positives

% ripple indices
idxRipples          = idxCandidates(bValidPow, :);

% false-positive indices (= candidates without a valid power spectrum)
idxFalsePositives   = idxCandidates(~bValidPow, :);

%% create output

% output structure
out                     = [];
out.idxRipples          = idxRipples; % start and end indices of the candidate ripples
out.idxFalsePositives   = idxFalsePositives; % start and end indices of the false positives

% report
fprintf('\nNumber of ripples: %d.\n', size(idxRipples, 1));
fprintf('Number of false positives: %d.\n', size(idxFalsePositives, 1));
fprintf('====================================================== Done: selection of ripples.\n');
