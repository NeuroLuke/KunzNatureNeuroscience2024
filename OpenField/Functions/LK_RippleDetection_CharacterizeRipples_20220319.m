function out = LK_RippleDetection_CharacterizeRipples_20220319(mycfg, myeeg)
%
% LK_RippleDetection_CharacterizeRipples_20220319 characterizes ripples.
%
% Use as: out = LK_RippleDetection_CharacterizeRipples_20220319(mycfg, myeeg);
%
% Input #1 is a structure with a matrix that contains the start and end
% indices of the ripples (i.e., an m x 2 matrix).
%
% Input #2 is a fieldtrip-like structure with fields:
%   label     	--> channel name;
%   fsample   	--> sampling frequency;
%   trial     	--> original data;
%   time     	--> original time;
%   index     	--> sample indices of the data;
%   dataB      	--> bandpass-filtered data;
%   dataEB    	--> smoothed envelope of the bandpass-filtered data, with
%                   nans during artifacts.
%
% Output is a structure "ripples" with multiple types of information about
% the ripples.
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Characterization of ripples (channel: %s) ...\n', myeeg.label{:});

%% characterize candidate ripples

% loop through candidates and characterize them
ripples(size(mycfg.idxRipples, 1), 1)   = struct();
for iR = 1:size(mycfg.idxRipples, 1)
    
    % data for this candidate
    thisIdxRipples  = mycfg.idxRipples(iR, 1):mycfg.idxRipples(iR, 2);
    thisIdx         = myeeg.index{1}(thisIdxRipples); % sample indices
    thisTime        = myeeg.time{1}(thisIdxRipples); % time points
    thisDataB       = myeeg.dataB{1}(thisIdxRipples); % bandpass signal
    thisDataEB      = myeeg.dataEB{1}(thisIdxRipples); % smoothed envelope of the bandpass signal
    
    % start, middle, and end sample index
    ripples(iR).startIdx        = min(thisIdx);
    ripples(iR).midIdx          = round(mean(thisIdx));
    ripples(iR).endIdx          = max(thisIdx);
    
    % start, middle, and end time
    ripples(iR).startTime       = min(thisTime);
    ripples(iR).midTime         = mean(thisTime);
    ripples(iR).endTime         = max(thisTime);
    
    % duration
    ripples(iR).duration        = ripples(iR).endTime - ripples(iR).startTime;
    
    % global peak of the band-pass signal
    [peakAmp, peakIdx]          = max(thisDataB, [], 2);
    ripples(iR).peakAmp         = peakAmp;
    ripples(iR).peakIdx         = thisIdx(peakIdx);
    ripples(iR).peakTime        = thisTime(peakIdx);
    
    % global trough of the band-pass signal
    [troughAmp, troughIdx]      = min(thisDataB, [], 2);
    ripples(iR).troughAmp       = troughAmp;
    ripples(iR).troughIdx       = thisIdx(troughIdx);
    ripples(iR).troughTime      = thisTime(troughIdx);
    
    % global peak of the envelope signal
    [peakAmp, peakIdx]          = max(thisDataEB, [], 2);
    ripples(iR).envPeakAmp      = peakAmp;
    ripples(iR).envPeakIdx      = thisIdx(peakIdx);
    ripples(iR).envPeakTime     = thisTime(peakIdx);
    
    % global trough of the envelope signal
    [troughAmp, troughIdx]   	= min(thisDataEB, [], 2);
    ripples(iR).envTroughAmp    = troughAmp;
    ripples(iR).envTroughIdx    = thisIdx(troughIdx);
    ripples(iR).envTroughTime   = thisTime(troughIdx);
    
    % sum of the envelope signal
    ripples(iR).envSum          = sum(thisDataEB, 2);
    
    % local peaks of the band-pass signal
    if numel(thisDataB) > 2
        [~, peaks]              	= findpeaks(thisDataB);
        if ~isempty(peaks)
            ripples(iR).peaksIdx    = thisIdx(peaks);
            ripples(iR).peaksTime   = thisTime(peaks);
        else
            ripples(iR).peaksIdx    = nan;
            ripples(iR).peaksTime   = nan;
        end
    else
        ripples(iR).peaksIdx        = nan;
        ripples(iR).peaksTime       = nan;
    end
    
    % local troughs of the band-pass signal
    if numel(thisDataB) > 2
        [~, troughs]                	= findpeaks(-1 .* thisDataB);
        if ~isempty(troughs)
            ripples(iR).troughsIdx      = thisIdx(troughs);
            ripples(iR).troughsTime     = thisTime(troughs);
        else
            ripples(iR).troughsIdx      = nan;
            ripples(iR).troughsTime     = nan;
        end
    else
        ripples(iR).troughsIdx          = nan;
        ripples(iR).troughsTime         = nan;
    end
    
    % ripple-specific frequency of peaks and troughs
    if all(~isnan(ripples(iR).peaksIdx)) && all(~isnan(ripples(iR).troughIdx))
        ripples(iR).freq    = myeeg.fsample / mean([diff(ripples(iR).peaksIdx), diff(ripples(iR).troughsIdx)]);
    else
        ripples(iR).freq    = nan;
    end
end

%% create output

% output structure
out         = [];
out.ripples = ripples;

% report
fprintf('\n====================================================== Done: characterization of ripples.\n');
