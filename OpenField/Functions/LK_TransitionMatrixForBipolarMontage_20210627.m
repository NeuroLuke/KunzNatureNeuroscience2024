function out = LK_TransitionMatrixForBipolarMontage_20210627(labelold)
%
% LK_TransitionMatrixForBipolarMontage_20210627 computes the transition
% matrix to convert the EEG data into a bipolar montage.
%
% The transition matrix follows the fieldtrip standard: the rows are the
% new labels, the columns are the old labels.
%
% Lukas Kunz, 2021

%% sanity checking

% cave: replace dashes in the old labels with an "X"
if any(contains(labelold, 'X'))
    error('An electrode label contains an "X".'); % ensure that "X" is not already used in any channel label
end

% replace dashes
if any(contains(labelold, '-'))
    labelold    = strrep(labelold, '-', 'X');
    repFlag     = '-';
end

% replace underscores
if any(contains(labelold, '_'))
    labelold    = strrep(labelold, '_', 'X');
    repFlag     = '_';
end

%% create new labels based on the old labels

% unique electrodes
uniqueElec  = unique(replace(labelold, digitsPattern, ''), 'stable');

% loop through unique elctrodes
labelnew    = [];
for iUE = 1:size(uniqueElec, 1)
    
    % channels from this electrode
    logIdx    	= strcmp(replace(labelold, digitsPattern, ''), uniqueElec{iUE});
    thisChans   = labelold(logIdx);
    
    % sort in ascending order
    [~, I]      = sort(cellfun(@str2double, replace(thisChans, lettersPattern, '')));
    thisChans   = thisChans(I);
    
    % loop through channels from this electrode
    for iChan = 1:size(thisChans, 1) - 1
        
        % contributing channels
        chan1   = thisChans{iChan};
        chan2   = thisChans{iChan + 1};
        
        % check whether they are neighboring channels
        if (str2double(replace(chan1, lettersPattern, '')) + 1) ~= str2double(replace(chan2, lettersPattern, ''))
            fprintf('\tSkipping this channel from the same electrode because the channels are not neighbors: %s - %s.\n', chan1, chan2);
            continue;
        end
        
        % new label for bipolar channel
        labelnew    = cat(1, labelnew, {strcat(chan1, '-', chan2)});
    end
end

%% create transition matrix

% transition matrix
tra     = zeros(size(labelnew, 1), size(labelold, 1)); % new bipolar electrodes * original electrodes
for iNewChan = 1:size(labelnew, 1)
    
    % old channels composing the new channel
    labelsOld12             = split(labelnew{iNewChan}, '-');
    
    % minuend
    idxOld1                 = strcmp(labelsOld12{1}, labelold);
    tra(iNewChan, idxOld1)  = 1;
    
    % subtrahend
    idxOld2                 = strcmp(labelsOld12{2}, labelold);
    tra(iNewChan, idxOld2)  = -1;
end

%% sanity checking

% cave: re-replace the "X" with dashes
if any(contains(labelold, 'X'))
    
    % report
    fprintf('\tReplacing the "X" in "labelold" and "labelnew" with its original sign (for %s).\n', labelold{contains(labelold, 'X')});
    
    % replace
    labelold    = strrep(labelold, 'X', repFlag);
    labelnew    = strrep(labelnew, 'X', repFlag);
end

%% output

% create output structure
out             = [];
out.tra         = tra;
out.labelold    = labelold;
out.labelnew    = labelnew;
