function idx = LK_log2idx(logVec)
%
% LK_log2idx converts a logical vector into an n x 2 matrix with indices
% that denote the start and end samples of true values in the logical
% vector.
%
% Lukas Kunz, 2021

% check whether input is logical
if ~islogical(logVec)
    error('Input is not logical.');
end

% check whether input has the correct dimensions
if size(logVec, 1) > size(logVec, 2)
    logVec  = transpose(logVec);
end

% convert logical vector into indices
idx = [transpose(find(diff([0, logVec]) == 1)), ...
    transpose(find(diff([logVec, 0]) == -1))];