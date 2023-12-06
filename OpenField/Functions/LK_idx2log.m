function logVec = LK_idx2log(idx, num)
%
% LK_idx2log converts an n x 2 matrix with indices that denote the start
% and end samples of true values into a logical vector of length num.
%
% Use as: logVec = LK_idx2log(idx, num);
%
% Lukas Kunz, 2021

% create a logical vector of length num
logVec  = false(1, num);

% if there are no indices, return
if isempty(idx)
    return;
end

% convert idx matrix into idx vector
idxVec  = cell2mat(transpose(arrayfun(@(x, y)  x:y, idx(:, 1), idx(:, 2), 'UniformOutput', 0)));

% set relevant entries of the logical vector to true
logVec(idxVec)  = true;