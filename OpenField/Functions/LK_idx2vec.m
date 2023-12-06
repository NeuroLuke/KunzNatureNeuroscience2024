function vec = LK_idx2vec(mat)
%
% LK_idx2vec converts an m x 2 matrix of start and end indices (first and
% second column of the matrix, respectively) into an 1 x n vector
% containing all indices between the respective start and end indices.
%
% LK_idx2vec is different from LK_idx2log.
%
% Lukas Kunz, 2021

% check input
if size(mat, 2) ~= 2
    error('The size of the input is not correct.');
end

% convert matrix into vector
vec = cell2mat(transpose(arrayfun(@(x, y) x:y, round(mat(:, 1)), round(mat(:, 2)), 'UniformOutput', false)));

% check output: is the length of the complete index vector reasonable based
% on the original start-end indices?
if size(vec, 2) ~= (sum(diff(mat, [], 2) + 1))
    error('The size of the output is not correct.');
end