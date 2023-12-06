function outVector = LK_resizev(inVector, augFac)
%
% LK_resizev resizes vector "inVector" using the augmentation factor
% "augFac". inVector is an n x 1 or 1 x n vector.
%
% Output is the resized vertical vector "outVector".
%
% Use as: outVector = LK_resizev(inVector, augFac)
%
% Lukas Kunz, 2021

% size of input vector
[n, m]          = size(inVector);

% sanity check
if ~(n == 1 || m == 1)
    error("The input is not a vector.");
end

% flip vector if it does not have the correct dimensions
if n < m
    inVector    = transpose(inVector);
end

% initialize output vector
outVector       = cell(size(inVector));

% fill output vector with values from input vector
for i = 1:size(inVector, 1)
    outVector{i, 1} = repmat(inVector(i, 1), augFac, 1);
end

% unfold output vector
outVector   = cell2mat(outVector);