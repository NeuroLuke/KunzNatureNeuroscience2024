function outMatrix = LK_resizem(inMatrix, augFac)
%
% LK_resizem resizes matrix "inMatrix" using the augmentation factor
% augFac. inMatrix is an n x m matrix.
%
% Output is the resized matrix "outMatrix".
%
% Lukas Kunz, 2021

% size of input matrix
[n, m]      = size(inMatrix);

% initialize output matrix
outMatrix   = cell(n, m);

% fill output matrix with values from 2D input matrix
for iN = 1:size(inMatrix, 1)
    for iM = 1:size(inMatrix, 2)
        outMatrix{iN, iM}   = repmat(inMatrix(iN, iM), augFac, augFac);
    end
end

% unfold output matrix
outMatrix   = cell2mat(outMatrix);