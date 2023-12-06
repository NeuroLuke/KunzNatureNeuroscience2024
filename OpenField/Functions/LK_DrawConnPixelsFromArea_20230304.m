function drawnPixels = LK_DrawConnPixelsFromArea_20230304(A, n)
%
% LK_DrawConnPixelsFromArea_20230304 draws n connected pixels from within
% the area of ones of the 2D matrix A.
%
% It can be used to create surrogate place fields, for example.
%
% Use as: drawnPixels = LK_DrawConnPixelsFromArea_20230304(A, n)
%
% Input:
%   A: 2D matrix of zeros and ones, where the ones indicate the possible
%      locations of the drawn pixels (e.g., a place-field map).
%   n: number of connected pixels to draw.
%
% Output:
%   drawnPixels: connected pixels randomly drawn from the area of ones
%                inside matrix A.
%
% Lukas Kunz, 2023

% check whether A has been binarized
if numel(unique(A(:))) ~= 2 || ~all(unique(A(:)) == [0; 1])
    error('The input matrix has not been binarized.');
end

% find the indices of the pixels within the area of ones in A
A_idx                   = find(A);

% initialize a mask of pixels to draw
drawnPixels             = zeros(size(A));

% choose a starting pixel at random within the area of ones in A
rand_idx                = datasample(A_idx, 1);
drawnPixels(rand_idx)   = 1;

% while there are fewer than n pixels in the drawn pixels mask, draw
% additional pixels
numIter = 0;
while nnz(drawnPixels) < n
    
    % keep track of the number of iterations
    numIter = numIter + 1;
    if numIter > (2 * numel(A))
        drawnPixels = [];
        warning('Drawing contiguous pixels from area: no solution found.');
        return;
    end
    
    % find the indices of the pixels in the drawn pixels mask
    drawn_idx                   = find(drawnPixels);
    [drawn_rows, drawn_cols]    = ind2sub(size(A), drawn_idx);
    
    % possible next pixels to draw
    poss_rows   = [drawn_rows - 1; drawn_rows + 1; drawn_rows + 0; drawn_rows + 0]; % above, below, this, this
    poss_cols   = [drawn_cols + 0; drawn_cols + 0; drawn_cols - 1; drawn_cols + 1]; % this, this, left, right
    bValid      = (poss_rows >= 1) & (poss_rows <= size(A, 1)) & (poss_cols >= 1) & (poss_cols <= size(A, 2));
    poss_idx    = sub2ind(size(A), poss_rows(bValid), poss_cols(bValid));
    
    % restrict possible pixels: should be part of the ones in A and should
    % not have already been drawn
    poss_idx    = poss_idx(ismember(poss_idx, A_idx) & ~ismember(poss_idx, drawn_idx));
    
    % skip if there are no possible indices
    if isempty(poss_idx)
        continue;
    end
    
    % add one possible pixel to the drawn pixels
    chosen_idx              = datasample(poss_idx, 1);
    drawnPixels(chosen_idx) = 1;
end
