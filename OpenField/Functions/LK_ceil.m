function outNum = LK_ceil(inNum, precision)
%
% LK_ceil is identical to MATLAB's ceil with the additional option to
% specify a precision to ceil towards a specific decimal place. The
% precision is specified as a power of 10 (with positive precision values
% leading to more precise ceiling).
%
% Lukas Kunz, 2021

prec    = 10 .^ precision;
outNum  = ceil(inNum .* prec) ./ prec;