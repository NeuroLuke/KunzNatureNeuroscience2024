function dropError = LK_CalculateDropError(x1, x2, y1, y2)
%
% LK_CalculateDropError calculates the drop error given the response
% locations (x1, y1) and the correct locations (x2, y2).
%
% Use as: dropError = LK_CalculateDropError(x1, x2, y1, y2);
% 
% Lukas Kunz, 2022

dropError = sqrt((x1 - x2) .^ 2 + (y1 - y2) .^ 2);