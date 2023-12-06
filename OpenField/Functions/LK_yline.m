function yl = LK_yline(y, myLineStyle, myColor, myStackLevel)
%
% LK_yline creates a constant line at a specified y value.
%
% LK_yline is similar to MATLAB's yline, but allows to stack the line to
% the bottom.
%
% Use as: yl = LK_yline(y, myLineStyle, myColor, myStackLevel)
%
% See also: yline, LK_xline
%
% Lukas Kunz, 2021

% information about current axes
tmpAx = get(gca);

% plot yline
yl = plot([min(tmpAx.XLim), max(tmpAx.XLim)], [y, y], myLineStyle, ...
    'Color', myColor);

% stack
uistack(yl, myStackLevel);