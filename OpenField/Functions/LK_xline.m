function xl = LK_xline(x, myLineStyle, myColor, myStackLevel)
%
% LK_xline creates a constant line at a specified x value.
%
% LK_xline is similar to MATLAB's xline, but allows to stack the line to
% the bottom.
%
% Use as: xl = LK_xline(x, myLineStyle, myColor, myStackLevel)
%
% Lukas Kunz, 2021

% information about current axes
tmpAx = get(gca);

% plot xline
xl = plot([x, x], [min(tmpAx.YLim), max(tmpAx.YLim)], myLineStyle, ...
    'Color', myColor);

% stack
uistack(xl, myStackLevel);