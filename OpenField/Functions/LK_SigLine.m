function LK_SigLine(x, y, p, varargin)
%
% LK_SigLine plots a significance line.
%
% Input:
%   x: all x values
%   y: two y values that indicate the position of the horizontal lines
%   p: p-values (if float, then the significance level is indicated via
%       different shades of gray; if logical, all significant time points
%       are plotted in black)
%
% Lukas Kunz, 2021

% sanity checks
if size(x, 1) >= size(x, 2)
    error('X data has the wrong dimensions.');
end
if size(y, 2) >= size(y, 1)
    error('Y data has the wrong dimensions.');
end
if size(p, 1) > size(p, 2)
    p   = transpose(p); % same dimensions as x data
end

% if p-values are specified as logical
if islogical(p)
    pNew                = ones(size(p));
    pNew(p == true)     = 0.0005; % set p-value below 0.001 to mark in black
    p                   = pNew;
end

% spacing between x values
xSpacing    = min(diff(x));

% differentiate between different alpha levels using different shades of
% gray
alphaLevels = {0.05, [0.6, 0.6, 0.6]; ...
    0.01, [0.3, 0.3, 0.3]; ...
    0.001, [0, 0, 0]};

% loop through alpha levels
for iAlpha = 1:size(alphaLevels, 1)
    
    % significant areas
    [L, NUM]    = bwlabel(p < alphaLevels{iAlpha, 1});
    
    % loop through significant areas
    for iNUM = 1:NUM
        
        % prepare patch
        xPatchForth     = [min(x(L == iNUM)) - xSpacing / 2, max(x(L == iNUM)) + xSpacing / 2];
        xPatchBack      = fliplr(xPatchForth);
        yPatchLow       = repmat(min(y), 1, numel(xPatchForth));
        yPatchHigh      = fliplr(repmat(max(y), 1, numel(xPatchBack)));
        
        % plot patch
        myP = patch([xPatchForth, xPatchBack], [yPatchLow, yPatchHigh], alphaLevels{iAlpha, 2}, ...
            'EdgeColor', 'none');
        
        % if you have specified a color for the patch
        if ~isempty(varargin)
            set(myP, 'FaceColor', varargin{1});
        end
        uistack(myP, 'top');
    end
    
    % plot rest in white
    [L, NUM]    = bwlabel(L == 0);
    
    % loop through unsignificant areas
    for iNUM = 1:NUM
        
        % prepare patch
        xPatchForth     = [min(x(L == iNUM)) - xSpacing / 2, max(x(L == iNUM)) + xSpacing / 2];
        xPatchBack      = fliplr(xPatchForth);
        yPatchLow       = repmat(min(y), 1, numel(xPatchForth));
        yPatchHigh      = fliplr(repmat(max(y), 1, numel(xPatchBack)));
        
        % plot patch
        myP = patch([xPatchForth, xPatchBack], [yPatchLow, yPatchHigh], [1, 1, 1], ...
            'EdgeColor', 'none');
        uistack(myP, 'top');
    end
end

% plot lower limit
myYL = yline(min(y), 'Color', [0, 0, 0]);
% if you have specified a color for the patch
if ~isempty(varargin)
    set(myYL, 'Color', varargin{1});
end

% plot upper limit
myYL = yline(max(y), 'Color', [0, 0, 0]);
% if you have specified a color for the patch
if ~isempty(varargin)
    set(myYL, 'Color', varargin{1});
end
