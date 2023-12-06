function bMask = LK_MaskBehavioralSampling(cfg)
%
% LK_MaskBehavioralSampling creates a mask based on the number of discrete
% observations per factor level. Each factor is checked separately and
% independently from the other factors.
%
% Input is a structure with fields
%   data        --> time x n matrix (for n factors to check)
%   factorNames --> factor names
%   minNumObs   --> minimum number of observations per factor level
%
% Lukas Kunz, 2021

%% process data

% ensure that you have enough observations per factor level
factorLevels    = cell(size(cfg.factorNames)); % e.g., size = 3 x 1; each cell entry contains the unique group-levels for each group
bMaskAll        = nan(size(cfg.data, 1), size(cfg.data, 2)); % timepoint-specific mask; e.g., size = 10,000 x 3
for iFactor = 1:size(cfg.data, 2)
    
    % unique levels for this factor
    bNan            = isnan(cfg.data(:, iFactor));
    uniqueLevels    = unique(cfg.data(~bNan, iFactor));
    
    % number of discrete observations per factor level
    numObs          = nan(size(uniqueLevels, 1), 1);
    for iLevel = 1:size(uniqueLevels, 1)
        
        % number of discrete observations for this factor level
        [~, thisNumObs]     = bwlabel(cfg.data(:, iFactor) == uniqueLevels(iLevel, 1));
        numObs(iLevel, 1)   = thisNumObs;
    end
    
    % exclude levels for which the number of observations is too low
    factorLevels{iFactor} 	= uniqueLevels(numObs >= cfg.minNumObs);
    bMaskAll(:, iFactor)   	= ismember(cfg.data(:, iFactor), factorLevels{iFactor});
end

%% mask

% create a mask for the data: enough observations and no nans
bMask   = all(bMaskAll, 2) & all(~isnan(cfg.data), 2);
