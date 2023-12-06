%==========================================================================
% This script provides a simulation to illustrate the coactivity z-score.
%
% References: Sosa et al., Neuron, 2020.
%
% Lukas Kunz, 2022
%==========================================================================

% start
clear; close all; clc;

% paths
paths           = [];
paths.functions = 'E:\OpenField\Functions\';
paths.save      = 'E:\OpenField\NeuronalCoactivityDuringRipples_20211025\CoactivityIllustration\';

% add paths
addpath(genpath(paths.functions));
mkdir(paths.save);

%% coactivity scores

% counts
N = 10; % number of ripples
nA = transpose(0:N); % number of ripples in which cell A is active
nB = transpose(0:N); % number of ripples in which cell B is active
nAB = transpose(0:N); % number of ripples in which cell A and B are active

% preallocate and loop through combinations
z = nan(size(nA, 1), size(nB, 1), size(nAB, 1));
for iA = 1:size(nA, 1)
    for iB = 1:size(nB, 1)
        for iAB = 1:size(nAB)
            
            % compute coactivity z-score
            try
                z(iA, iB, iAB)  = LK_CoactivityZScore_20220615(N, nA(iA), nB(iB), nAB(iAB));
            catch
            end
            
            % set inf and -inf to nan
            if z(iA, iB, iAB) == Inf || z(iA, iB, iAB) == -Inf
                z(iA, iB, iAB)  = nan;
            end
        end
    end
end

%% figure

% create figure and extract video frames
f = figure('units', 'centimeters', 'position', [2, 2, 6.2, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.3, 4, 4]);
myCLim = [-4, 4];
for iAB = 1:size(z, 3)
    
    % create imagesc of z-values
    imagesc(nB, nA, z(:, :, iAB), ...
        'AlphaData', ~isnan(z(:, :, iAB)));
    xl = xlabel('# ripples cell B active');
    yl = ylabel('# ripples cell A active');
    tl = title(sprintf('# ripples both cells active: %d', nAB(iAB)), ...
        'units', 'normalized', 'position', [0.5, 1.04]);
    cb = colorbar('units', 'centimeters', 'position', [5.4, 2.5, 0.25, 1]);
    cb.Ticks = myCLim;
    caxis(myCLim);
    colormap(parula);
    title(cb, 'z');
    set(gca, 'xlim', [min(nB) - 0.5, max(nB) + 0.5], 'xtick', min(nB):2:max(nB), 'xticklabelrotation', 0, ...
        'ylim', [min(nA) - 0.5, max(nA) + 0.5], 'ytick', min(nA):2:max(nA));
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    grid minor;
    
    % extract frame data
    F(iAB) = getframe(gcf);
    
    % draw and save current figure and clear window
    drawnow;
    LK_print(f, strcat(paths.save, 'LK_SmallIllustration_zScoreCoactivity_', num2str(iAB), '_20230330'), '-dtiff', '-r300');
    cla;
end

%% video

% create the video writer with 1 fps
writerObj = VideoWriter(strcat(paths.save, 'LK_SmallIllustration_zScoreCoactivity_20230330.avi'));
writerObj.FrameRate = 1;

% open the video writer and write the frames to the video
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);
