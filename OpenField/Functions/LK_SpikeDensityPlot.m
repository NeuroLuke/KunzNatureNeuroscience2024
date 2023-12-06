function out = LK_SpikeDensityPlot(wavs, sr)
%
% LK_SpikeDensityPlot plots the spike waveforms of a given unit as a
% density plot.
%
% Input:
%   wavs: n x m matrix with waveforms from n spikes
%   sr: sampling rate
%
% Output:
%   spikeDensity: p x m matrix with waveform density
%   spikeTime: 1 x m vector with time axis
%   yCenters: p x 1 vector with y-values
%   yMin: minimum y-value
%   yMax: maximum y-value
%
% Function adapted from Reber et al., PLOS Biology, 2019.
%
% Lukas Kunz, 2021

% x-axis: time
spikeIdx    = 1:size(wavs, 2);
spikeTime   = (spikeIdx - 1) ./ sr .* 1000; % in milliseconds

% y-axis: voltages
yMin        = floor(min(wavs(:))); % lowest point of the figure
yMax        = ceil(max(wavs(:))); % highest point of the figure
yRes        = 150; % number of bins in the y-direction
yCenters    = transpose(linspace(yMin, yMax, yRes));
yBinSize    = yCenters(2) - yCenters(1);
yEdges      = [yCenters - yBinSize / 2, yCenters + yBinSize / 2];

% preallocate 2D histogram that contains bin counts
spikeDensity    = zeros(numel(yCenters), numel(spikeTime)); % voltage * time (e.g., 150 * 64)
% loop through spikeTime bins
for iTime = 1:numel(spikeIdx)
    % loop through voltage-bins
    for iY = 1:numel(yCenters)
        spikeDensity(iY, iTime) = sum(wavs(:, iTime) <= yEdges(iY, 2) & wavs(:, iTime) > yEdges(iY, 1));
    end
end

% normalization
spikeDensity    = spikeDensity ./ max(spikeDensity(:));

% remove extreme outliers in order to keep color resolution
cutoff                                  = 5 * std(spikeDensity(:)); % cutoff for values that are too high
spikeDensity(spikeDensity > cutoff)     = cutoff; % replace n with n without too high bin values

% plot the 2D histogram
pcolor(spikeTime, yCenters, spikeDensity);
shading interp;

%% output
out                 = [];
out.spikeDensity    = spikeDensity;
out.spikeTime       = spikeTime;
out.yCenters        = yCenters;
out.yMin            = yMin;
out.yMax            = yMax;

