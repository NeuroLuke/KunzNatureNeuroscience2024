function out = LK_ArtifactDetection_20220319(mycfg, myeeg)
%
% LK_ArtifactDetection_20220319 detects timepoint-specific artifacts based
% on several criteria.
%
% Use as: out = LK_ArtifactDetection_20220319(mycfg, myeeg);
%
% Input #1 is a structure with optional analysis settings:
%   time2exclude    --> time to exclude around each artifact;
%   amp            	--> artifact detection based on amplitudes;
%   gra           	--> artifact detection based on gradients;
%   hfa          	--> artifact detection based on HFA;
%   pow          	--> artifact detection based on power (summed across
%                       multiple frequency bands).
% 
% Input #2 is a fieldtrip structure containing data from one channel:
%   label        	--> channel label;
%   fsample       	--> sampling frequency;
%   time          	--> trial time;
%   trial       	--> trial data.
%
% Output is a structure with fields
%   artifactMat             --> number-of-artifacts x number-of-timepoints
%                               matrix with artifact indices;
%   bArtifact               --> 1 x number-of-timepoints boolean vector
%                               (derived from artifactMat);
%   artifactOnsEndsTime     --> number-of-artifacts x 2 matrix (seconds);
%   artifactOnsEndsSample   --> number-of-artifacts x 2 matrix (samples);
%   artifactType            --> number-of-artifacts x artifact-type matrix.
%
% References:
%   Staresina et al., Nat Neurosci, 2015;
%   Ngo et al., eLife, 2020.
%
% Lukas Kunz, 2022

% report
fprintf('\n====================================================== Performing artifact detection (channel: %s) ...\n', myeeg.label{:});

%% broadband filtering of the data

% apply low high-pass filter to get rid of drifts
cfg             = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.hpfiltord   = 5;
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 150;
bbData          = ft_preprocessing(cfg, myeeg); % broadband data

%% artifact detection using amplitude threshold
% rationale: epileptic spikes exhibit large amplitudes

% if configuration structure contains amplitude settings
if isfield(mycfg, 'amp')
    
    % report
    fprintf('\n------------------------------------------------------ using amplitudes.\n');
    
    % estimate mean and standard deviation of values
    dataCenter      = median(bbData.trial{1}, 2, 'omitnan');
    dataSpread      = iqr(bbData.trial{1}, 2);
    
    % amplitude threshold
    ampThresh       = [dataCenter - dataSpread * mycfg.amp.threshFac, dataCenter + dataSpread * mycfg.amp.threshFac];
    
    % find values exceeding the threshold
    bAmpArtifact    = bbData.trial{1} < min(ampThresh) | bbData.trial{1} > max(ampThresh);
end

%% artifact detection using gradient threshold
% rationale: epileptic spikes exhibit sharp changes in amplitude

% if configuration structure contains gradient settings
if isfield(mycfg, 'gra')
    
    % report
    fprintf('\n------------------------------------------------------ using gradients.\n');
    
    % estimate gradients
    gradients       = diff(bbData.trial{1}, [], 2);
    gradients       = [gradients, median(gradients, 'omitnan')]; % add a median gradient at the end
    
    % estimate median and inter-quartile-range of gradients
    gradientsCenter = median(gradients, 2, 'omitnan');
    gradientsSpread = iqr(gradients, 2);
    
    % gradient threshold
    graThresh       = [gradientsCenter - gradientsSpread * mycfg.gra.threshFac, gradientsCenter + gradientsSpread * mycfg.gra.threshFac];
    
    % find gradients exceeding the threshold
    bGraArtifact    = gradients < min(graThresh) | gradients > max(graThresh);
end

%% artifact detection using absolute ampitude after applying a high-pass-filter of 250Hz
% rationale: remove pathological high-frequency oscillations

% if configuration structure contains HFA settings
if isfield(mycfg, 'hfa')
    
    % report
    fprintf('\n------------------------------------------------------ using HFA amplitudes.\n');
    
    % high-pass filter
    cfg             = [];
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = mycfg.hfa.freq;
    hfData          = ft_preprocessing(cfg, myeeg);
    
    % estimate median and inter-quartile-range
    hfCenter        = median(hfData.trial{1}, 2, 'omitnan');
    hfSpread        = iqr(hfData.trial{1}, 2);
    
    % high-frequency amplitude threshold
    hfThresh        = [hfCenter - hfSpread * mycfg.hfa.threshFac, hfCenter + hfSpread * mycfg.hfa.threshFac];
    
    % find values exceeding the threshold
    bHFArtifact     = hfData.trial{1} > max(hfThresh); % only consider significantly increased HFA
end

%% artifact detection using TF power
% rationale: epileptic spikes exhibit power increases across a broad
% frequency spectrum

% if configuration structure contains power settings
if isfield(mycfg, 'pow')

    % report
    fprintf('\n------------------------------------------------------ using power.\n');
    
    % compute time-frequency power spectrum
    cfg             = [];
    cfg.method      = 'wavelet';
    cfg.width       = 7;
    cfg.toi         = 'all';
    cfg.foi         = logspace(log10(mycfg.pow.lowerFreq), log10(mycfg.pow.upperFreq), mycfg.pow.numFreqs); % n logarithmically spaced frequencies
    cfg.pad         = 'nextpow2';
    cfg.output      = 'pow';
    tfData          = ft_freqanalysis(cfg, myeeg);
    
    % log-transformation and frequency-specific z-scoring
    squeezedPower   = squeeze(tfData.powspctrm); % freq x time
    logPower        = log(squeezedPower); % freq x time
    zPower          = normalize(logPower, 2, 'zscore'); % normalize across time; ignore nans
    
    % sum z-scored power across frequencies
    zPowerSum       = sum(zPower, 1, 'omitnan'); % sum across freqs
    
    % estimate center and spread
    powerCenter     = median(zPowerSum, 2, 'omitnan');
    powerSpread     = iqr(zPowerSum, 2);
    
    % define power threshold
    powerThresh     = powerCenter + powerSpread * mycfg.pow.threshFac;
    
    % find values exceeding the threshold
    bPowerArtifact  = zPowerSum > powerThresh; % only consider significantly increased power
end

%% identification of all artifacts (1 = artifact, 0 = no artifact)

% combine the different artifact types
artifactMat     = [];
if isfield(mycfg, 'amp')
    artifactMat = cat(1, artifactMat, bAmpArtifact .* mycfg.amp.idx);
end
if isfield(mycfg, 'gra')
    artifactMat = cat(1, artifactMat, bGraArtifact .* mycfg.gra.idx);
end
if isfield(mycfg, 'hfa')
    artifactMat = cat(1, artifactMat, bHFArtifact .* mycfg.hfa.idx);
end
if isfield(mycfg, 'pow')
    artifactMat = cat(1, artifactMat, bPowerArtifact .* mycfg.pow.idx);
end

% identify samples with any type of artifact
bArtifact   = sum(artifactMat, 1) > 0; % along first dimension ==> 1 x timepoints

%% exclusion of additional time around each artifact

% report
fprintf('\n------------------------------------------------------ excluding additional %.3f seconds around each artifact.\n', mycfg.time2Exclude);

% exclude +/- n seconds around each artifact
convVecLength      	= mycfg.time2Exclude * myeeg.fsample * 2; % multiply by 2 because convolution vector is centered on midpoint
if mod(convVecLength, 2) == 0
    convVecLength  	= convVecLength + 1; % length of convolution vector be odd
end
bArtifact        	= conv(bArtifact, ones(1, convVecLength), 'same') > 0;

%% convert the boolean artifact vector into onset- and end-timepoints

% clusters of contiguous timepoints identified as an artifact
[L, NUM]    = bwlabel(bArtifact);

% loop through artifact clusters
artifactOnsEndsTime     = nan(NUM, 2);
artifactOnsEndsSample   = nan(NUM, 2);
numArt                  = max([size(artifactMat, 1); max(artifactMat(:))]);
artifactType            = nan(NUM, numArt); % number-of-artifacts x types-of-artifacts
for iNUM = 1:NUM
    
    % onsets and ends (seconds)
    artifactOnsEndsTime(iNUM, :)    = [min(myeeg.time{1}(L == iNUM)), max(myeeg.time{1}(L == iNUM))];
    
    % onsets and ends (samples)
    artifactOnsEndsSample(iNUM, :)  = [find(L == iNUM, 1, 'first'), find(L == iNUM, 1, 'last')];
    
    % type of artifact
    tmpIdx                      = artifactMat(:, L == iNUM);
    tmpIdx                      = unique(tmpIdx(tmpIdx ~= 0));
    artifactType(iNUM, tmpIdx)  = tmpIdx;
end

%% report

% percentage of samples removed because of artifact
percArtifact    = 100 * sum(bArtifact, 2) / size(bArtifact, 2);
fprintf('\n====================================================== Done. Percentage of artifact samples: %.2f%%.\n\n', percArtifact);

%% output

% create output
out                         = [];
out.artifactMat             = artifactMat; % types-of-artifacts x time
out.bArtifact               = bArtifact; % 1 x time
out.artifactOnsEndsTime     = artifactOnsEndsTime; % number-of-artifacts x types-of-artifacts
out.artifactOnsEndsSample   = artifactOnsEndsSample;
out.artifactType            = artifactType;
