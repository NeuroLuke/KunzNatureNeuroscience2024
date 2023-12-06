%==========================================================================
% This script examines whether ripples are phase-locked to slow
% oscillations.
%
% References: Logothetis et al., Nature, 2012.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;

% settings
param                       = [];
param.myRNG                 = 1111;
param.numSurrogates         = 1001; % number of surrogates
%- for subject information
param.info.name             = 'subjectdata_20210616';
%- for type of preprocessing
param.preproc.type          = 'FBS'; % filtered, bipolar, select
param.preproc.fileName      = strcat('eeg', param.preproc.type, '.mat');
%- for ripples
param.ripple.bpfreq         = [80, 140];
param.ripple.threshMinFac   = 2;
%- for channels
param.channel.choi          = {'HCleftBipolar'; 'HCrightBipolar'}; % channels of interest
%- for phase estimation using Hilbert transformation
param.phase.freqBand        = [0.5, 2]; % frequency band
param.phase.filtOrdFac      = 2; % filter-order factor
%- for power
param.power.numLevels       = 3; % number of power levels to differentiate between

% paths
paths           = [];
paths.subjects  = 'E:\OpenField\SubjectData_20210616\';
paths.macro     = 'E:\OpenField\MacroPreprocessing_20210627\';
paths.ripple  	= strcat('E:\OpenField\RipplesMacro_20210614\20220320_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '\');
paths.save      = strcat('E:\OpenField\RipplesMacroPhaseLocking_20210806\20230521_', param.preproc.type, '_', ...
    num2str(min(param.ripple.bpfreq)), 'to', num2str(max(param.ripple.bpfreq)), 'Hz_tMF', num2str(param.ripple.threshMinFac), '\', ...
    num2str(min(param.phase.freqBand)), 'to', num2str(max(param.phase.freqBand)), 'Hz\');
paths.fieldtrip = 'E:\fieldtrip\fieldtrip-20210614\';
paths.functions = 'E:\OpenField\Functions\';

% add paths
addpath(genpath(paths.functions));
addpath(paths.fieldtrip);
ft_defaults;
ft_warning('off');
mkdir(paths.save);

% subjects
subjects    = load(strcat(paths.subjects, 'subjects.mat'));
subjects    = subjects.subjects(~strcmp(subjects.subjects(:, 2), 'Excluded'), :);

%% preallocate

% save settings
save(strcat(paths.save, 'settings'));

% main results for each channel
allRes      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % report and reset rng
    fprintf('\nSUBJECT: %s.\n\n', subjects{iSub});
    rng(param.myRNG);
    
    %% channels of interest
    
    % load subject information
    subjectdata = load(fullfile(paths.subjects, subjects{iSub}, param.info.name));
    subjectdata = subjectdata.subjectdata;
    
    % channels of interest in this subject
    choi        = [];
    if any(strcmp(param.channel.choi, 'HCleftBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCleftBipolar);
    end
    if any(strcmp(param.channel.choi, 'HCrightBipolar'))
        choi    = cat(1, choi, subjectdata.macro.HCrightBipolar);
    end
    
    % skip if this subject does not have channels of interest
    if isempty(choi)
        warning('This subject does not have any channels of interest.');
        continue;
    end
    
    %% eeg data
    
    % load preprocessed iEEG data for estimating the phase of slow
    % oscillations
    if strcmp(param.preproc.type, 'FBS')
        eegAll  = load(fullfile(paths.macro, subjects{iSub}, param.preproc.fileName));
        eegAll  = eegAll.eegFBS;
    end
    
    %% loop through channels of interest
    for iChan = 1:numel(choi)
        
        %% eeg and ripple data from this channel of interest
        
        % preprocessed eeg data from this channel
        cfg         = [];
        cfg.channel = choi{iChan};
        eeg         = ft_selectdata(cfg, eegAll);
        
        % ripple data from this channel
        ripples     = load(fullfile(paths.ripple, subjects{iSub}, choi{iChan}, 'ripples.mat'));
        ripples     = ripples.ripples;
        
        % sanity check
        if numel(ripples.bRipple) ~= numel(eeg.trial{1})
            error('Problem with the congruence of the eeg and the ripple data (regarding their size).');
        end
        
        % report
        fprintf('\n\n------------------------------------------------------\n');
        fprintf('Working on "%s", channel "%s".\n', subjects{iSub}, eeg.label{:});
        fprintf('------------------------------------------------------\n\n');
        
        %% LFP phase and power estimation using Matlab's filter design
        
        % bandpass filtering
        bpFilter = designfilt('bandpassfir', 'FilterOrder', param.phase.filtOrdFac * (eeg.fsample / min(param.phase.freqBand)), ...
            'CutoffFrequency1', min(param.phase.freqBand), 'CutoffFrequency2', max(param.phase.freqBand), 'SampleRate', eeg.fsample);
        bpData  = filtfilt(bpFilter, eeg.trial{1});
        
        % phase and power
        phase       = angle(hilbert(bpData - mean(bpData, 'omitnan')));
        powerData   = abs(bpData) .^ 2;
        
        %% extract phases at ripples
        
        % phase at ripples
        rippleIdx       = cell2mat({ripples.ripples.peakIdx}');
        ripplePhases    = transpose(phase(rippleIdx));
        
        % mean phase at ripples
        meanRipplePhase = circ_mean(ripplePhases);
        
        % Rayleigh's test
        [Rp, Rz]        = circ_rtest(ripplePhases);
        
        %% surrogate statistics using permutations of the inter-ripple-intervals
        
        % inter-ripple-intervals (in indices)
        IRI     = diff([0; rippleIdx]);
        
        % loop through surrogates
        meanRipplePhaseSurro    = nan(param.numSurrogates, 1);
        RzSurro                 = nan(param.numSurrogates, 1);
        for iSurro = 1:param.numSurrogates
            
            % phases at surrogate ripples
            IRISurro                        = datasample(IRI, numel(IRI), 'replace', false);
            rippleIdxSurro                  = cumsum(IRISurro);
            ripplePhasesSurro               = transpose(phase(rippleIdxSurro));
            
            % mean phase at ripples
            meanRipplePhaseSurro(iSurro, 1) = circ_mean(ripplePhasesSurro);
            
            % Rayleigh's test
            [~, RzSurro(iSurro, 1)]         = circ_rtest(ripplePhasesSurro);
        end
        
        %% differentiate between different power levels
        
        % power at ripples
        ripplePowers        = transpose(powerData(rippleIdx));
        
        % different power levels
        ripplePowerLevels   = nan(size(ripplePowers, 1), 1);
        for iLevel = 1:param.power.numLevels
            bThisLevel                      = ripplePowers >= quantile(ripplePowers, (iLevel - 1) / param.power.numLevels);
            ripplePowerLevels(bThisLevel)   = iLevel;
        end
        
        %% results
        
        % basic information
        unitRes                         = [];
        unitRes.idx                     = [iSub, iChan];
        unitRes.subjectName             = subjects{iSub};
        unitRes.chanName                = choi{iChan};
        
        % phase information
        unitRes.ripplePhases            = ripplePhases;
        unitRes.meanRipplePhase         = meanRipplePhase;
        unitRes.Rp                      = Rp;
        unitRes.Rz                      = Rz;
        
        % power information to examine phase locking as a function of power
        unitRes.ripplePowers            = ripplePowers;
        unitRes.ripplePowerLevels       = ripplePowerLevels;
        
        % surrogate phase information
        unitRes.RzSurro                 = RzSurro;
        unitRes.meanRipplePhaseSurro    = meanRipplePhaseSurro;
        
        % collapse across channels
        allRes  = cat(1, allRes, unitRes);
        
        %% figure: distribution of ripple phases
        
        % bin edges for the histogram
        xEdges = linspace(-pi, pi, 25);
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'visible', 'off');
        axes('units', 'centimeters', 'position', [1.6, 1.6, 3.5, 3.5]);
        hold on;
        histogram(ripplePhases, xEdges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
        xline(meanRipplePhase, '-', 'Color', [1, 0, 0], 'LineWidth', 2);
        xl = xlabel('Phase (rad)');
        yl = ylabel('Count');
        tl = title(['z = ', num2str(Rz, '%.2f'), ', \itP\rm = ', num2str(Rp, '%.2f')]);
        if Rp < 0.01
            tl = title(['z = ', num2str(Rz, '%.2f'), ', \itP\rm < 0.01']);
        end
        set(gca, ...
            'xlim', [min(xEdges), max(xEdges)], 'xtick', [-pi, 0, pi], 'xticklabel', {'-\pi', '0', '\pi'}, ...
            'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
        set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
        % save figure
        LK_print(f, strcat(paths.save, subjects{iSub}, '_', choi{iChan}, '_PhasesAtRipples_20210806'), '-dpng', '-r150');
        close(f);
        
        %% figure: filter response
        
        % estimate filter response
        [H, F]  = freqz(bpFilter, 8192, eeg.fsample);
        M       = abs(H); % get magnitude
        dB      = 20 * log10(M); % convert to dB
        
        % figure: filter response
        f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'visible', 'off');
        axes('units', 'centimeters', 'position', [1.6, 1.4, 4, 4]);
        hold on;
        plot(F, dB, 'Color', [1, 0, 0]);
        set(gca, 'xlim', [0, 10], 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
        tmpAx = get(gca);
        patch([min(param.phase.freqBand), max(param.phase.freqBand), max(param.phase.freqBand), min(param.phase.freqBand), min(param.phase.freqBand)], ...
            [min(tmpAx.YLim), min(tmpAx.YLim), max(tmpAx.YLim), max(tmpAx.YLim), min(tmpAx.YLim)], [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        grid minor;
        xl = xlabel('Frequency (Hz)');
        yl = ylabel('Magnitude (dB)', 'units', 'normalized', 'position', [-0.25, 0.5]);
        tl = title('Filter response');
        set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
        % save figure
        LK_print(f, strcat(paths.save, subjects{iSub}, '_', choi{iChan}, '_FilterResponse_20230405'), '-dpng', '-r300');
        close(f);
        
        %% figure: empirical data, filtered data, and ripple timepoints
        
        % ripple timepoints
        rippleTimepoints    = [ripples.ripples.peakTime]';
        
        % some example segments
        exStart = transpose(linspace(min(eeg.time{1}), max(eeg.time{1}), 51));
        exEnd   = exStart + 8; % seconds
        
        % loop through example segments
        for iEx = 1:(size(exStart, 1) - 1)
            
            % data from this segment
            bTwoi                   = eeg.time{1} >= exStart(iEx, 1) & eeg.time{1} <= exEnd(iEx, 1);
            thisRippleTimepoints    = rippleTimepoints(rippleTimepoints >= exStart(iEx, 1) & rippleTimepoints <= exEnd(iEx, 1));
            
            % create figure
            f = figure('units', 'centimeters', 'position', [2, 2, 16, 5], 'visible', 'off');
            axes('units', 'centimeters', 'position', [1.5, 1.5, 14, 3]);
            hold on;
            plot(eeg.time{1}(bTwoi), eeg.trial{1}(bTwoi) - mean(eeg.trial{1}(bTwoi)), 'Color', [0.5, 0.5, 0.5]);
            plot(eeg.time{1}(bTwoi), bpData(bTwoi), 'Color', [0, 0, 1, 0.5], 'LineWidth', 2);
            if ~isempty(thisRippleTimepoints)
                for iR = 1:numel(thisRippleTimepoints)
                    xline(thisRippleTimepoints(iR), '-', 'Color', [1, 0, 0], 'LineWidth', 1);
                end
            end
            xl = xlabel('Time (s)');
            yl = ylabel('Voltage (\muV)', 'units', 'normalized', 'position', [-0.018, 0.5]);
            set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
            tmpAx = get(gca);
            set(gca, 'xlim', [exStart(iEx, 1), exEnd(iEx, 1)], 'xtick', exStart(iEx, 1):exEnd(iEx, 1), 'xticklabel', {0:(exEnd(iEx, 1) - exStart(iEx, 1))}, ...
                'ylim', round([min(tmpAx.YLim), max(tmpAx.YLim)]), 'ytick', round([min(tmpAx.YLim), max(tmpAx.YLim)]), 'tickdir', 'out', 'ticklength', [0.01, 0.01]);
            title(sprintf('%d, %s, %d', iSub, choi{iChan}, iEx), 'fontunits', 'centimeters', 'fontsize', 0.2, 'fontweight', 'normal');
            % save figure
            LK_print(f, strcat(paths.save, subjects{iSub}, '_', choi{iChan}, '_Ripples2Delta_Ex', num2str(iEx), '_20230405'), '-dpng', '-r300');
            close(f);
        end
        
        %% close all open figures
        close all;
    end
end

%% save results
save(strcat(paths.save, 'results'));

%% load results
r   = load(strcat(paths.save, 'results.mat'));

%% distribution of all mean phases

% all ripple phases (one per ripple)
allRipplePhases     = cell2mat({r.allRes.ripplePhases}');

% all mean ripple phases (one per channel)
allMeanRipplePhase  = cell2mat({r.allRes.meanRipplePhase}');

% all raw Rayleigh p-values
allRp               = cell2mat({r.allRes.Rp}');

% report
fprintf('\nResults.\n');
fprintf('Number of all ripples with a preferred ripple phase: %d.\n', sum(~isnan(allRipplePhases)));
fprintf('Number of channels with a significant distribution (based on the raw Rayleigh p-value): %d (out of %d).\n', sum(allRp < 0.05), size(allRp, 1));

%% distribution of mean ripple phases across channels - statistics

% grand-average mean ripple phase, across channels
grandMeanRipplePhase        = circ_mean(allMeanRipplePhase);
grandSTDRipplePhase         = circ_std(allMeanRipplePhase);

% Rayleigh's test on the mean ripple phases, across channels
[Rp, Rz]                    = circ_rtest(allMeanRipplePhase);
fprintf('Grand-average ripple phase, across channels: %.3f +/- %.3f°.\n', rad2deg(grandMeanRipplePhase), rad2deg(grandSTDRipplePhase));
fprintf('Rayleigh test on the mean ripple phases, across channels: z = %.3f, p = %.3f.\n', Rz, Rp);

% Rayleigh's test on the surrogate mean ripple phases, across channels
allMeanRipplePhaseSurro     = cell2mat({r.allRes.meanRipplePhaseSurro})'; % e.g., size = 62 x 1001
grandMeanRipplePhaseSurro   = nan(param.numSurrogates, 1); % e.g., size = 1001 x 1
RpSurro                     = nan(param.numSurrogates, 1);
RzSurro                     = nan(param.numSurrogates, 1);
for iSurro = 1:r.param.numSurrogates
    
    % grand-average mean ripple phase, across channels
    grandMeanRipplePhaseSurro(iSurro, 1)        = circ_mean(allMeanRipplePhaseSurro(:, iSurro));
    
    % Rayleigh's test on the surrogate mean ripple phases, across channels
    [RpSurro(iSurro, 1), RzSurro(iSurro, 1)]    = circ_rtest(allMeanRipplePhaseSurro(:, iSurro));
end

% compare empirical Rayleigh's test with surrogate Rayleigh's tests
RzRank  = sum(Rz > RzSurro) / sum(~isnan(RzSurro));
fprintf('The Rayleigh z value on the empirical mean ripple phases is larger than %.2f%% of the surrogate Rayleigh z values (P = %.3f).\n', 100 * RzRank, 1 - RzRank);

% compare the empirical mean phases with the surrogate mean phases
[kuiperP, kuiperk] = circ_kuipertest(allMeanRipplePhase, allMeanRipplePhaseSurro(:));
fprintf('Grand-average surrogate ripple phase, across channels and surrogate rounds: %.3f°.\n', rad2deg(circ_mean(grandMeanRipplePhaseSurro(:))));
fprintf('Are "allMeanRipplePhase" different from "allMeanRipplePhaseSurro"? Kuiper''s k = %.3f, Kuiper''s p = %.3f.\n', kuiperk, kuiperP);

%% figure: distribution of mean ripple phases across channels

% create figure
dt                          = [];
dt.allMeanRipplePhase       = allMeanRipplePhase;
dt.allMeanRipplePhaseSurro  = allMeanRipplePhaseSurro;
dt.xEdges                   = linspace(-pi, pi, 21); % bin edges for the histogram
f = LK_PlotRipplePhaseLockingAcrossChannels_20231030(dt);

% save figure
LK_print(f, strcat(paths.save, 'allMeanRipplePhase_20210806'), '-dtiff', '-r300');

% save figure data
save(strcat(paths.save, 'Fig_2i'), 'dt');

%% phase locking as a function of power

% preallocate
RzLevel                 = nan(size(r.allRes, 1), r.param.power.numLevels);
allMeanRipplePhaseLevel = nan(size(r.allRes, 1), r.param.power.numLevels);

% loop through channels
for iChan = 1:size(r.allRes, 1)
    
    % loop through power levels
    for iLevel = min(r.allRes(iChan).ripplePowerLevels):max(r.allRes(iChan).ripplePowerLevels)
        
        % mean ripple phase for this power level
        thisRipplePhases                        = r.allRes(iChan).ripplePhases(r.allRes(iChan).ripplePowerLevels == iLevel);
        allMeanRipplePhaseLevel(iChan, iLevel)  = circ_mean(thisRipplePhases);
        
        % Rayleigh's test
        [~, RzLevel(iChan, iLevel)] = circ_rtest(thisRipplePhases);
    end
end

% stats: test for consistency of preferred phase across channels,
% separately for each power level
RpMeanPhasesLevel   = nan(r.param.power.numLevels, 1);
RzMeanPhasesLevel   = nan(r.param.power.numLevels, 1);
for iLevel = 1:r.param.power.numLevels
    [RpMeanPhasesLevel(iLevel, 1), RzMeanPhasesLevel(iLevel, 1)]    = circ_rtest(allMeanRipplePhaseLevel(:, iLevel));
end

% stats: repeated measures ANOVA to test whether the channel-wise phase
% locking is stronger for high as compared to medium and low power levels
t                       = table(RzLevel(:, 1), RzLevel(:, 2), RzLevel(:, 3), 'VariableNames', {'Level1', 'Level2', 'Level3'});
withinDesign            = table([1; 2; 3], 'VariableNames', {'PowerLevel'});
withinDesign.PowerLevel = categorical(withinDesign.PowerLevel);
rm                      = fitrm(t, 'Level1-Level3 ~ 1', 'WithinDesign', withinDesign);
ranovatbl               = ranova(rm, 'WithinModel', 'PowerLevel');

% post-hoc multiple comparisons
pPowerLevels            = nan(size(RzLevel, 2));
tPowerLevels            = nan(size(RzLevel, 2));
for i = 1:size(RzLevel, 2)
    for j = 1:size(RzLevel, 2)
        [~, p, ~, stats]    = ttest(RzLevel(:, i), RzLevel(:, j));
        pPowerLevels(i, j)  = p;
        tPowerLevels(i, j)  = stats.tstat;
    end
end
% Bonferroni correction for multiple comparisons
nComp           = size(RzLevel, 2) * (size(RzLevel, 2) - 1) / 2; % number of unique comparisons that can be made
pPowerLevels    = pPowerLevels * nComp;

% report
fprintf('\nComparison of phase-locking strength across different power levels (low, medium, high):\n');
fprintf('t-values:\n');
disp(tPowerLevels);
fprintf('p-values (Bonferroni-corrected)\n');
disp(pPowerLevels);

%% figure: distribution of preferred phases and consistency of preferred phases across channels

% bin edges for the histogram
xEdges  = linspace(-pi, pi, 21);

% loop through different power levels
for iLevel = 1:r.param.power.numLevels
    
    % report
    fprintf('Rayleigh test on the mean ripple phases, across channels, power level #%d: z = %.3f, p = %.3f.\n', iLevel, RzMeanPhasesLevel(iLevel, 1), RpMeanPhasesLevel(iLevel, 1));
    
    % adjust title
    if iLevel == 1
        titleText = 'Low delta power';
    elseif iLevel == 2
        titleText = 'Medium delta power';
    elseif iLevel == 3
        titleText = 'High delta power';
    end
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6], 'Color', [1, 1, 1]);
    axes('units', 'centimeters', 'position', [1.75, 1.4, 3.5, 3.5]);
    hold on;
    histogram(allMeanRipplePhaseLevel(:, iLevel), xEdges, 'normalization', 'probability', 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0, 0, 0]);
    xl = xlabel('Phase (rad)');
    yl = ylabel('Probability');
    t1 = text(0.05, 0.85, {['z = ', num2str(RzMeanPhasesLevel(iLevel, 1), '%.3f')], ['\itP\rm = ', num2str(RpMeanPhasesLevel(iLevel, 1), '%.3f')]}, 'units', 'normalized', ...
        'fontunits', 'centimeters', 'fontsize', 0.35);
    tl = title({titleText, ['(', LK_IndicateThousands(size(allMeanRipplePhase, 1)), ' channels)']});
    set(gca, ...
        'xlim', [min(xEdges), max(xEdges)], 'xtick', [-pi, 0, pi], 'xticklabel', {'-\pi', '0', '\pi'}, ...
        'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off', 'Layer', 'top');
    set([gca, xl, yl, tl], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    % save figure
    LK_print(f, strcat(paths.save, 'allMeanRipplePhasePowerLevel', num2str(iLevel), '_20230406'), '-dtiff', '-r300');
end

%% figure: phase-locking strength per channel per power level

% reset rng
rng(r.param.myRNG);

% mean and SEM
m   = mean(RzLevel, 1);
ste = LK_ste(RzLevel);

% create figure
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 4]);
hold on;
% plot Rayleigh z-values
for iLevel = 1:size(RzLevel, 2)
    bar(iLevel, m(iLevel), 'FaceColor', 'none', 'EdgeColor', [0, 0, 0]); % mean
    plot(iLevel + rand(size(RzLevel, 1), 1) * 0.5 - 0.25, RzLevel(:, iLevel), '.', 'Color', [0.7, 0.7, 0.7]); % individual data points
    plot([iLevel, iLevel], [m(iLevel) - ste(iLevel), m(iLevel) + ste(iLevel)], '-', 'Color', [0, 0, 0], 'LineWidth', 4); % SEM
end
% significance between group 1 and group 3
if pPowerLevels(1, 3) < 0.001
    plot([1, 1, 3, 3], [80, 100, 100, 80], '-', 'Color', [0, 0, 0]);
    text(2, 120, '***', 'Color', [0, 0, 0], 'fontunits', 'centimeters', 'fontsize', 0.4, 'horizontalalignment', 'center');
end
% significance between group 2 and group 3
if pPowerLevels(2, 3) < 0.001
    plot([2, 2, 3, 3], [220, 300, 300, 220], '-', 'Color', [0, 0, 0]);
    text(2.5, 380, '***', 'Color', [0, 0, 0], 'fontunits', 'centimeters', 'fontsize', 0.4, 'horizontalalignment', 'center');
end
xl = xlabel('Power');
yl = ylabel('Rayleigh''s \itz');
set([gca, xl, yl], 'FontUnits', 'centimeters', 'FontSize', 0.4);
set(gca, 'xlim', [0.4, 3.6], 'xtick', 1:size(RzLevel, 2), 'xticklabel', {'Low', 'Medium', 'High'}, 'xticklabelrotation', 0, ...
    'ylim', [0, 300], 'ytick', [1, 10, 100], 'yticklabel', {'1', '10', '100'}, 'yscale', 'log', 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
% save figure
LK_print(f, strcat(paths.save, 'phaseLockingByPower_20230406'), '-dpng', '-r300');
