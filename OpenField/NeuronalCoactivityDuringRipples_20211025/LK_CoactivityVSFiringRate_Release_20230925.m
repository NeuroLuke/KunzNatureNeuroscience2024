%==========================================================================
% This script examines the relationship between coactivity scores, firing
% rates, and brain regions.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; close all hidden; clc;
addpath(genpath('E:\OpenField\Functions\'));

% settings
param               = [];
param.conditions    = {'retPrefObjRLwPF'; 'retHalf2PrefObjRLwPF'; 'encPrefObjOLwPF'; 'encHalf2PrefObjOLwPF'};
param.coacName      = 'resCoac_z_';
param.coacStat      = 'permTestZero_z_'; % statistics to use for the temporal region of interest
% different groups
param.corrGroups    = {'meanFRObjCell', 'meanZTwoi'; 'meanFRPlaceCell', 'meanZTwoi'}; % correlation groups
param.corrType      = 'Pearson';
param.regions       = {'AMY'; 'EC'; 'HC'; 'PHC'; 'TP'; 'X'}; % X means all regions
param.cellGroups    = {'ObjectCells', 'Object cells'; 'PlaceCells', 'Place cells'}; % cell groups
param.FR            = {'Low'; 'High'};

% paths
paths               = [];
paths.coactivity    = 'E:\OpenField\NeuronalCoactivityDuringRipples_20211025\20230424_RipplesMacro_20210614_20220320_FBS_80to140Hz_tMF2\ripples_z_tR100_tS5\'; % coactivations
paths.units         = 'E:\OpenField\UnitQuality_20220324\Evaluation_20230330\'; % firing rates
paths.save          = 'E:\OpenField\NeuronalCoactivityDuringRipples_20211025\CoactivityVSFiringRates_20230522\';
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% coactivity settings from the coactivity analysis
coacParam           = load(strcat(paths.coactivity, 'settings.mat'));

% text file for storing the output
txtFile             = fopen(strcat(paths.save, 'statistics.txt'), 'w');

% all results
allRes              = [];

%% save settings
save(strcat(paths.save, 'settings'));

%% loop through conditions
for iCond = 1:size(param.conditions, 1)
    
    % report
    fprintf('\n\nCondition: %s.\n', param.conditions{iCond, 1});
    
    % load coactivity results, coactivity statistics, and firing rates
    resCoac     = load(strcat(paths.coactivity, param.coacName, param.conditions{iCond, 1}, '.mat')); % coactivity
    statsCoac 	= load(strcat(paths.coactivity, param.coacStat, param.conditions{iCond, 1}, '.mat')); % statistics
    resUnits    = load(strcat(paths.units, 'results.mat'), 'allRes');
    
    % define target time window
    logIdxTarget    = coacParam.param.coac.time >= min(coacParam.param.coac.target) & coacParam.param.coac.time <= max(coacParam.param.coac.target); % logical indexing of the target window
    tTarget         = coacParam.param.coac.time(logIdxTarget); % time during the target window
    
    % unfold data from this condition
    unfCoacZ2D          = cell2mat(resCoac.coacZ2D); % cell*cell*rippleChannel x twoiO x twoiP
    unfCoacObjCellIdx   = cell2mat(resCoac.coacObjCellIdx); % object-cell indices
    unfCoacPlaceCellIdx = cell2mat(resCoac.coacPlaceCellIdx); % place-cell indices
    
    % permute so that samples are the 3rd dimension
    unfCoacZ2D          = permute(unfCoacZ2D, [2, 3, 1]); % twoiO x twoiP x cell*cell*rippleChannel
    
    % restrict to target window
    coacTarget          = unfCoacZ2D(logIdxTarget, logIdxTarget, :);
    warning('The coacTarget has these dimensions: %d x %d x %d.', size(coacTarget, 1), size(coacTarget, 2), size(coacTarget, 3));
    
    % statistical results of this coactivity map to define time window of
    % interest (twoi) for extracting coactivity values
    bTwoi               = statsCoac.permTestZero.logIdxSigClus; % logIdxSigClus | logIdxMaxClus
    warning('The bTwoi has these dimensions: %d x %d.', size(bTwoi, 1), size(bTwoi, 2));
        
    % regions of object and place cells
    unfCoacObjCellRegion    = LK_unwrapNestedCell(resCoac.coacObjCellRegion);
    unfCoacPlaceCellRegion  = LK_unwrapNestedCell(resCoac.coacPlaceCellRegion);
    
    %% data: coactivity values and firing rates
    
    % preallocate
    meanZTwoi       = nan(size(unfCoacZ2D, 3), 1); % mean coactivity z-value in the twoi
    meanFRObjCell   = nan(size(unfCoacZ2D, 3), 1); % average firing rate of object cells
    meanFRPlaceCell = nan(size(unfCoacZ2D, 3), 1); % average firing rate of place cells
    
    % loop through coactivity maps
    for iCoac = 1:size(unfCoacZ2D, 3)
        
        % coactivity map for this cell combination
        thisCoacTarget      = coacTarget(:, :, iCoac);
                
        % average coactivity in the time window of interest (twoi)
        meanZTwoi(iCoac, 1) = mean(thisCoacTarget(bTwoi), 'omitnan');
                
        % mean firing rate of the object cell
        thisObjCellIdx              = unfCoacObjCellIdx(iCoac, :);
        logIdx                      = all(thisObjCellIdx == cell2mat({resUnits.allRes.idx}'), 2);
        meanFRObjCell(iCoac, 1)     = resUnits.allRes(logIdx).meanFR;
        
        % mean firing rate of the place cell
        thisPlaceCellIdx            = unfCoacPlaceCellIdx(iCoac, :);
        logIdx                      = all(thisPlaceCellIdx == cell2mat({resUnits.allRes.idx}'), 2);
        meanFRPlaceCell(iCoac, 1)   = resUnits.allRes(logIdx).meanFR;
    end
    
    %% correlations: coactivity values and firing rates
    
    % loop through the different groups
    for iG = 1:size(param.corrGroups, 1)
        
        % select data: cell type and type of coactivity
        if strcmp(param.corrGroups{iG, 1}, 'meanFRObjCell')
            x       = meanFRObjCell;
            myColor = rgb('orange');
        elseif strcmp(param.corrGroups{iG, 1}, 'meanFRPlaceCell')
            x       = meanFRPlaceCell;
            myColor = rgb('blue');
        end
        if strcmp(param.corrGroups{iG, 2}, 'meanZTwoi')
            y   = meanZTwoi;
        end
        
        % stats: correlation between firing rate and coactivity
        bValid      = ~isnan(x) & ~isnan(y);
        [rho, pval] = corr(x, y, 'rows', 'complete', 'type', param.corrType);
        fprintf('Correlation between "%s" and "%s": rho = %.3f, p = %.3f.\n', param.corrGroups{iG, 1}, param.corrGroups{iG, 2}, rho, pval);
        fprintf(txtFile, '%s, %s, %s, rho = %.6f, p = %.6f\n', param.conditions{iCond, 1}, param.corrGroups{iG, 1}, param.corrGroups{iG, 2}, rho, pval);
        
        % sanity check: test the coactivity values against 0
        if sum(~isnan(y)) >= 2
            [~, pY, ~, statsY]  = ttest(y);
            fprintf('\tAre the coactivity values above 0? t(%d) = %.3f, p = %.3f.\n', statsY.df, statsY.tstat, pY);
        end
        
        %% figure: correlation between firing rate and coactivity
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 6.5, 6.5]);
        
        % scatter plot with linear fit
        axes('units', 'centimeters', 'position', [1.3, 1.4, 4, 4]);
        hold on;
        plot(x, y, 'o', 'Color', myColor, 'MarkerSize', 3);
        if sum(bValid) >= 2
            fitobject = fit(x(bValid), y(bValid), 'poly1');
            plot(x, fitobject(x), 'LineWidth', 2, 'Color', [0, 0, 0]);
        end
        xl = xlabel('Firing rate (Hz)');
        yl = ylabel('Coactivity (z)', 'units', 'normalized', 'position', [-0.175, 0.5]);
        t1 = text(0.05, 0.025, sprintf('\\itr\\rm = %.2f, \\itP\\rm = %.3f', rho, pval), 'units', 'normalized', 'VerticalAlignment', 'bottom');
        set([gca, xl, yl, t1], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
        set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'on', 'layer', 'top');
        if sum(bValid) >= 2
            set(gca, 'xlim', [min(x(bValid)) - 1, max(x(bValid)) + 1], 'ylim', [min(y(bValid)) - 0.5, max(y(bValid)) + 0.5]); % adjust axis limits
        end
        ax1 = get(gca);
        
        % histogram: x
        axes('units', 'centimeters', 'position', [1.3, 5.4, 4, 1]);
        histogram(x(bValid), linspace(min(ax1.XLim), max(ax1.XLim), 11), 'FaceColor', myColor, 'EdgeColor', myColor, 'orientation', 'vertical');
        set(gca, 'xlim', ax1.XLim);
        axis off;
        
        % histogram: y
        axes('units', 'centimeters', 'position', [5.3, 1.4, 1, 4]);
        histogram(y(bValid), linspace(min(ax1.YLim), max(ax1.YLim), 11), 'FaceColor', myColor, 'EdgeColor', myColor, 'orientation', 'horizontal');
        set(gca, 'ylim', ax1.YLim);
        axis off;
        
        % save figure
        LK_print(f, strcat(paths.save, 'Corr_', param.conditions{iCond, 1}, '_', param.corrGroups{iG, 1}, '_', param.corrGroups{iG, 2}, '_20230925'), '-dpng', '-r600');
        
        %% inset: coactivity as a function of low vs. high firing rate
        
        % separate data
        i       = [1, 2];
        cutoff  = median(x, 'omitnan');
        d       = {y(x < cutoff), y(x >= cutoff)};
        
        % t-test
        if sum(~isnan(d{1})) >= 2 && sum(~isnan(d{2})) >= 2
            [~, pCoacFR]    = ttest2(d{1}, d{2});
        else
            pCoacFR = 1;
        end
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 3, 3]);
        axes('units', 'centimeters', 'position', [1.3, 1.2, 1.5, 1.4]);
        hold on;
        for iM = 1:size(d, 2)
            % report mean and SEM
            thisM   = mean(d{iM}, 1, 'omitnan');
            thisSTE = LK_ste(d{iM});
            fprintf('Mean and SEM: %.3f and %.3f.\n', thisM, thisSTE);
            % plot outliers
            bOutlier = isoutlier(d{iM}, 'quartiles'); % elements more than 1.5 IQRs above the upper quartile or below the lower quartile
            if any(bOutlier)
                plot(iM, d{iM}(bOutlier), '.', 'Color', [0, 0, 0]);
            end            
            % plot minimum and maximum
            thisMin = min(d{iM}(~bOutlier));
            thisMax = max(d{iM}(~bOutlier));
            plot(iM + [-0.1, 0.1, 0, 0, -0.1, 0.1], [thisMin, thisMin, thisMin, thisMax, thisMax, thisMax], '-', 'Color', [0, 0, 0]);
            % plot 25th and 75th percentile
            thisCI = prctile(d{iM}, [25, 75]);
            thisPatch = patch(iM + [-0.2, 0.2, 0.2, -0.2, -0.2], [thisCI(1), thisCI(1), thisCI(2), thisCI(2), thisCI(1)], myColor, 'EdgeColor', [0, 0, 0], 'LineWidth', 1);
            if iM == 2
                set(thisPatch, 'FaceColor', [1, 1, 1]);
            end
            % plot median
            thisMed = median(d{iM}, 'omitnan');
            plot(iM + [-0.2, 0.2], [thisMed, thisMed], '-', 'LineWidth', 1, 'Color', [0, 0, 0]);
        end
        if pCoacFR > 0.05
            t1 = text(0.5, 1.2, 'n.s.', 'units', 'normalized', 'HorizontalAlignment', 'center');
        else
            t1 = text(0.5, 1.1, '*', 'units', 'normalized', 'HorizontalAlignment', 'center');
        end
        tmpAx = get(gca);
        xl = xlabel('Firing rate', 'units', 'normalized', 'position', [0.5, -0.45]);
        yl = ylabel('z', 'units', 'normalized', 'position', [-0.2, 0.5], 'rotation', 0, 'VerticalAlignment', 'middle');
        set([gca, xl, yl, t1], 'FontUnits', 'centimeters', 'FontSize', 0.4);
        set(gca, 'xlim', [min(i) - 0.6, max(i) + 0.6], 'xtick', i, 'xticklabel', {'L', 'H'}, ...
            'ylim', tmpAx.YLim, 'ytick', tmpAx.YLim, 'yticklabel', {num2str(tmpAx.YLim', 2)}, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
        % save figure
        LK_print(f, strcat(paths.save, 'BarGraph_', param.conditions{iCond, 1}, '_', param.corrGroups{iG, 1}, '_', param.corrGroups{iG, 2}, '_20230925'), '-dpng', '-r300');
    end
    
    %% collect results
    
    % results from this condition
    thisRes                         = [];
    thisRes.condition               = param.conditions{iCond, 1};
    thisRes.meanZTwoi               = meanZTwoi;
    thisRes.meanFRObjCell           = meanFRObjCell;
    thisRes.meanFRPlaceCell         = meanFRPlaceCell;
    thisRes.unfCoacObjCellIdx       = unfCoacObjCellIdx;
    thisRes.unfCoacPlaceCellIdx     = unfCoacPlaceCellIdx;
    thisRes.unfCoacObjCellRegion    = unfCoacObjCellRegion;
    thisRes.unfCoacPlaceCellRegion  = unfCoacPlaceCellRegion;
    
    % collect across conditions
    allRes  = cat(1, allRes, thisRes);
    
    % close figures
    close all;
end

%% comparison of coactivity during retrieval and re-encoding

% comparisons
comparisons     = {...
    'retPrefObjRLwPF', 'encPrefObjOLwPF', 'meanZTwoi', 'Retrieval', 'Re-encoding'; ... % retrieval vs. re-encoding, all
    'retHalf2PrefObjRLwPF', 'encHalf2PrefObjOLwPF', 'meanZTwoi', 'Retrieval', 'Re-encoding'; ... % retrieval vs. re-encoding, 2nd half
    };

% loop through comparisons
for iComp = 1:size(comparisons, 1)
    
    % report
    fprintf('\n\nComparison: "%s" vs. "%s" for "%s".\n', comparisons{iComp, 1}, comparisons{iComp, 2}, comparisons{iComp, 3});    
    
    % logical indices to the conditions
    idx1    = strcmp({allRes.condition}', comparisons{iComp, 1}); % first condition
    idx2    = strcmp({allRes.condition}', comparisons{iComp, 2}); % second condition
    
    % z-values of the conditions
    if strcmp(comparisons{iComp, 3}, 'meanZTwoi')
        z1  = allRes(idx1).meanZTwoi; % e.g., size = 1104 x 1
        z2 	= allRes(idx2).meanZTwoi;
    else
        error('Type of coactivity has not been specified.');
    end
    
    %% relationship between coactivations during retrieval and during re-encoding
    
    % correlation between both conditions
    bValid      = ~isnan(z1) & ~isnan(z2);
    if sum(bValid) <= 1
        warning('Not enough data for correlating the coactivity values between the two conditions, thus skipping.');
        continue;
    end
    [rho, pval] = corr(z1(bValid), z2(bValid), 'type', param.corrType);
    fprintf('Correlation between "%s" and "%s" for "%s": rho = %.3f, p = %.3f.\n', comparisons{iComp, 1}, comparisons{iComp, 2}, comparisons{iComp, 3}, rho, pval);
    fprintf(txtFile, '%s, %s, %s, rho = %.6f, p = %.6f\n', comparisons{iComp, 1}, comparisons{iComp, 2}, comparisons{iComp, 3}, rho, pval);
        
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 7, 7]);
    % scatter plot
    axes('units', 'centimeters', 'position', [1.9, 1.9, 4, 4]);
    hold on;
    plot(z1, z2, 'o', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 3);
    fitobject = fit(z1(bValid), z2(bValid), 'poly1');
    p = plot([-10, 10], fitobject([-10, 10]), '-', 'Color', [0, 0, 0], 'LineWidth', 1.5);
    xl = xlabel({'Coactivity (z)', comparisons{iComp, 4}});
    yl = ylabel({comparisons{iComp, 5}, 'Coactivity (z)'});
    t1 = text(0.05, 0.025, sprintf('\\itr\\rm = %.2f, \\itP\\rm = %.3f', rho, pval), 'units', 'normalized', 'VerticalAlignment', 'bottom');
    set([gca, xl, yl, t1], 'fontunits', 'centimeters', 'fontsize', 0.4, 'fontweight', 'normal');
    set(gca, 'xlim', [min(z1(bValid)) - 0.5, max(z1(bValid)) + 0.5], 'ylim', [min(z2(bValid)) - 0.5, max(z2(bValid)) + 0.5], 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'on');
    ax1 = get(gca);
    % histogram: z1
    axes('units', 'centimeters', 'position', [1.9, 5.9, 4, 1]);
    histogram(z1(bValid), linspace(min(ax1.XLim), max(ax1.XLim), 11), 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', [0.7, 0.7, 0.7], 'orientation', 'vertical');
    set(gca, 'xlim', ax1.XLim);
    axis off;
    % histogram: z2
    axes('units', 'centimeters', 'position', [5.9, 1.9, 1, 4]);
    histogram(z2(bValid), linspace(min(ax1.YLim), max(ax1.YLim), 11), 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', [0.7, 0.7, 0.7], 'orientation', 'horizontal');
    set(gca, 'ylim', ax1.YLim);
    axis off;
    % save figure
    LK_print(f, strcat(paths.save, 'Corr_', comparisons{iComp, 1}, '_', comparisons{iComp, 2}, '_', comparisons{iComp, 3}, '_20230925'), '-dpng', '-r300');
    
    %% compare coactivity values between retrieval and re-encoding as a function of brain region
    
    % separately for place and object cells
    for iCG = 1:size(param.cellGroups, 1)
        
        %% stats
        
        % brain regions
        if strcmp(param.cellGroups{iCG, 1}, 'ObjectCells')
            reg1    = allRes(idx1).unfCoacObjCellRegion;
            reg2    = allRes(idx2).unfCoacObjCellRegion;
        elseif strcmp(param.cellGroups{iCG, 1}, 'PlaceCells')
            reg1    = allRes(idx1).unfCoacPlaceCellRegion;
            reg2    = allRes(idx2).unfCoacPlaceCellRegion;
        end
        
        % dependent and independent variables
        DV  = [z1; z2]; % coactivity scores
        IV1 = categorical([ones(size(z1)); 2 * ones(size(z2))]); % trial phase
        IV2 = [reg1; reg2]; % brain region
                
        % restrict the analysis to regions with enough data
        bNotNan             = ~isnan(DV);
        [N, ~, ~, labels]   = crosstab(IV1(bNotNan), IV2(bNotNan));
        validRegions        = labels(all(N > 0, 1), 2); % regions with at least one observation per trial phase
        
        % restrict the analysis to main and valid regions
        bMainReg            = ismember(IV2, param.regions) & ismember(IV2, validRegions);
        
        % two-way ANOVA to estimate interaction
        [pCoacPhaseReg, tblCoacPhaseReg]    = anovan(DV(bMainReg), {IV1(bMainReg), IV2(bMainReg)}, 'model', 'full', 'varnames', {'Phase', 'Region'});
        pInterCoacPhaseReg                  = tblCoacPhaseReg{strcmp(tblCoacPhaseReg(:, 1), 'Phase*Region'), contains(tblCoacPhaseReg(1, :), 'Prob')};
        
        % report ANOVA results
        fprintf('\n%s. Two-way ANOVA with factors "Phase" and "Region" on "%s":\n', param.cellGroups{iCG, 2}, comparisons{iComp, 3});
        disp(tblCoacPhaseReg);
        
        %% figure
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
        axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 4]);
        hold on;
        for iReg = 1:size(param.regions, 1) % loop through regions
            
            % select cells from this region
            if strcmp(param.cellGroups{iCG, 1}, 'ObjectCells')
                if strcmp(param.regions{iReg}, 'X')
                    bReg1   = true(size(allRes(idx1).unfCoacObjCellRegion));
                    bReg2   = true(size(allRes(idx2).unfCoacObjCellRegion));
                else 
                    bReg1   = strcmp(allRes(idx1).unfCoacObjCellRegion, param.regions{iReg});
                    bReg2   = strcmp(allRes(idx2).unfCoacObjCellRegion, param.regions{iReg});
                end
            elseif strcmp(param.cellGroups{iCG, 1}, 'PlaceCells')
                if strcmp(param.regions{iReg}, 'X')
                    bReg1   = true(size(allRes(idx1).unfCoacPlaceCellRegion));
                    bReg2   = true(size(allRes(idx2).unfCoacPlaceCellRegion));
                else
                    bReg1   = strcmp(allRes(idx1).unfCoacPlaceCellRegion, param.regions{iReg});
                    bReg2   = strcmp(allRes(idx2).unfCoacPlaceCellRegion, param.regions{iReg});
                end
            end
            
            % coactivity scores from the first condition
            d1          = z1(bReg1);
            m1          = mean(d1, 1, 'omitnan');
            ste1        = LK_ste(d1);
            med1        = median(d1, 1, 'omitnan');
            ci1         = prctile(d1, [25, 75]);
            bOutlier1   = isoutlier(d1, 'quartiles');
            min1        = min(d1(~bOutlier1));
            max1        = max(d1(~bOutlier1));

            % coactivity scores from the second condition
            d2          = z2(bReg2);
            m2   	    = mean(d2, 1, 'omitnan');
            ste2        = LK_ste(d2);
            med2        = median(d2, 1, 'omitnan');
            ci2         = prctile(d2, [25, 75]);
            bOutlier2   = isoutlier(d2, 'quartiles');
            min2        = min(d2(~bOutlier2));
            max2        = max(d2(~bOutlier2));
            
            % stats: comparison against 0
            if sum(~isnan(z1(bReg1))) >= 2
                [~, p1] = ttest(z1(bReg1)); % first condition
            else
                p1 = nan;
            end
            if sum(~isnan(z2(bReg2))) >= 2
                [~, p2] = ttest(z2(bReg2)); % second condition
            else
                p2 = nan;
            end
            
            % stats: comparison between conditions
            if sum(~isnan(z1(bReg1))) >= 2 && sum(~isnan(z2(bReg2))) >= 2
                [~, p1VS2, ~, stats1VS2]    = ttest2(z1(bReg1), z2(bReg2));
            else
                p1VS2 = nan;
                stats1VS2 = nan;
            end
            
            % boxplot 1
            if any(bOutlier1)
                plot(iReg - 0.2, d1(bOutlier1), '.', 'Color', rgb('orange')); % outliers
            end
            plot(iReg - 0.2 + 0.75 * [-0.1, 0.1, 0, 0, -0.1, 0.1], [min1, min1, min1, max1, max1, max1], 'Color', [0, 0, 0], 'LineWidth', 1); % min and max
            patch1 = patch(iReg - 0.2 + 0.75 * [-0.2, 0.2, 0.2, -0.2, -0.2], [ci1(1), ci1(1), ci1(2), ci1(2), ci1(1)], rgb('orange')); % quartiles
            plot1 = plot(iReg - 0.2 + 0.75 * [-0.2, 0.2], [med1, med1], '-', 'Color', [0, 0, 0], 'LineWidth', 2);
            % boxplot 2
            if any(bOutlier2)
                plot(iReg + 0.2, d2(bOutlier2), '.', 'Color', rgb('darkgreen')); % outliers
            end
            plot(iReg + 0.2 + 0.75 * [-0.1, 0.1, 0, 0, -0.1, 0.1], [min2, min2, min2, max2, max2, max2], 'Color', [0, 0, 0], 'LineWidth', 1); % min and max
            patch2 = patch(iReg + 0.2 + 0.75 * [-0.2, 0.2, 0.2, -0.2, -0.2], [ci2(1), ci2(1), ci2(2), ci2(2), ci2(1)], rgb('darkgreen')); % quartiles
            plot2 = plot(iReg + 0.2 + 0.75 * [-0.2, 0.2], [med2, med2], '-', 'Color', [0, 0, 0], 'LineWidth', 2);

            % significance
            if p1 < 0.05
                text(iReg - 0.2, max1, '*', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'baseline', 'FontSize', 12);
            end
            if p2 < 0.05
                text(iReg + 0.2, max2, '*', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'baseline', 'FontSize', 12);
            end
            if p1VS2 < 0.05
                set([plot1, plot2], 'Color', [1, 0, 0]);
            end
        end
        % enhance axes
        xl = xlabel('Region');
        yl = ylabel('Coactivity (z)', 'units', 'normalized', 'position', [-0.25, 0.5]);
        tl = title([param.cellGroups{iCG, 2}, sprintf(' (\\itP\\rm = %.3f)', pInterCoacPhaseReg)], 'units', 'normalized', 'position', [0.5, 1.02]);
        tmpAx = get(gca);
        set(gca, 'xlim', [0.4, numel(param.regions) + 0.6], 'xtick', 1:numel(param.regions), 'xticklabel', cellfun(@(x) x(1), param.regions), ...
            'ylim', tmpAx.YLim .* 1.1, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
        set([gca, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
        % save figure
        LK_print(f, strcat(paths.save, 'RegionWise_', comparisons{iComp, 1}, '_', comparisons{iComp, 2}, '_', comparisons{iComp, 3}, '_', param.cellGroups{iCG, 1}, '_20230925'), '-dpng', '-r300');
    end
    
    %% compare coactivity values between retrieval and re-encoding as a function of low vs. high firing rates (median split)
    
    % separately for place and object cells
    for iCG = 1:size(param.cellGroups, 1)
        
        %% stats
        
        % firing-rate groups
        if strcmp(param.cellGroups{iCG, 1}, 'ObjectCells')
            FR1 = (allRes(idx1).meanFRObjCell < median(allRes(idx1).meanFRObjCell, 'omitnan')) + 1; % 1 = low FR; 2 = high FR
            FR2 = (allRes(idx2).meanFRObjCell < median(allRes(idx2).meanFRObjCell, 'omitnan')) + 1;
        elseif strcmp(param.cellGroups{iCG, 1}, 'PlaceCells')
            FR1 = (allRes(idx1).meanFRPlaceCell < median(allRes(idx1).meanFRPlaceCell, 'omitnan')) + 1;
            FR2 = (allRes(idx2).meanFRPlaceCell < median(allRes(idx2).meanFRPlaceCell, 'omitnan')) + 1;
        end
        
        % dependent and independent variables
        DV  = [z1; z2];
        IV1 = categorical([ones(size(z1)); 2 * ones(size(z2))]); % e.g., 1 = retrieval, 2 = re-encoding
        IV2 = categorical([FR1; FR2]); % high and low firing rates
        
        % two-way ANOVA
        [pCoacPhaseFR, tblCoacPhaseFR]  = anovan(DV, {IV1, IV2}, 'model', 'full', 'varnames', {'Phase', 'FR'});
        pInterCoacPhaseFR               = tblCoacPhaseFR{strcmp(tblCoacPhaseFR(:, 1), 'Phase*FR'), contains(tblCoacPhaseFR(1, :), 'Prob')};
        
        % report ANOVA results
        fprintf('\n%s. Two-way ANOVA with factors "Phase" and "FR" on "%s":\n', param.cellGroups{iCG, 2}, comparisons{iComp, 3});
        disp(tblCoacPhaseFR);
        
        %% figure
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
        axes('units', 'centimeters', 'position', [1.7, 1.4, 4, 4]);
        hold on;
        for iFR = min(FR1):max(FR1) % loop through different firing-rate groups
                        
            % coactivity scores from the first condition
            d1          = z1(FR1 == iFR);
            m1          = mean(d1, 1, 'omitnan');
            ste1        = LK_ste(d1);
            med1        = median(d1, 1, 'omitnan');
            ci1         = prctile(d1, [25, 75]);
            bOutlier1   = isoutlier(d1, 'quartiles');
            min1        = min(d1(~bOutlier1));
            max1        = max(d1(~bOutlier1));
            
            % coactivity scores from the second condition
            d2          = z2(FR2 == iFR);
            m2   	    = mean(d2, 1, 'omitnan');
            ste2        = LK_ste(d2);
            med2        = median(d2, 1, 'omitnan');
            ci2         = prctile(d2, [25, 75]);
            bOutlier2   = isoutlier(d2, 'quartiles');
            min2        = min(d2(~bOutlier2));
            max2        = max(d2(~bOutlier2));
            
            % stats: comparison against 0
            if sum(~isnan(z1(FR1 == iFR))) >= 2
                [~, p1] = ttest(z1(FR1 == iFR));
            else
                p1 = nan;
            end
            if sum(~isnan(z2(FR2 == iFR))) >= 2
                [~, p2] = ttest(z2(FR2 == iFR));
            else
                p2 = nan;
            end
            
            % stats: comparison between conditions
            if sum(~isnan(z1(FR1 == iFR))) >= 2 && sum(~isnan(z2(FR2 == iFR))) >= 2
                [~, p1VS2, ~, stats1VS2]    = ttest2(z1(FR1 == iFR), z2(FR2 == iFR));
            else
                p1VS2       = nan;
                stats1VS2   = nan;
            end
            
            % boxplot 1
            if any(bOutlier1)
                plot(iFR - 0.2, d1(bOutlier1), '.', 'Color', rgb('orange')); % outliers
            end
            plot(iFR - 0.2 + 0.75 * [-0.1, 0.1, 0, 0, -0.1, 0.1], [min1, min1, min1, max1, max1, max1], 'Color', [0, 0, 0], 'LineWidth', 1); % min and max
            patch1 = patch(iFR - 0.2 + 0.75 * [-0.2, 0.2, 0.2, -0.2, -0.2], [ci1(1), ci1(1), ci1(2), ci1(2), ci1(1)], rgb('orange')); % quartiles
            plot1 = plot(iFR - 0.2 + 0.75 * [-0.2, 0.2], [med1, med1], '-', 'Color', [0, 0, 0], 'LineWidth', 2);
            % boxplot 2
            if any(bOutlier2)
                plot(iFR + 0.2, d2(bOutlier2), '.', 'Color', rgb('darkgreen')); % outliers
            end
            plot(iFR + 0.2 + 0.75 * [-0.1, 0.1, 0, 0, -0.1, 0.1], [min2, min2, min2, max2, max2, max2], 'Color', [0, 0, 0], 'LineWidth', 1); % min and max
            patch2 = patch(iFR + 0.2 + 0.75 * [-0.2, 0.2, 0.2, -0.2, -0.2], [ci2(1), ci2(1), ci2(2), ci2(2), ci2(1)], rgb('darkgreen')); % quartiles
            plot2 = plot(iFR + 0.2 + 0.75 * [-0.2, 0.2], [med2, med2], '-', 'Color', [0, 0, 0], 'LineWidth', 2);
            
            % significance
            if p1 < 0.05
                text(iFR - 0.2, max1, '*', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'baseline', 'FontSize', 12);
            end
            if p2 < 0.05
                text(iFR + 0.2, max2, '*', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'baseline', 'FontSize', 12);
            end
            if p1VS2 < 0.05
                set([plot1, plot2], 'Color', [1, 0, 0]); % indicate significance with red medians
            end
        end
        % enhance axes
        xl = xlabel('Firing rate');
        yl = ylabel('Coactivity (z)', 'units', 'normalized', 'position', [-0.25, 0.5]);
        tl = title([param.cellGroups{iCG, 2}, sprintf(' (\\itP\\rm = %.3f)', pInterCoacPhaseFR)], 'units', 'normalized', 'position', [0.5, 1.02]);
        tmpAx = get(gca);
        set(gca, 'xlim', [0.4, numel(param.FR) + 0.6], 'xtick', 1:numel(param.FR), 'xticklabel', param.FR, 'ylim', tmpAx.YLim .* 1.1, 'tickdir', 'out', 'ticklength', [0.02, 0.02]);
        set([gca, xl, yl, tl], 'FontUnits', 'centimeters', 'FontSize', 0.4, 'FontWeight', 'normal');
        % save figure
        LK_print(f, strcat(paths.save, 'FRWise_', comparisons{iComp, 1}, '_', comparisons{iComp, 2}, '_', comparisons{iComp, 3}, '_', param.cellGroups{iCG, 1}, '_20230925'), '-dpng', '-r300');
    end
end

% close txt file
fclose all;