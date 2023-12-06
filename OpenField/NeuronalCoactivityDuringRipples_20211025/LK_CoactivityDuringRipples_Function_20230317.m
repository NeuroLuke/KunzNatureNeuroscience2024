function out = LK_CoactivityDuringRipples_Function_20230317(cfg)
%
% LK_CoactivityDuringRipples_Function_20230317 computes the coactivity
% between object cells and place cells during ripples. The coactivity is
% tested for different ripple conditions.
%
% The coactivity is estimated via a z-score (Sosa et al., Neuron, 2020) or
% via a Pearson correlation.
%
% Input and output are structures with multiple fields.
%
% Lukas Kunz, 2023

% report
fprintf('\n========================================================================\n');
fprintf('Coactivity analysis between object cells and place cells during ripples.\n');

% output folders
mkdir(cfg.r.paths.save);
mkdir(cfg.r.paths.pics);

%% ripple conditions

% different ripple conditions for calculating the coactivities
conditions  = {...
    'all'; ... % all ripples (for sanity checks)
    
    'cueRet'; ... % ripples during cue & retrieval
    'cueRetPrefObjRLwPF'; 'cueRetHalf1PrefObjRLwPF'; 'cueRetHalf2PrefObjRLwPF'; ... % preferred object, response location within place field
    'cueRetPrefObjRLoPF'; 'cueRetHalf1PrefObjRLoPF'; 'cueRetHalf2PrefObjRLoPF'; ... % preferred object, response location outside place field
    
    'ret'; ... % ripples during retrieval
    'retPrefObjRLwPF'; 'retHalf1PrefObjRLwPF'; 'retHalf2PrefObjRLwPF'; ... % preferred object, response location within place field
    'retPrefObjRLoPF'; 'retHalf1PrefObjRLoPF'; 'retHalf2PrefObjRLoPF'; ... % preferred object, response location outside place field
    'retMovPrefObjRLwPF'; 'retMovHalf1PrefObjRLwPF'; 'retMovHalf2PrefObjRLwPF'; ... % movement, preferred object, response location within place field
    'retMovPrefObjRLoPF'; 'retMovHalf1PrefObjRLoPF'; 'retMovHalf2PrefObjRLoPF'; ... % movement, preferred object, response location outside place field
    'retNoMovPrefObjRLwPF'; 'retNoMovHalf1PrefObjRLwPF'; 'retNoMovHalf2PrefObjRLwPF'; ... % no-movement, preferred object, response location within place field
    'retNoMovPrefObjRLoPF'; 'retNoMovHalf1PrefObjRLoPF'; 'retNoMovHalf2PrefObjRLoPF'; ... % no-movement, preferred object, response location outside place field
    'retBeLePrefObjRLwPF'; 'retAfLePrefObjRLwPF'; ... % learning, preferred object, response location within place field
    'retBeLePrefObjRLoPF'; 'retAfLePrefObjRLoPF'; ... % learning, preferred object, response location outside place field
    
    'enc'; ... % ripples during (re)encoding
    'encPrefObjOLwPF'; 'encHalf1PrefObjOLwPF'; 'encHalf2PrefObjOLwPF'; ... % preferred object, object location within place field
    'encPrefObjOLoPF'; 'encHalf1PrefObjOLoPF'; 'encHalf2PrefObjOLoPF'; ... % preferred object, object location outside place field
    'encMovPrefObjOLwPF'; 'encMovHalf1PrefObjOLwPF'; 'encMovHalf2PrefObjOLwPF'; ... % movement, preferred object, object location within place field
    'encMovPrefObjOLoPF'; 'encMovHalf1PrefObjOLoPF'; 'encMovHalf2PrefObjOLoPF'; ... % movement, preferred object, object location outside place field
    'encNoMovPrefObjOLwPF'; 'encNoMovHalf1PrefObjOLwPF'; 'encNoMovHalf2PrefObjOLwPF'; ... % no-movement, preferred object, object location within place field
    'encNoMovPrefObjOLoPF'; 'encNoMovHalf1PrefObjOLoPF'; 'encNoMovHalf2PrefObjOLoPF'; ... % no-movement, preferred object, object location outside place field
    'encBeLePrefObjOLwPF'; 'encAfLePrefObjOLwPF'; ... % learning, preferred object, object location within place field
    'encBeLePrefObjOLoPF'; 'encAfLePrefObjOLoPF'; ... % learning, preferred object, object location outside place field
    
    'encNITI'; ... % ripples during (re)encoding and next ITI
    'encNITIPrefObjOLwPF'; 'encNITIHalf1PrefObjOLwPF'; 'encNITIHalf2PrefObjOLwPF'; ... % preferred object, object location within place field
    'encNITIPrefObjOLoPF'; 'encNITIHalf1PrefObjOLoPF'; 'encNITIHalf2PrefObjOLoPF'; ... % preferred object, object location outside place field
    };

%% coactivity of object cells and place cells during ripples from different conditions

% loop through the different ripple conditions
for iCond = 1:numel(conditions)
    
    % reset rng for reproducibility
    rng(numel(conditions{iCond}));
    
    % report condition
    fprintf('\n\n========================================================== CONDITION #%d: %s.\n', iCond, conditions{iCond});
    
    % name of the file for saving
    saveFile = strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', conditions{iCond}, '.mat');
    if exist(saveFile, 'file') > 0
        fprintf('\tResults were previously saved, thus skipping.\n');
        continue;
    end
    
    % coactivity z-scores
    coacZ2D                	= cell(size(cfg.r.allRes, 1), 1); % 2D coactivity
    coacZ2DSurro            = cell(size(cfg.r.allRes, 1), 1); % surrogate 2D coactivity
    coacNumRipples          = cell(size(cfg.r.allRes, 1), 1); % number of ripples contributing to the coactivity values
    coacNumUniqueTrials     = cell(size(cfg.r.allRes, 1), 1); % number of unique trials contributing to the coactivity values
    
    % cell indices and brain regions relevant to the coactivity analysis
    coacObjCellIdx        	= cell(size(cfg.r.allRes, 1), 1); % indices of the object cells
    coacPlaceCellIdx       	= cell(size(cfg.r.allRes, 1), 1); % indices of the place cells
    coacObjCellRegion     	= cell(size(cfg.r.allRes, 1), 1); % brain region of the object cells
    coacPlaceCellRegion    	= cell(size(cfg.r.allRes, 1), 1); % brain region of the place cells
    coacRippleChanName     	= cell(size(cfg.r.allRes, 1), 1); % ripple channels
    coacBSameHemiObjRipple  = cell(size(cfg.r.allRes, 1), 1); % whether the object cell was recorded from the same hemisphere as the ripples
    coacBSameHemiObjPlace   = cell(size(cfg.r.allRes, 1), 1); % whether the object cell was recorded from the same hemisphere as the place cell
    
    % loop through all cells
    for iRes = 1:size(cfg.r.allRes, 1)
        
        % report status
        if mod(iRes, 200) == 1
            fprintf('\tWorking on object cell-ripple channel-combination #%d.\n', iRes);
        end
        
        % skip if this is not an object cell
        if ~strcmp(cfg.allCellTypes(iRes).objectCell, 'ObjectCell')
            continue;
        end
        
        %% cellular indexing
        
        % indexing of this object cell
        thisObjCellIdx      = cfg.allIdx(iRes, :);
        bThisObjCell        = all(cfg.allIdx == thisObjCellIdx, 2) & strcmp(cfg.allRippleChanName, cfg.allRippleChanName{iRes});
        
        % place cells that were simultaneously recorded with this object
        % cell in this subject relative to this ripple channel
        bThisPlaceCells = strcmp({cfg.allCellTypes.placeCell}', 'PlaceCell') & ...
            cfg.allIdx(:, 1) == cfg.allIdx(iRes, 1) & strcmp(cfg.allRippleChanName, cfg.allRippleChanName{iRes}) & ~bThisObjCell;
        
        % indexing of these place cells
        thisPlaceCellsIdx   = cfg.allIdx(bThisPlaceCells, :);
        
        % skip if there are no simultaneously recorded place cells
        if size(thisPlaceCellsIdx, 1) < 1
            fprintf('\tSkipping object cell [%d, %d, %d] because there are no co-recorded place cells.\n', thisObjCellIdx);
            continue;
        end
        
        % cellular indices
        coacObjCellIdx{iRes, 1}         = repmat(thisObjCellIdx, size(thisPlaceCellsIdx, 1), 1);
        coacPlaceCellIdx{iRes, 1}       = thisPlaceCellsIdx;
        
        % brain regions
        coacObjCellRegion{iRes, 1}      = repmat(cfg.allSpikeChanRegion(iRes), size(thisPlaceCellsIdx, 1), 1); % brain region of this object cell
        coacPlaceCellRegion{iRes, 1}    = cfg.allSpikeChanRegion(bThisPlaceCells); % brain regions of the simultaneously recorded place cells
        
        % info whether the object cell and the ripple were recorded from
        % the same hemisphere
        coacRippleChanName{iRes, 1}     = repmat(cfg.allRippleChanName(iRes), size(thisPlaceCellsIdx, 1), 1); % name of the ripple channel
        coacBSameHemiObjRipple{iRes, 1} = double(repmat(cfg.allBSameHemiSpikeRipple(iRes), size(thisPlaceCellsIdx, 1), 1));
        
        % whether the object cell and the place cell were recorded from the
        % same hemisphere
        coacBSameHemiObjPlace{iRes, 1}  = double(strcmp(cfg.allSpikeChanHemisphere(iRes), cfg.allSpikeChanHemisphere(bThisPlaceCells)));
        
        %% unpack ripple information
        
        % trial phase and trial index per ripple
        trialPhasePerRipple = cfg.r.allRes(iRes).trialPhasePerRipple;
        trialIdxPerRipple   = cfg.r.allRes(iRes).trialIdxPerRipple;
        
        % ripples during trials with the preferred object of the object
        % cell
        tmpObjIdx           = all(cell2mat({cfg.r.objRes.allRes.idx}') == thisObjCellIdx, 2);
        prefObjName         = cfg.r.objRes.allRes(tmpObjIdx).prefObjName; % preferred object of this object cell (ranges between 0 and 7)
        bPrefObjRipples     = cfg.r.allRes(iRes).objNamePerRipple == prefObjName; % ripples from trials with the preferred object of the object cell
                
        % first half of all ripples
        bHalf1PerRipple     = cfg.r.allRes(iRes).bHalf1PerRipple;
        
        % ripples before vs. after learning
        bB4LearnPerRipple   = cfg.r.allRes(iRes).bB4LearnPerRipple;
        
        %% loop through place cells and estimate coactivity
        
        % preallocate
        zScore2D        = nan(size(thisPlaceCellsIdx, 1), size(cfg.r.param.coac.twoi, 1), size(cfg.r.param.coac.twoi, 1)); % coactivity z-scores (place-cells x timepoints x timepoints)
        zScore2DSurro   = nan(size(zScore2D)); % surrogate coactivity z-scores
        numRipples      = nan(size(zScore2D)); % number of ripples contributing to the coactivity values
        numUniqueTrials = nan(size(zScore2D)); % number of unique trials contributing to the coactivity values
        
        % loop through place cells
        for iPlaceCell = 1:size(thisPlaceCellsIdx, 1)
            
            %% indexing and place field of this place cell
            
            % indexing of this place cell
            thisPlaceCellIdx    = thisPlaceCellsIdx(iPlaceCell, :);
            bThisPlaceCell      = all(cfg.allIdx == thisPlaceCellIdx, 2) & strcmp(cfg.allRippleChanName, cfg.allRippleChanName{iRes});
            
            % place field of this place cell
            tmpPlaceIdx         = all(cell2mat({cfg.r.placeRes.allRes.idx}') == thisPlaceCellIdx, 2);
            thisPF              = cfg.r.placeRes.allRes(tmpPlaceIdx).PF;
            thisPFIdx           = cfg.r.placeRes.place.idxTemplate(thisPF);
            
            %% ripples when the object location / response location is within the place field of the place cell
            
            % ripples during trials in which the object location is within
            % the place field
            tmpXBins        = discretize(cfg.r.allRes(iRes).objLocPerRipple(:, 1), cfg.r.placeRes.place.xEdges);
            tmpYBins        = discretize(cfg.r.allRes(iRes).objLocPerRipple(:, 2), cfg.r.placeRes.place.yEdges);
            tmpXYBins       = (tmpXBins - 1) .* numel(cfg.r.placeRes.place.yCenters) + tmpYBins;
            bOLwPFRipples   = ismember(tmpXYBins, thisPFIdx);
            
            % ripples during trials in which the response location is
            % within the place field
            tmpXBins        = discretize(cfg.r.allRes(iRes).respLocPerRipple(:, 1), cfg.r.placeRes.place.xEdges);
            tmpYBins        = discretize(cfg.r.allRes(iRes).respLocPerRipple(:, 2), cfg.r.placeRes.place.yEdges);
            tmpXYBins       = (tmpXBins - 1) .* numel(cfg.r.placeRes.place.yCenters) + tmpYBins;
            bRLwPFRipples   = ismember(tmpXYBins, thisPFIdx);
                        
            %% ripple-locked coactivity for specific ripple condition
            
            % reset the ripple mask
            bMaskRipple     = false(size(trialPhasePerRipple, 1), 1);
            
            % create mask for the ripples
            if strcmp(conditions{iCond}, 'all') % all ripples
                bMaskRipple = true(size(trialPhasePerRipple, 1), 1);
                
            elseif strcmp(conditions{iCond}, 'cueRet') % cue and retrieval
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'});
            elseif strcmp(conditions{iCond}, 'cueRetPrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'cueRetHalf1PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'cueRetHalf2PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & ~bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'cueRetPrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'cueRetHalf1PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'cueRetHalf2PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'C', 'R'}) & ~bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
                
            elseif strcmp(conditions{iCond}, 'ret') % retrieval
                bMaskRipple = contains(trialPhasePerRipple, 'R');
            elseif strcmp(conditions{iCond}, 'retPrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retHalf1PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retHalf2PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & ~bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retPrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retHalf1PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retHalf2PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & ~bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovPrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovHalf1PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovHalf2PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & ~bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovPrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovHalf1PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retMovHalf2PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_Movement') & ~bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovPrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovHalf1PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovHalf2PrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & ~bHalf1PerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovPrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovHalf1PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retNoMovHalf2PrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R_NoMovement') & ~bHalf1PerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retBeLePrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bB4LearnPerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retAfLePrefObjRLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & ~bB4LearnPerRipple & bPrefObjRipples & bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retBeLePrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & bB4LearnPerRipple & bPrefObjRipples & ~bRLwPFRipples;
            elseif strcmp(conditions{iCond}, 'retAfLePrefObjRLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'R') & ~bB4LearnPerRipple & bPrefObjRipples & ~bRLwPFRipples;
                            
            elseif strcmp(conditions{iCond}, 'enc') % (re)encoding
                bMaskRipple = contains(trialPhasePerRipple, 'E');
            elseif strcmp(conditions{iCond}, 'encPrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encHalf1PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encHalf2PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & ~bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encPrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encHalf1PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encHalf2PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & ~bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovPrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovHalf1PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovHalf2PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & ~bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovPrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovHalf1PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encMovHalf2PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_Movement') & ~bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovPrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovHalf1PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovHalf2PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & ~bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovPrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovHalf1PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNoMovHalf2PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E_NoMovement') & ~bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encBeLePrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bB4LearnPerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encAfLePrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & ~bB4LearnPerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encBeLePrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & bB4LearnPerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encAfLePrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, 'E') & ~bB4LearnPerRipple & bPrefObjRipples & ~bOLwPFRipples;
                
        	elseif strcmp(conditions{iCond}, 'encNITI') % (re)encoding and next ITI
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'});
            elseif strcmp(conditions{iCond}, 'encNITIPrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNITIHalf1PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNITIHalf2PrefObjOLwPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & ~bHalf1PerRipple & bPrefObjRipples & bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNITIPrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNITIHalf1PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;
            elseif strcmp(conditions{iCond}, 'encNITIHalf2PrefObjOLoPF')
                bMaskRipple = contains(trialPhasePerRipple, {'E', 'NI'}) & ~bHalf1PerRipple & bPrefObjRipples & ~bOLwPFRipples;

            end
            
            %% 2D coactivity
            
            % skip the coactivity computations, if there are less than two
            % ripples
            if sum(bMaskRipple) < 2
                continue;
            end
            
            % unpack object-cell activity and place-cell activity and
            % reduce to ripples of interest
            objCellActivity     = cell2mat(cfg.allBActiveDuringRipples(bThisObjCell, :)); % number-of-ripples x timepoints
            objCellActivity     = objCellActivity(bMaskRipple, :); % reduce to ripples-of-interest
            placeCellActivity   = cell2mat(cfg.allBActiveDuringRipples(bThisPlaceCell, :));
            placeCellActivity   = placeCellActivity(bMaskRipple, :);
            
            % surrogate object-cell activity based on a circular shift
            randShift               = datasample(1:(size(objCellActivity, 1) - 1), 1);
            objCellActivitySurro    = circshift(objCellActivity, randShift, 1); % shift along the first dimension (i.e., along ripples)
            
            % loop through the different time windows for the object cell
            for iTwoiO = 1:size(cfg.r.param.coac.twoi, 1)
                
                % loop through the different time windows for the place cell
                for iTwoiP = 1:size(cfg.r.param.coac.twoi, 1)
                                        
                    % compute the coactivity
                    if strcmp(cfg.r.param.coac.type, 'z')
                        
                        % z-score
                        zScore2D(iPlaceCell, iTwoiO, iTwoiP)        = LK_CoactivityZScore_20220615(objCellActivity(:, iTwoiO), placeCellActivity(:, iTwoiP));
                        zScore2DSurro(iPlaceCell, iTwoiO, iTwoiP)   = LK_CoactivityZScore_20220615(objCellActivitySurro(:, iTwoiO), placeCellActivity(:, iTwoiP));
                        
                    elseif strcmp(cfg.r.param.coac.type, 'r')
                        
                        % Pearson correlation
                        zScore2D(iPlaceCell, iTwoiO, iTwoiP)        = LK_PearsonCorrelation_20220816(objCellActivity(:, iTwoiO), placeCellActivity(:, iTwoiP));
                        zScore2DSurro(iPlaceCell, iTwoiO, iTwoiP)   = LK_PearsonCorrelation_20220816(objCellActivitySurro(:, iTwoiO), placeCellActivity(:, iTwoiP));
                        
                    end
                    
                    % count the number of ripples and the number of unique
                    % trials
                    numRipples(iPlaceCell, iTwoiO, iTwoiP)      = sum(bMaskRipple);
                    numUniqueTrials(iPlaceCell, iTwoiO, iTwoiP) = numel(unique(trialIdxPerRipple(bMaskRipple)));
                end
            end
        end
        
        %% collect across cells
        
        % collect across time windows
        coacZ2D{iRes}               = zScore2D;
        coacZ2DSurro{iRes}          = zScore2DSurro;
        coacNumRipples{iRes}        = numRipples;
        coacNumUniqueTrials{iRes}   = numUniqueTrials;
    end
    
    %% save
    
    % save important variables
    save(saveFile, 'coacZ2D', 'coacZ2DSurro', 'coacNumRipples', 'coacNumUniqueTrials', ...
        'coacObjCellIdx', 'coacPlaceCellIdx', 'coacObjCellRegion', 'coacPlaceCellRegion', 'coacRippleChanName', 'coacBSameHemiObjRipple', 'coacBSameHemiObjPlace', ...
        '-v7.3');
end

%% stats and figures: coactivity between object cells and place cells during ripples

% close open figures and reset rng for reproducibility
close all;

% different ripple conditions for statistical evaluation
conditions  = {...    
    'cueRetPrefObjRLwPF', 'cueRetPrefObjRLoPF'; ... % cue & retrieval
    'cueRetHalf1PrefObjRLwPF', 'cueRetHalf1PrefObjRLoPF'; ...
    'cueRetHalf2PrefObjRLwPF', 'cueRetHalf2PrefObjRLoPF'; ...
    
    'retPrefObjRLwPF', 'retPrefObjRLoPF'; ... % retrieval
    'retHalf1PrefObjRLwPF', 'retHalf1PrefObjRLoPF'; ...
    'retHalf2PrefObjRLwPF', 'retHalf2PrefObjRLoPF'; ...
    'retMovHalf2PrefObjRLwPF', 'retMovHalf2PrefObjRLoPF'; ... % retrieval movement
    'retNoMovHalf2PrefObjRLwPF', 'retNoMovHalf2PrefObjRLoPF'; ...
    'retBeLePrefObjRLwPF', 'retBeLePrefObjRLoPF'; ... % retrieval, learning
    'retAfLePrefObjRLwPF', 'retAfLePrefObjRLoPF'; ...
    
    'encPrefObjOLwPF', 'encPrefObjOLoPF'; ... % encoding
    'encHalf1PrefObjOLwPF', 'encHalf1PrefObjOLoPF'; ...
    'encHalf2PrefObjOLwPF', 'encHalf2PrefObjOLoPF'; ...
    'encMovHalf2PrefObjOLwPF', 'encMovHalf2PrefObjOLoPF'; ... % encoding movement
    'encNoMovHalf2PrefObjOLwPF', 'encNoMovHalf2PrefObjOLoPF'; ...
    'encBeLePrefObjOLwPF', 'encBeLePrefObjOLoPF'; ... % encoding, learning
    'encAfLePrefObjOLwPF', 'encAfLePrefObjOLoPF'; ...
    
    'encNITIPrefObjOLwPF', 'encNITIPrefObjOLoPF'; ... % encoding and next ITI
    'encNITIHalf1PrefObjOLwPF', 'encNITIHalf1PrefObjOLoPF'; ...
    'encNITIHalf2PrefObjOLwPF', 'encNITIHalf2PrefObjOLoPF'; ...
    };

% text file with results
txtFileName = strcat(cfg.r.paths.save, 'statsCoac_', cfg.r.param.coac.type, '_objCellPlaceCell_20220511');
if exist(txtFileName, 'file')
    txtFile = fopen(txtFileName, 'a+');
else
    txtFile = fopen(txtFileName, 'wt');
end

% statistics and plots for coactivity per condition
for iCond = 1:size(conditions, 1)
    
    % report condition
    fprintf('\n\n========================================================== CONDITION #%d: "%s" (vs. "%s").\n', iCond, conditions{iCond, 1}, conditions{iCond, 2});
    
    % reset rng for reproducibility
    rng(numel(conditions{iCond, 1}));
    
    % load coactivity results from this condition
    resCoac = load(strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', conditions{iCond, 1}), 'coacZ2D', 'coacZ2DSurro', ...
        'coacNumRipples', 'coacNumUniqueTrials', 'coacObjCellIdx');
    
    % load coacitivity results from the corresponding condition in which
    % the location of the preferred object is outside the place field
    outCoac = load(strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', conditions{iCond, 2}), 'coacZ2D');
    
    % unfold 2D coactivity z-scores from this condition
    unfCoacZ2D              = cell2mat(resCoac.coacZ2D); % cell*cell*rippleChannel x twoiO x twoiP
    unfCoacZ2DSurro         = cell2mat(resCoac.coacZ2DSurro);
    unfCoacNumRipples       = cell2mat(resCoac.coacNumRipples); % how many ripples contribute to the coactivity scores
    unfCoacObjCellIdx       = cell2mat(resCoac.coacObjCellIdx); % indices of the object cells contributing to the coactivity scores
    
    % unfold 2D coactivity z-scores from the contrast condition
    unfOutCoacZ2D           = cell2mat(outCoac.coacZ2D);
    
    % sanity check
    if size(unfCoacZ2D, 1) ~= size(unfCoacObjCellIdx, 1)
        error('Problem with the data dimensions in "unfCoacZ2D".');
    end
    
    % permute the dimensions so that samples are the 3rd dimension
    unfCoacZ2D              = permute(unfCoacZ2D, [2, 3, 1]); % twoiO x twoiP x cell*cell*rippleChannel
    unfCoacZ2DSurro         = permute(unfCoacZ2DSurro, [2, 3, 1]);
    unfCoacNumRipples       = permute(unfCoacNumRipples, [2, 3, 1]);
    unfOutCoacZ2D           = permute(unfOutCoacZ2D, [2, 3, 1]); % contrast condition
    
    %% target and baseline coactivity
    
    % coactivity during the target period
    logIdxTarget    	= cfg.r.param.coac.time >= min(cfg.r.param.coac.target) & cfg.r.param.coac.time <= max(cfg.r.param.coac.target);
    coacTarget          = unfCoacZ2D(logIdxTarget, logIdxTarget, :);
    coacTargetSurro     = unfCoacZ2DSurro(logIdxTarget, logIdxTarget, :);
    
    % coactivity during the pre-baseline period
    logIdxBasePre     	= cfg.r.param.coac.time >= min(cfg.r.param.coac.basePre) & cfg.r.param.coac.time <= max(cfg.r.param.coac.basePre);
    coacBasePre         = unfCoacZ2D(logIdxBasePre, logIdxBasePre, :);
    
    % coactivity during the post-baseline period
    logIdxBasePost      = cfg.r.param.coac.time >= min(cfg.r.param.coac.basePost) & cfg.r.param.coac.time <= max(cfg.r.param.coac.basePost);
    coacBasePost     	= unfCoacZ2D(logIdxBasePost, logIdxBasePost, :);
    
    % mean coactivity during the pre- and post-baseline period
    coacBase        	= mean(cat(4, coacBasePre, coacBasePost), 4, 'omitnan'); % coactivity during the pre- and post-period (averaged)
    
    % contrast coactivity during the target period
    outCoacTarget       = unfOutCoacZ2D(logIdxTarget, logIdxTarget, :);
    
    %% statistics: target coactivity > 0
    
    % test the coactivity z-values from the target period against zero
    dat                 = [];
    dat.mat             = coacTarget;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8; % type of connectivity (4 or 8)
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTestZero        = LK_3DPermutationTest_OneSampleAgainstZero_20220302(dat);
    
    % save stats
    fprintf(txtFile, '%s,targetVS0,stat=%.6f,rank=%.6f,p=%.6f\n', conditions{iCond, 1}, max(permTestZero.clusStat), max(permTestZero.clusStatRank), min(permTestZero.clusStatP));
    save(strcat(cfg.r.paths.save, 'permTestZero_', cfg.r.param.coac.type, '_', conditions{iCond, 1}, '.mat'), 'permTestZero');
    
    %% statistics: target coactivity > surrogate target coactivity
    
    % sanity check
    if ~all(isnan(coacTarget(:)) == isnan(coacTargetSurro(:)))
        error('Different distribution of NaNs in "coacTarget" and "coacTargetSurro".');
    end
    
    % test the coactivity z-values from the target period minus the
    % surrogate coactivity z-values from the target period against zero
    dat                 = [];
    dat.mat             = coacTarget - coacTargetSurro;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8;
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTestSurro       = LK_3DPermutationTest_OneSampleAgainstZero_20220302(dat);
    
    % save stats
    fprintf(txtFile, '%s,targetVSsurro,stat=%.6f,rank=%.6f,p=%.6f\n', conditions{iCond, 1}, max(permTestSurro.clusStat), max(permTestSurro.clusStatRank), min(permTestSurro.clusStatP));
    save(strcat(cfg.r.paths.save, 'permTestSurro_', cfg.r.param.coac.type, '_', conditions{iCond, 1}, '.mat'), 'permTestSurro');
    
    %% statistics: target coactivity > baseline coactivity
    
    % test the target coactivity z-values against the baseline coactivity
    % z-values
    dat                 = [];
    dat.mat1            = coacTarget;
    dat.mat2            = coacBase;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8;
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTestBase        = LK_3DPermutationTest_TwoConnectedSamples_20220901(dat);
    
    % save stats
    fprintf(txtFile, '%s,targetVSbaseline,stat=%.6f,rank=%.6f,p=%.6f\n', conditions{iCond, 1}, max(permTestBase.clusStat), max(permTestBase.clusStatRank), min(permTestBase.clusStatP));
    save(strcat(cfg.r.paths.save, 'permTestBase_', cfg.r.param.coac.type, '_', conditions{iCond, 1}, '.mat'), 'permTestBase');
    
    %% statistics: within place field vs. outside place field
    
    % test the coactivities-of-interest against the contrast-coactivities
    dat                 = [];
    dat.mat1            = coacTarget;
    dat.mat2            = outCoacTarget;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8;
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTestOut         = LK_3DPermutationTest_TwoSamples_20220302(dat);
    
    % save stats
    fprintf(txtFile, '%s,targetVSoutTarget,stat=%.6f,rank=%.6f,p=%.6f\n', conditions{iCond, 1}, max(permTestOut.clusStat), max(permTestOut.clusStatRank), min(permTestOut.clusStatP));
    save(strcat(cfg.r.paths.save, 'permTestOut_', cfg.r.param.coac.type, '_', conditions{iCond, 1}, '.mat'), 'permTestOut');
    
    %% figure: 2D coactivity (evaluted using different tests)
    
    % mean across cell-cell-rippleChannels and other variables
    t       = cfg.r.param.coac.time(logIdxTarget); % time
    groups  = {'target', [0, 0, 0]; 'targetVSSurro', [0, 0.5, 0]; 'targetVSBase', [1, 0, 0]; 'targetVSoutTarget', [0, 0, 1]};
    m       = {mean(coacTarget, 3, 'omitnan'), mean(coacTarget - coacTargetSurro, 3, 'omitnan'), ...
        mean(coacTarget, 3, 'omitnan') - mean(coacBase, 3, 'omitnan'), mean(coacTarget, 3, 'omitnan') - mean(outCoacTarget, 3, 'omitnan')};
    sig     = {permTestZero.logIdxSigClus, permTestSurro.logIdxSigClus, permTestBase.logIdxSigClus, permTestOut.logIdxSigClus};
    pValue  = {min(permTestZero.clusStatP), min(permTestSurro.clusStatP), min(permTestBase.clusStatP), min(permTestOut.clusStatP)};
    
    % loop through the different data groups
    for iGroup = 1:size(groups, 1)
        
        % create figure
        dt              = [];
        dt.t            = t;
        dt.m            = m{iGroup};
        dt.sig          = sig{iGroup};
        dt.pValue       = pValue{iGroup};
        dt.color        = groups{iGroup, 2};
        dt.condition    = conditions{iCond, 1};
        dt.coac         = cfg.r.param.coac;
        f = LK_PlotCoactivity_20231015(dt);
        
        % save figure
        LK_print(f, strcat(cfg.r.paths.save, 'coacZ2D_', cfg.r.param.coac.type, '_', groups{iGroup, 1}, '_objCellPlaceCell_', conditions{iCond, 1}, '_20220616'), '-dpng', '-r300');
    
        % save figure data
        if strcmp(conditions{iCond, 1}, 'retPrefObjRLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7c_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7c_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7c_right'), 'dt');
            end
        elseif strcmp(conditions{iCond, 1}, 'retHalf1PrefObjRLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7d_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7d_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7d_right'), 'dt');
            end
        elseif strcmp(conditions{iCond, 1}, 'retHalf2PrefObjRLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7e_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7e_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7e_right'), 'dt');
            end
        elseif strcmp(conditions{iCond, 1}, 'encPrefObjOLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7f_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7f_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7f_right'), 'dt');
            end
        elseif strcmp(conditions{iCond, 1}, 'encHalf1PrefObjOLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7g_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7g_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7g_right'), 'dt');
            end
        elseif strcmp(conditions{iCond, 1}, 'encHalf2PrefObjOLwPF')
            if strcmp(groups{iGroup, 1}, 'target')
                save(strcat(cfg.r.paths.save, 'Fig_7h_left'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSBase')
                save(strcat(cfg.r.paths.save, 'Fig_7h_middle'), 'dt');
            elseif strcmp(groups{iGroup, 1}, 'targetVSoutTarget')
                save(strcat(cfg.r.paths.save, 'Fig_7h_right'), 'dt');
            end
        end
    end
    
    %% figure: number of cell-cell-ripple channel-combinations and number of ripples

    % loop through different data groups
    groups  = {'numCombinations'; 'numRipples'};
    for iGroup = 1:numel(groups)
        
        % select data depending on data group
        if strcmp(groups{iGroup}, 'numCombinations')
            data                    = sum(~isnan(coacTarget), 3, 'omitnan'); % number of cell-cell-ripple channel-combinations
        elseif strcmp(groups{iGroup}, 'numRipples')
            data                    = unfCoacNumRipples; % twoiO x twoiP x cell*cell*rippleChannel combinations
            data(isnan(unfCoacZ2D)) = nan;
            data                    = sum(data(logIdxTarget, logIdxTarget, :), 3, 'omitnan'); % number of ripples
        end
        
        % create figure
        dt      = [];
        dt.t    = t;
        dt.data = data;
        f = LK_PlotCoactivityAdditionalInformation_20231015(dt);
        
        % save the figure
        LK_print(f, strcat(cfg.r.paths.save, groups{iGroup}, '_', cfg.r.param.coac.type, '_target_objCellPlaceCell_', conditions{iCond}, '_20220511'), '-dpng', '-r300');
    end
    
    % close open figures
    close all;
end

% close text file
fclose('all');

%% comparison between retrieval and re-encoding with normalization

% comparisons to make
comparisons = {...
    'retHalf2PrefObjRLwPF', 'encHalf2PrefObjOLwPF'; ...
    };

% target time period
logIdxTarget    = cfg.r.param.coac.time >= min(cfg.r.param.coac.target) & cfg.r.param.coac.time <= max(cfg.r.param.coac.target);
t               = cfg.r.param.coac.time(logIdxTarget);

% loop through comparisons
for iComp = 1:size(comparisons, 1)
    
    % report comparison
    fprintf('\n\n========================================================== COMPARISON #%d: "%s" vs "%s".\n', iComp, comparisons{iComp, 1}, comparisons{iComp, 2});
    
    % reset rng for reproducibility
    rng(numel(comparisons{iComp, 1}));
    
    % load coactivity results
    resCoac1    = load(strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', comparisons{iComp, 1}), 'coacZ2D');
    resCoac2    = load(strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', comparisons{iComp, 2}), 'coacZ2D');
    
    % unfold 2D coactivity
    unfCoac1  	= cell2mat(resCoac1.coacZ2D); % cell*cell*rippleChannel x twoiO x twoiP
    unfCoac2  	= cell2mat(resCoac2.coacZ2D);
        
    % permute the dimensions so that samples are the 3rd dimension
    unfCoac1   	= permute(unfCoac1, [2, 3, 1]); % twoiO x twoiP x cell*cell*rippleChannel
    unfCoac2    = permute(unfCoac2, [2, 3, 1]);
    
    % restrict to target time window
    unfCoac1    = unfCoac1(logIdxTarget, logIdxTarget, :);
    unfCoac2    = unfCoac2(logIdxTarget, logIdxTarget, :);
    
    %% normalize the coactivity maps
    % the coactivity map of each cell pair then ranges between 0 and 1
        
    % normalize the first coactivity map
    unfCoacNorm1      	= unfCoac1 - repmat(min(min(unfCoac1, [], 1), [], 2), size(unfCoac1, 1), size(unfCoac1, 2), 1); % subtract minimum
    unfCoacNorm1     	= unfCoacNorm1 ./ repmat(max(max(unfCoacNorm1, [], 1), [], 2), size(unfCoacNorm1, 1), size(unfCoacNorm1, 2), 1);
    
    % normalize the second coactivity map
    unfCoacNorm2      	= unfCoac2 - repmat(min(min(unfCoac2, [], 1), [], 2), size(unfCoac2, 1), size(unfCoac2, 2), 1);
    unfCoacNorm2      	= unfCoacNorm2 ./ repmat(max(max(unfCoacNorm2, [], 1), [], 2), size(unfCoacNorm2, 1), size(unfCoacNorm2, 2), 1);
    
    % cluster-based permutation test for condition 1 > 2
    dat                 = [];
    dat.mat1            = unfCoacNorm1;
    dat.mat2            = unfCoacNorm2;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8;
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTest1vs2        = LK_3DPermutationTest_TwoConnectedSamples_20220901(dat);
    
    % cluster-based permutation test for condition 2 > 1
    dat                 = [];
    dat.mat1            = unfCoacNorm2;
    dat.mat2            = unfCoacNorm1;
    dat.firstalpha      = 0.05;
    dat.secondalpha     = 0.05;
    dat.direction       = 'pos';
    dat.conn            = 8;
    dat.numSurrogates   = cfg.r.param.coac.numSurrogates;
    permTest2vs1        = LK_3DPermutationTest_TwoConnectedSamples_20220901(dat);
    
    %% figure: 2D coactivity for the comparisons
    
    % mean across cell-cell-rippleChannels
    groups  = {'1vs2', [0, 0, 0]; '2vs1', [0, 0, 0]}; % contrast, color
    m       = {mean(unfCoacNorm1, 3, 'omitnan') - mean(unfCoacNorm2, 3, 'omitnan'); ...
        mean(unfCoacNorm2, 3, 'omitnan') - mean(unfCoacNorm1, 3, 'omitnan')}; % mean coactivity maps
    sig     = {permTest1vs2.logIdxSigClus; permTest2vs1.logIdxSigClus}; % significant time areas
    pValue  = {min(permTest1vs2.clusStatP); min(permTest2vs1.clusStatP)}; % p-values from cluster-based permutation tests
    
    % loop through the different data groups
    for iGroup = 1:size(groups, 1)
        
        %% create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 7.4, 6]);
        
        % plot data for this comparison
        axes('units', 'centimeters', 'position', [2, 1.8, 4, 4]);
        hold on;
        imagesc(t, t, m{iGroup}); % mean
        contour(t, t, sig{iGroup}, 1, '-', 'Color', [1, 1, 1], 'LineWidth', 1); % significance
        colormap parula;
        xline(0, '--', 'Color', [0, 0, 0]);
        yline(0, '--', 'Color', [0, 0, 0]);
        set(gca, 'xtick', [-0.2, 0, 0.2], 'xticklabel', {'-0.2', '0', '0.2'}, 'ytick', [-0.2, 0, 0.2], 'yticklabel', {'-0.2', '0', '0.2'}, 'box', 'on');
        axis equal;
        % indicate p-value
        if pValue{iGroup} < 0.001
            t1 = text(0.03, 0.925, '\itP\rm<0.001', 'units', 'normalized', 'Color', [1, 1, 1]);
        else
            t1 = text(0.03, 0.925, ['\itP\rm=', num2str(pValue{iGroup}, '%.3f')], 'units', 'normalized', 'Color', [1, 1, 1]);
        end
        % indicate type of group contrast
        tmpAx = get(gca);
        plot([min(tmpAx.XLim), max(tmpAx.XLim), max(tmpAx.XLim), min(tmpAx.XLim), min(tmpAx.XLim)], ...
            [min(tmpAx.YLim), min(tmpAx.YLim), max(tmpAx.YLim), max(tmpAx.YLim), min(tmpAx.YLim)], '-', 'linewidth', 1.5, 'Color', groups{iGroup, 2});
        % colorbar limits
        myCLim = round(max(abs(m{iGroup}(:))), 1) * [-1, 1];
        caxis(myCLim);
        % axis labeling
        xl = xlabel({'Time during ripple (s)', 'Place cell'});
        yl = ylabel({'Object cell', 'Time during ripple (s)'});
        cb = colorbar('eastoutside', 'units', 'centimeters', 'position', [6.15, 3.05, 0.35, 1.5]);
        cb.Ticks = myCLim;
        cb.TickLabels = {num2str(min(myCLim), '%.1f'), num2str(max(myCLim), '%.1f')};
        cb.LineWidth = 1;
        title(cb, [cfg.r.param.coac.type, '_{rel}']);
        set([gca, xl, yl, t1], 'fontunits', 'centimeters', 'fontsize', 0.4);
        drawnow;
        
        % save the figure
        LK_print(f, strcat(cfg.r.paths.save, 'coacZ2D_', cfg.r.param.coac.type, '_', groups{iGroup, 1}, '_objCellPlaceCell_', comparisons{iComp, 1}, '_', comparisons{iComp, 2}, '_20220616'), '-dpng', '-r300');
    end
end

%% detailed information about coactivity maps for select conditions
% - cell pair-specific z-scores underlying the peaks in the coactivity maps
% - brain regions of the cell pairs
% - cell pair-specific coactivity maps

% select conditions for which to produce combination-specific figures
selConditions  	= {'retPrefObjRLwPF'; 'retHalf1PrefObjRLwPF'; 'retHalf2PrefObjRLwPF'; 'encPrefObjOLwPF'; 'encHalf1PrefObjOLwPF'; 'encHalf2PrefObjOLwPF'};

% target time period
logIdxTarget  	= cfg.r.param.coac.time >= min(cfg.r.param.coac.target) & cfg.r.param.coac.time <= max(cfg.r.param.coac.target);
t               = cfg.r.param.coac.time(logIdxTarget);

% brain regions for the figure
regions4Fig     = {'AMY'; 'EC'; 'FG'; 'HC'; 'Insula'; 'PHC'; 'TP'; 'VC'};

% loop through select conditions
for iSC = 1:numel(selConditions)
    
    % load coactivity results from this condition
    resCoac                 = load(strcat(cfg.r.paths.save, 'resCoac_', cfg.r.param.coac.type, '_', selConditions{iSC}), ...
        'coacZ2D', 'coacObjCellIdx', 'coacPlaceCellIdx', 'coacRippleChanName', 'coacObjCellRegion', 'coacPlaceCellRegion');
    
    % unfold the coactivity results
    unfCoacZ2D              = cell2mat(resCoac.coacZ2D);
    unfCoacObjCellIdx       = cell2mat(resCoac.coacObjCellIdx);
    unfCoacPlaceCellIdx     = cell2mat(resCoac.coacPlaceCellIdx);
    unfCoacRippleChanName   = LK_unwrapNestedCell(resCoac.coacRippleChanName);
    unfCoacObjCellRegion    = LK_unwrapNestedCell(resCoac.coacObjCellRegion);
    unfCoacPlaceCellRegion  = LK_unwrapNestedCell(resCoac.coacPlaceCellRegion);
    
    % permute the coactivity results and cut to target time window
    unfCoacZ2D   = permute(unfCoacZ2D, [2, 3, 1]); % twoiO x twoiP x cell*cell*rippleChannel
    unfCoacZ2D   = unfCoacZ2D(logIdxTarget, logIdxTarget, :);
    
    %% individual z-scores underlying the peaks
    
    % mean of the coactivity maps and its peak
    m                           = mean(unfCoacZ2D, 3, 'omitnan');
    [peakVal, peakIdx]          = max(m(:));
    [peakIdxRow, peakIdxColumn] = ind2sub(size(m), peakIdx);
    
    % individual z-values of the peak
    indZ    = squeeze(unfCoacZ2D(peakIdxRow, peakIdxColumn, :));
    
    % sanity check
    if peakVal ~= m(peakIdxRow, peakIdxColumn)
        error('Incorrect extraction of peak index.');
    end
    
    % histogram of individual z-values
    dat             = [];
    dat.figPosition = [2, 2, 6, 6];
    dat.axPosition  = [1.5, 1.5, 4, 4];
    dat.data        = indZ;
    dat.xline       = mean(indZ, 'omitnan');
    dat.binEdges    = -4:0.5:4;
    dat.xlabel      = 'Coactivity z';
    dat.ylabel      = 'Count';
    dat.title       = ['(p=', num2str(t(peakIdxColumn) * 1000, '%.0f'), 'ms/o=', num2str(t(peakIdxRow) * 1000, '%.0f'), 'ms)'];
    dat.fontSize    = 0.4;
    f               = LK_HistogramFigure_20220321(dat);
    
    % save figure
    LK_print(f, strcat(cfg.r.paths.save, 'coacZ2D_', cfg.r.param.coac.type, '_indZAtPeak_objCellPlaceCell_', selConditions{iSC, 1}, '_20220616'), '-dpng', '-r300');
    
    %% brain regions of the cell pairs contributing to the peaks
    
    % regions of object and place cells contributing to the main result
    numCellsPerRegion   = nan(numel(regions4Fig), numel(regions4Fig));
    for iO = 1:numel(regions4Fig)
        for iP = 1:numel(regions4Fig)
            numCellsPerRegion(iO, iP)   = sum(strcmp(unfCoacObjCellRegion(~isnan(indZ)), regions4Fig{iO}) & strcmp(unfCoacPlaceCellRegion(~isnan(indZ)), regions4Fig{iP}));
        end
    end
    
    % create figure
    f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
    axes('units', 'centimeters', 'position', [1.1, 1.6, 3.75, 3.75]);
    imagesc(numCellsPerRegion);
    colormap([[1, 1, 1]; cool]);
    cb = colorbar('units', 'centimeters', 'position', [5, 2.5, 0.2, 2]);
    cb.FontSize = 12;
    cb.Ticks = [0, max(numCellsPerRegion(:))];
    xl = xlabel('Place-cell region');
    yl = ylabel('Object-cell region');
    set(gca, 'xtick', 1:numel(regions4Fig), 'xticklabel', cellfun(@(x) x(1), cellstr(regions4Fig)), 'xticklabelrotation', 0, ...
        'ytick', 1:numel(regions4Fig), 'yticklabel', cellfun(@(x) x(1), cellstr(regions4Fig)));
    set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
    % save figure
    LK_print(f, strcat(cfg.r.paths.save, 'numCellsPerRegionAtPeak_', selConditions{iSC}, '_20220619'), '-dpng', '-r300');
    
    %% individual coactivity maps
    
    % loop through combinations and plot the coactivity maps
    for iComb = 1:size(unfCoacZ2D, 3)
        
        % skip if you do not want to plot individual coactivity maps
        if cfg.r.param.coac.bPlotIndividuals == false || all(all(isnan(unfCoacZ2D(:, :, iComb))))
            continue;
        end
        
        % create figure
        f = figure('units', 'centimeters', 'position', [2, 2, 7.4, 6.2], 'visible', 'off');
        
        % plot data for this combination
        axes('units', 'centimeters', 'position', [2, 1.8, 4, 4]);
        hold on;
        imagesc(t, t, unfCoacZ2D(:, :, iComb));
        xline(0, '--', 'Color', [0, 0, 0]);
        yline(0, '--', 'Color', [0, 0, 0]);
        set(gca, 'xtick', [-0.2, 0, 0.2], 'xticklabel', {'-0.2', '0', '0.2'}, 'ytick', [-0.2, 0, 0.2], 'yticklabel', {'-0.2', '0', '0.2'}, 'box', 'on', 'Layer', 'Top');
        axis equal;
        % axis labeling
        xl = xlabel({'Time during ripple (s)', 'Place cell'});
        yl = ylabel({'Object cell', 'Time during ripple (s)'});
        % colorbar
        colormap parula;
        cLim = [min(min(unfCoacZ2D(:, :, iComb))), max(max(unfCoacZ2D(:, :, iComb)))];
        if all(isnan(cLim))
            cLim    = [-1, 1];
        elseif min(cLim) == max(cLim)
            cLim    = [min(cLim) - 1, max(cLim) + 1];
        end
        caxis(cLim);
        cb = colorbar('eastoutside', 'units', 'centimeters', 'position', [6.15, 3.05, 0.35, 1.5]);
        cb.Ticks = cLim;
        cb.TickLabels = {num2str(min(cLim), '%.1f'), num2str(max(cLim), '%.1f')};
        cb.LineWidth = 1;
        title(cb, cfg.r.param.coac.type);
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
        title(sprintf('%d / %d / %d - %d / %d / %d - %s', unfCoacObjCellIdx(iComb, :), unfCoacPlaceCellIdx(iComb, :), unfCoacRippleChanName{iComb}), ...
            'fontunits', 'centimeters', 'fontsize', 0.2, 'fontweight', 'normal', 'units', 'normalized', 'position', [0.5, 1.025]);
        drawnow;
        
        % save the figure
        LK_print(f, strcat(cfg.r.paths.pics, 'coacZ2D_', cfg.r.param.coac.type, ...
            '_objCell_', num2str(unfCoacObjCellIdx(iComb, 1)), '_', num2str(unfCoacObjCellIdx(iComb, 2)), '_', num2str(unfCoacObjCellIdx(iComb, 3)), ...
            '_placeCell_', num2str(unfCoacPlaceCellIdx(iComb, 1)), '_', num2str(unfCoacPlaceCellIdx(iComb, 2)), '_', num2str(unfCoacPlaceCellIdx(iComb, 3)), ...
            '_', unfCoacRippleChanName{iComb}, '_', selConditions{iSC}), ...
            '-dpng', '-r300');
        close all;
    end
end

%% output

% create output
out     = [];
