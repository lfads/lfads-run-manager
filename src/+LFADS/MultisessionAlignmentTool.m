classdef MultisessionAlignmentTool < handle
    
    properties(SetAccess=protected)
        nFactors
        spikeBinMs % spike bins taken from run.params.spikeBinMs
        countsByDataset % nDatasets x 1 cell of nChannels_ds x nTimepoints x nTrials
        nDatasets
        conditions % full list of unique conditions across datasets 
        nChannelsByDataset % nDatasets x 1 vector of channel counts
        nTrialsByDataset % nDatasets x 1 vector of trial counts
        conditionIdxByDataset % nDatasets x vector of condition idenity
        nConditions
        nTime
        conditionAvgsByDataset % nDatasets x 1 cell of nChannels_ds x nTimepoints x nConditions tensors of condition averages
        
        datasetNames % nDatasets x 1 cell of dataset names
        
        alignmentMatrices % nDatasets x 1 cell of nChannelsByDataset(iDS) x nFactors      
        alignmentBiases % nDatasets x 1 cell of nChannelsByDataset(iDS) x 1
        
        % results of prepareAlignmentMatricesUsingTrialAveragedPCR
        pcAvg_allDatasets % nFactors x nTimepoints x nConditions tensor of PCs using all datasets
        pcAvg_reconstructionByDataset % nFactors x nTimepoints x nConditions x nDatasets tensor of best reconstruction of pcAvg_allDatasets through the alignment matrices
    end

    methods
        function tool = MultisessionAlignmentTool(run, seqData, datasetNames)
            %
            % Parameters
            % ------------
            % seqData : `nDatasets` cell of struct arrays of sequence data
            %   Sequence data for each dataset as returned by `convertDatasetToSequenceStruct`
            %
            
            tool.spikeBinMs = run.params.spikeBinMs;
            tool.nDatasets = numel(seqData);   
            tool.nFactors = run.params.c_factors_dim;
            tool.datasetNames = LFADS.Utils.makecol(datasetNames);
            
            % compute all unique conditions across datasets
            condField = 'conditionId';
            c = cell(tool.nDatasets, 1);
            conditionsEachTrial = cell(tool.nDatasets, 1);
            for iDS = 1:tool.nDatasets
               conditionsEachTrial{iDS} = {seqData{iDS}.(condField)}';
               if isscalar(conditionsEachTrial{iDS}{1})
                   conditionsEachTrial{iDS} = cell2mat(conditionsEachTrial{iDS});
               end
               
               if iscell(conditionsEachTrial{iDS})
                   c{iDS} = removeempty(unique(conditionsEachTrial{iDS}));
               else
                   c{iDS} = unique(removenan(conditionsEachTrial{iDS}));
               end
            end
            conditions = unique(cat(1, c{:}));
            if ~isnumeric(conditions) && ~iscellstr(conditions)
                error('conditionId field must contain numeric or string values');
            end
            tool.conditions = conditions;
            tool.nConditions = numel(conditions);
            
            % lookup table of each trial to the unique-ified condition id list
            tool.conditionIdxByDataset = cell(tool.nDatasets, 1);
            for iDS = 1:tool.nDatasets
                [~, tool.conditionIdxByDataset{iDS}] = ismember(conditionsEachTrial{iDS}, conditions);
            end
            
            % seqData is a struct array over trials, where y contains
            % nChannels x nTimepoints
            % convert each to a nTrials x nChannels x nTimepoints tensor,
            % rebinned at run.params.spikeBinMs
            tool.countsByDataset = cell(tool.nDatasets, 1);
            [tool.nChannelsByDataset, tool.nChannelsByDataset] = deal(nan(tool.nDatasets, 1));
            
            % populate single-trial data in countsByDataset
            for iDS = 1:tool.nDatasets
                origBinMs = seqData{iDS}.binWidthMs;
                rebinBy = tool.spikeBinMs / origBinMs;
                assert(abs(rebinBy - round(rebinBy)) < 1e-6, 'Ratio of new spike bin ms to original spike bin ms must be integer');
                rebinBy = round(rebinBy);
                countsOriginalBinning = cat(3, seqData{iDS}.y); % C x T x tRials
                
                tool.nChannelsByDataset(iDS) = size(countsOriginalBinning, 1);
                tool.nTrialsByDataset(iDS) = size(countsOriginalBinning, 3);
                
                % rebin to tool.spikeBinMs
                C = tool.nChannelsByDataset(iDS);
                R = tool.nTrialsByDataset(iDS);
                Ttruncated = floor(size(countsOriginalBinning, 2) / rebinBy) * rebinBy;
                
                countsRebinned = reshape(sum(reshape(...
                    countsOriginalBinning(:, 1:Ttruncated, :), ...
                    [C rebinBy Ttruncated/rebinBy R]), 2), [C Ttruncated/rebinBy R]);
                    
                if iDS == 1
                    tool.nTime = size(countsRebinned, 2);
                else
                    assert(size(countsRebinned, 2) == tool.nTime, 'Number of timepoints (after rebinning) of all datasets must match');
                end
                
                tool.countsByDataset{iDS} = countsRebinned;
            end  
            
            % compute trial-averages for each dataset into conditionAvgsByDataset
            tool.conditionAvgsByDataset = cell(tool.nDatasets, 1);
            for iDS = 1:tool.nDatasets
                condIdx = tool.conditionIdxByDataset{iDS};
                tool.conditionAvgsByDataset{iDS} = nan(tool.nChannelsByDataset(iDS), tool.nTime, tool.nConditions);
                
                for iC = 1:tool.nConditions
                    tool.conditionAvgsByDataset{iDS}(:, :, iC) = nanmean(tool.countsByDataset{iDS}(:, :, condIdx == iC), 3);
                end
            end
            
            function a = removenan(a)
                a = a(~isnan(a));
            end
            
            function a = removeempty(a)
                a = a(~cellfun(@isempty, a));
            end
        end
        
        function [alignmentMatrices, alignmentBiases] = computeAlignmentMatricesUsingTrialAveragedPCR(tool)
            % Prepares alignment matrices to seed the stitching process when 
            % using multiple days of sequence data for LFADS input file generation. 
            % Generate alignment matrices which specify the initial guess at the 
            % encoder matrix that converts neural activity from each dataset 
            % to a common set of factors (for stitching). Each alignment matrix
            % should be nNeurons (that session) x nFactors.
            % 
            % Theis implementation computes trial-averages (averaging
            % all trials with the same conditionId label) for each neuron
            % in each session. The trial-averages are then assembled into a
            % large nNeuronsTotal x (nConditions x time) matrix. The top
            % nFactors PCs of this matrix are computed (as linear
            % combinations of neurons). For each session, we then regress
            % the nNeuronsThisSession neurons against the top nFactors PCs. 
            % The alignment matrix is the matrix of regression
            % coefficients.
            %
            % If you wish to exclude a trial from the alignment matrix
            % calculations, set conditionId to NaN or ''
            %
            % Parameters
            % ------------
            % seqData : `nDatasets` cell of struct arrays of sequence data
            %   Sequence data for each dataset as returned by `convertDatasetToSequenceStruct`
            %
            % Returns
            % ----------
            % alignmentMatrices : `nDatasets` cell of `nNeuronsThisSession` x `nFactors` matrices
            %   For each dataset, an initial guess at the encoder matrices which maps `nNeuronsThisSession` (for that dataset) to a
            %   common set of `nFactors` (up to you to pick this). Seeding this well helps the stitching process. Typically,
            %   PC regression can provide a reasonable set of guesses.
            
            all_data_tensor = cat(1, tool.conditionAvgsByDataset{:}); % nChannelsTotal x nTime x nConditions
            all_data = all_data_tensor(:, :);
            which_dataset = cell2mat(arrayfun(@(n, idx) idx*ones(n, 1), tool.nChannelsByDataset, (1:tool.nDatasets)', 'UniformOutput', false));
            
            all_data_means = nanmean(all_data, 2);
            all_data_centered = bsxfun(@minus, all_data, all_data_means);
            
            % apply PCA
            nFactorsRequested = min(tool.nFactors, size(all_data_centered, 1)); % in case nFactors > dataset num channels
            try
                [pca_proj_mat, pc_data] = pca(all_data_centered', 'Rows', 'pairwise', 'NumComponents', nFactorsRequested);
            catch
                [pca_proj_mat, pc_data] = pca(all_data_centered', 'Rows', 'complete', 'NumComponents', nFactorsRequested);
            end
            % pca_proj_mat is nNeuronsTotal x nFactorsRequested
            % pc_data is nTime*nConditions x nFactorsRequested
            
            [alignmentMatrices, alignmentBiases] = deal(cell(tool.nDatasets, 1));
            tool.pcAvg_reconstructionByDataset = zeros(tool.nFactors, tool.nTime, tool.nConditions, tool.nDatasets);
            if tool.nDatasets == 1
                % single dataset alignment matrix just for dimensionality
                % reduction
                % just use the pca projection directly
                alignmentMatrices{1} = pca_proj_mat;
                alignmentBiases{1} = all_data_means;
                
                % pca_proj_mat is (nTime*nConditions) x nFactorsRequested
                % pcAvg_reconstructionByDataset is nFactors x nTime x nConditions x 1 
                tool.pcAvg_reconstructionByDataset(1:nFactorsRequested, :, :, 1) = reshape(pc_data', [nFactorsRequested, tool.nTime, tool.nConditions]);
                
            else
                % multi-dataset PC regression to reconstruct the global PCs
                % from each dataset's neurons alone
                
                % project all data into pca space
                % dim_reduced_data will be nFactors x nTime*nConditions
                dim_reduced_data = pca_proj_mat' * all_data_centered;
                tool.pcAvg_allDatasets = reshape(dim_reduced_data, [tool.nFactors, tool.nTime, tool.nConditions]);

                % get a mapping from each day to the lowD space
                for iDS = 1:tool.nDatasets
                    % nChannelsByDataset(iDS) x (nTime*nConditions)
                    this_dataset_data = all_data(which_dataset==iDS, :);
                    % nChannelsByDataset(iDS) x 1
                    this_dataset_means = all_data_means(which_dataset==iDS);

                    % figure out which timepoints are valid
                    tMask = ~any(isnan(this_dataset_data), 1) & ~any(isnan(dim_reduced_data), 1);
                    % nFactors x (nTime*nConditions)
                    dim_reduced_data_this = bsxfun(@minus, dim_reduced_data(:, tMask), ...
                        nanmean(dim_reduced_data(:, tMask), 2));
                    % nChannelsByDataset(iDS) x (nTime*nConditions)
                    this_dataset_centered = bsxfun(@minus, this_dataset_data(:, tMask), ...
                        nanmean(this_dataset_data(:, tMask), 2));

                    % regress this day's data against the global PCs -
                    % nChannelsByDataset(iDS) x nFactors
                   alignmentMatrices{iDS} = (this_dataset_centered' \ dim_reduced_data_this');
                    % and set mean as it will be subtracted from the data
                    % before projecting by the alignment matrix
                    % nChannelsByDataset(iDS) x 1
                   alignmentBiases{iDS} = squeeze(this_dataset_means);

                    if any(isnan(alignmentMatrices{iDS}(:)))
                        error('NaNs in the the alignment matrix');
                    end
                    if any(isnan(alignmentBiases{iDS}(:)))
                        error('NaNs in the the alignment biases');
                    end
                
                    % expand alignment matrices to include additional factors if
                    % too many requested
                    nCol = size(alignmentMatrices{iDS}, 2);
                    if nCol < tool.nFactors
                        alignmentMatrices{iDS} = cat(2, alignmentMatrices{iDS}, zeros(size(alignmentMatrices{iDS}, 1), tool.nFactors - nCol));
                    end

                    % prediction is nFactors x (nTime*nConditions), though deal
                    % with tMask masking the second dimension and re-inflate
                    % along this dimension
                    prediction_masked = alignmentMatrices{iDS}' * this_dataset_centered;
                    prediction = nan(tool.nFactors, tool.nTime*tool.nConditions);
                    prediction(:, tMask) = prediction_masked;

                    tool.pcAvg_reconstructionByDataset(:, :, :, iDS) = reshape(prediction, [tool.nFactors, tool.nTime, tool.nConditions]);
                end
            end
            
            tool.alignmentMatrices = alignmentMatrices;
            tool.alignmentBiases = alignmentBiases;
        end
    
        function plotAlignmentReconstruction(tool, factorIdx, conditionIdx, varargin)
            p = inputParser();
            p.addParameter('flip', false, @islogical); % stack conditions vertically and factors horizontally instead of vice versa
            p.addParameter('datasetIdx', 1:tool.nDatasets, @isvector);
            p.parse(varargin{:});
            flip = p.Results.flip;
            
            clf;
            if isempty(tool.pcAvg_allDatasets) || isempty(tool.pcAvg_reconstructionByDataset)
                tool.computeAlignmentMatricesUsingTrialAveragedPCR();
            end
            
            if nargin < 2 || isempty(factorIdx)
                factorIdx = 1:tool.nFactors;
            elseif isscalar(factorIdx)
                factorIdx = 1:factorIdx;
            end
            
            if nargin < 3 || isempty(conditionIdx)
                conditionIdx = 1:min(tool.nConditions, 6);
            elseif isscalar(conditionIdx)
                conditionIdx = 1:conditionIdx;
            end
            
            datasetIdx = p.Results.datasetIdx;
            if islogical(datasetIdx)
                datasetIdx = find(datasetIdx);
            end
                            
            nFactorsPlot = numel(factorIdx);
            nConditionsPlot = numel(conditionIdx);

            cmap = LFADS.Utils.hslmap(numel(datasetIdx));
            
            for iF = 1:nFactorsPlot
                f = factorIdx(iF);
                for iC = 1:nConditionsPlot
                    c = conditionIdx(iC);
                    if ~iscell(tool.conditions)
                        this_cond_desc = sprintf('Condition %d', tool.conditions(c));
                    else
                        this_cond_desc = tool.conditions{c};
                        assert(ischar(this_cond_desc));
                    end
                    
                    if flip
                        iSub = sub2ind([nFactorsPlot, nConditionsPlot], iF, iC);
                        LFADS.Utils.subtightplot(nConditionsPlot, nFactorsPlot, iSub, 0.01, [0.01 0.05], [0.05 0.01]);
                    else
                        iSub = sub2ind([nConditionsPlot, nFactorsPlot], iC, iF);
                        LFADS.Utils.subtightplot(nFactorsPlot, nConditionsPlot, iSub, 0.01, [0.01 0.05], [0.05 0.01]);
                    end

                    % plot single dataset reconstructions
                    target = squeeze(tool.pcAvg_allDatasets(f, :, c));
                    data = squeeze(tool.pcAvg_reconstructionByDataset(f, :, c, datasetIdx));
                    
                    for iDS = 1:numel(datasetIdx)
                        h = plot(data(:, iDS), 'LineWidth', 0.5, 'Color', cmap(iDS, :));
                        %LFADS.Utils.setLineOpacity(h, 0.5); 
                        LFADS.Utils.showInLegend(h, tool.datasetNames{datasetIdx(iDS)});
                        hold on;
                    end
                    
                    hTarget = plot(target, 'k-', 'LineWidth', 0.5);
                    LFADS.Utils.showInLegend(hTarget, 'Global');
                    
                    hold on;

                    axis tight;
                    box off;
                    set(gca, 'YTick', [], 'XTick', []);
                    axis off;
                    
                    if (flip && iC == 1)
                        h = title(sprintf('Factor %d', f));
                        h.Visible = 'on';
                        h.FontWeight = 'normal';
                    elseif (~flip && iF == 1)
                        h = title(this_cond_desc);
                        h.Visible = 'on';
                        h.FontWeight = 'normal';
                    end
                    
                    if (flip && iF == 1)
                        h = ylabel(this_cond_desc);
                        h.Visible = 'on';
                    elseif (~flip && iC == 1)
                        h = ylabel(sprintf('Factor %d', f));
                        h.Visible = 'on';
                    end
                    if (~flip && iF == 1 && iC == nConditionsPlot) || (flip && iC ==1 && iF == nFactorsPlot)
                        legend(gca, 'show', 'Location', 'Best');
                        legend boxoff;
                    end
                end
            end
        end
    end
end