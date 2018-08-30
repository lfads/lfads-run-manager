classdef TwoStageMultisessionAlignmentTool < handle

    properties %(SetAccess=protected)
        seqData % received from run, stored to facilitate scaling to keep single trial in factors in range
        
        nFactors
        spikeBinMs % spike bins taken from run.params.spikeBinMs
        datasetNames % nDatasets x 1 cell of dataset names
        conditions % full list of unique conditions across datasets
        conditionIdxByDataset % nDatasets x vector of condition idenity
        conditionAvgsByDataset % nDatasets x 1 cell of nChannels_ds x nTimepoints x nConditions tensors of condition averages

        datasetHasSubpopulation % nDatasets x nSubpopulations logical
        subpopulationIdsByDataset % nDatasets x 1 cell of nChannels_ds x 1 vector of subpopulation ids
        numSubfactorsBySubpopulation % nSubpopulations x 1

        readinMatrices_cxs % nDatasets x 1 cell of nChannelsByDataset(iDS) x nSubFactors(iDS)
        readinMatrices_cxs_mask % nDatasets x 1 cell of nChannelsByDataset(iDS) x nFactors logical masks
        readin_subfactor_subpopulation_ids % nDatasets x 1 cell of nSubfactors(iDS) x 1 vectors indicating which subpopulation each subfactor belongs to
        biases_c % nDatasets x 1 cell of nChannelsByDataset(iDS) x 1
        readinMatrices_sxf % nDatasets x 1 cell of nSubfactors (iDS) x nFactors
        out_biases_c % nDatasets x 1 cell of nChannelsByDatasets(iDS) x 1
        
        shared_readoutMatrix_fxs
        readoutMatrices_sxc
        readoutMatrices_sxc_mask
        
        % results of prepareAlignmentMatricesUsingTrialAveragedPCR
        usingRidgeRegression logical = false;
        pcAvg_bySubpopulation % nSubpopulation cell array of nSubFactors(iS) x T x C tensors of PCs for each subpop alone including data from all sessions
        pcAvg_global % nFactors x nTimepoints x nConditions tensor of PCs using all datasets
        
        global_sigma % normalizer applied to pcAvg_global
        subpop_sigma % nSubpop normalizer applied to pcAvg_bySubpopulation
        
        recon_in_subfac % nDatasets x nSubpop cell of nSubfactors(ds, is) x T x C
        target_in_subfac % nSubpop cell of nSubfactors(ds, is) x T x C
        recon_in_fac % nDatasets cell of nFactors x T x C
        target_in_fac % nFactors x T x C
        
        recon_out_subfac % nSubpop cell of nSubfactors x T x C
        target_out_subfac % nSubpop cell of nSubfactors x T x C
        recon_out_rates % nDatasets cell of nChannels(ds) x T x C
        target_out_rates % nDatasets cell of nChannels(ds) x T x C
    end

    properties(Dependent)
        nTime
        nDatasets
        nConditions
        nSubpopulations
        nSubfactors
        
        nChannelsByDataset % nDatasets x 1 vector of channel counts
    end

    methods % Dependent properties
        function v = get.nDatasets(t)
            v = numel(t.conditionAvgsByDataset);
        end

        function v = get.nConditions(t)
            v = numel(t.conditions);
        end

        function v = get.nSubpopulations(t)
            v = numel(t.numSubfactorsBySubpopulation);
        end

        function v = get.nSubfactors(t)
            v = sum(t.numSubfactorsBySubpopulation);
        end
        
        function v = get.nTime(t)
            v = size(t.conditionAvgsByDataset{1}, 2);
        end
        
        function v = get.nChannelsByDataset(t)
            v = cellfun(@(x) size(x, 1), t.conditionAvgsByDataset);
        end
    end

    methods
        function selectDatasets(tool, mask)
            % temporary
            tool.datasetNames = tool.datasetNames(mask);
            tool.conditionIdxByDataset = tool.conditionIdxByDataset(mask);
            tool.conditionAvgsByDataset = tool.conditionAvgsByDataset(mask);
            tool.subpopulationIdsByDataset = tool.subpopulationIdsByDataset(mask);
        end
        
        function tool = TwoStageMultisessionAlignmentTool(run, seqData)
            %
            % Parameters
            % ------------
            % seqData : `nDatasets` cell of struct arrays of sequence data
            %   Sequence data for each dataset as returned by `convertDatasetToSequenceStruct`
            %

            tool.seqData = seqData;
            tool.spikeBinMs = run.params.spikeBinMs;
            tool.nFactors = run.params.c_factors_dim;
            tool.datasetNames = {run.datasets.name}';
            tool.subpopulationIdsByDataset = arrayfun(@(ds) LFADS.Utils.makecol(ds.subpopulationIds), run.datasets, 'UniformOutput', false);
            tool.numSubfactorsBySubpopulation = run.params.numSubfactorsBySubpopulation;
            nDatasets = numel(seqData);
            
            % figure out number of subpopulations from the max id in all datasets
            nSubpops = max(cat(1, tool.subpopulationIdsByDataset{:}));
            if isscalar(tool.numSubfactorsBySubpopulation)
                tool.numSubfactorsBySubpopulation = repmat(tool.numSubfactorsBySubpopulation, nSubpops, 1);
            else
                assert(numel(tool.numSubfactorsBySubpopulation) == nSubpops, 'numSubfactorsBySubpopulation must have %d subpopulations', nSubpops);
            end
            
            tool.datasetHasSubpopulation = false(nDatasets, nSubpops);
            
            for iDS = 1:nDatasets
                tool.datasetHasSubpopulation(iDS, :) = ismember(1:nSubpops, tool.subpopulationIdsByDataset{iDS});
            end
            
            % compute all unique conditions across datasets
            condField = 'conditionId';
            c = cell(nDatasets, 1);
            conditionsEachTrial = cell(nDatasets, 1);
            prog = ProgressBar(nDatasets, 'Finding unique conditions');
            for iDS = 1:nDatasets
                prog.update(iDS);
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
            if ~isnumeric(conditions) && ~iscellstr(conditions) && ~iscategorical(conditions)
                error('conditionId field must contain numeric or string values');
            end
            tool.conditions = conditions;
            nConditions = numel(conditions);
            prog.finish();
            debug('%d unique conditions found\n', nConditions);
            
            % lookup table of each trial to the unique-ified condition id list
            tool.conditionIdxByDataset = cell(nDatasets, 1);
            for iDS = 1:nDatasets
                [~, tool.conditionIdxByDataset{iDS}] = ismember(conditionsEachTrial{iDS}, conditions);
            end

            % seqData is a struct array over trials, where y contains
            % nChannels x nTimepoints
            % convert each to a nTrials x nChannels x nTimepoints tensor,
            % rebinned at run.params.spikeBinMs
            tool.conditionAvgsByDataset = cell(tool.nDatasets, 1);
            
            % populate single-trial data in countsByDataset
            prog = ProgressBar(nDatasets, 'Rebinning and condition-averaging single trial data');
            for iDS = 1:nDatasets
                prog.update(iDS);
                
                countsRebinned = tool.computeCountsRebinned(iDS);
                if iDS == 1
                    nTime = size(countsRebinned, 2);
                else
                    assert(size(countsRebinned, 2) == nTime, 'Number of timepoints (after rebinning) of all datasets must match');
                end
                
                % compute trial averaged data
                condIdx = tool.conditionIdxByDataset{iDS};
                nCh = size(countsRebinned, 1);
                tool.conditionAvgsByDataset{iDS} = nan(nCh, nTime, nConditions);

                for iC = 1:nConditions
                    tool.conditionAvgsByDataset{iDS}(:, :, iC) = nanmean(countsRebinned(:, :, condIdx == iC), 3);
                end

            end
            prog.finish();

            function a = removenan(a)
                a = a(~isnan(a));
            end

            function a = removeempty(a)
                a = a(~cellfun(@isempty, a));
            end
        end
        
        function countsRebinned = computeCountsRebinned(tool, iDS)
            % countsRebinned is nCh x T (rebinned) x R (nTrials)
            origBinMs = tool.seqData{iDS}.binWidthMs;
            rebinBy = tool.spikeBinMs / origBinMs;
            assert(abs(rebinBy - round(rebinBy)) < 1e-6, 'Ratio of new spike bin ms to original spike bin ms must be integer');
            rebinBy = round(rebinBy);
            countsOriginalBinning = cat(3, tool.seqData{iDS}.y); % Ch x T x tRials

            % rebin to tool.spikeBinMs
            N = size(countsOriginalBinning, 1);
            R = size(countsOriginalBinning, 3);
            Ttruncated = floor(size(countsOriginalBinning, 2) / rebinBy) * rebinBy;

            countsRebinned = reshape(sum(reshape(...
                countsOriginalBinning(:, 1:Ttruncated, :), ...
                [N rebinBy Ttruncated/rebinBy R]), 2), [N Ttruncated/rebinBy R]);
        end

        function computeAlignmentMatrices(tool, varargin)
            % Theis implementation computes trial-averages (averaging
            % all trials with the same conditionId label) for each neuron
            % in each session. The trial-averages are then assembled into a
            % large nNeuronsTotal x (nConditions x time) matrix. The top
            % nFactors PCs of this matrix are computed (as linear
            % combinations of neurons). For each session, we then regress
            % the nNeuronsThisSession neurons against the top nFactors PCs.
            % The alignment matrix is the matrix of regression
            % coefficients.
 
            % Need to compute:
            % readinMatrices_cxs % nDatasets x 1 cell of nChannelsByDataset(iDS) x nSubFactors(iDS)
            % readinMatrices_cxs_mask % nDatasets x 1 cell of nChannelsByDataset(iDS) x nFactors logical masks
            % biases_c % nDatasets x 1 cell of nChannelsByDataset(iDS) x 1
            % readinMatrices_sxf % nDatasets x 1 cell of nSubfactorsUsed(iDS) x nFactors
            % shared_readoutMatrix_fxs % nFactors x nSubfactorsTotal
            % readoutMatrices_sxc % nDatasets x 1 cell of nSubfactorsTotal x nChannelsByDataset(iDS)
           
           function [pc_data_tensor, sigma] = do_pca_tensor_normalize(data_tensor, nPCs)
                % do pca along dim 1 of tensor and return projected scores as tensor with same shape along dims 2:3
                data_mat = data_tensor(:, :);
                %data_centered = bsxfun(@minus, data_mat, nanmean(data_mat, 2));

                try
                    [~, pc_mat, ~] = pca(data_mat', 'Rows', 'pairwise', 'NumComponents', nPCs, 'Centered', true);
                catch
                    [~, pc_mat] = pca(data_mat', 'Rows', 'complete', 'NumComponents', nPCs, 'Centered', true);
                end
                
                % normalize globally such that first pc has unit variance
                sigma = std(pc_mat(:, 1));
                if sigma == 0, sigma = 1; end
                pc_data_tensor = reshape(pc_mat' ./ sigma, [nPCs, size(data_tensor, 2), size(data_tensor, 3)]);
            end
            
            function [W, Y_hat] = do_regression(X_tensor, Y_tensor)
                % regress data in Y_tensor on data in X_tensor where the linear combination is done on dim 1 of X
                % W will be size(X, 1) x size(Y, 1)
                % alignmentMatrix is nChannelsByDataset(iDS) x nFactors
                % this_dataset_centered' is (TC x nChannels)
                % dim_reduced_data_this' is (TC x nFactors
                
                X_mat = X_tensor(:, :)';
                Y_mat = Y_tensor(:, :)';
                W = LFADS.Utils.ridge_cv(Y_mat, X_mat, ...
                    'normalizeEach', false, 'KFold', 10);
                
                if nargout > 1
                    Y_hat = do_prediction(X_tensor, W);
                end
            end
            
            function Y_hat = do_prediction(X_tensor, W)
                % X is channels x time x conditions
                % W is channels x K
                % Y_hat is K x time x conditions
                
                Y_hat = TensorUtils.linearCombinationAlongDimension(X_tensor, 1, W');
            end
                
            tool.usingRidgeRegression = true; % altered below if we end up actually aligning
            
            %% 1: COMPUTE GLOBAL PCS
            
            all_data_tensor = cat(1, tool.conditionAvgsByDataset{:}); % nChannelsTotal x nTime x nConditions
            [~, T, C] = size(all_data_tensor);
            [global_pc_data, tool.global_sigma] = do_pca_tensor_normalize(all_data_tensor, tool.nFactors);
            
            % nPCs x time x conditions
            tool.pcAvg_global = global_pc_data;
            
            all_means = nanmean(all_data_tensor(:, :), 2);
            all_data_centered = all_data_tensor - all_means;
            
            %% 2: COMPUTE SUBPOPULATION PCS
            
            subpopulationIds = cat(1, tool.subpopulationIdsByDataset{:});
            nSubpop = tool.nSubpopulations;
            
            [subpop_pc_data, tool.target_in_subfac, tool.target_out_subfac] = deal(cell(nSubpop, 1));
            tool.subpop_sigma = nan(nSubpop, 1);
            
            prog = ProgressBar(nSubpop, 'Computing subpopulation target PCs');
            for iS = 1:nSubpop
                prog.update(iS);
                
                % gather data from neurons within this subpopulation
                mask = subpopulationIds == iS;
                if ~any(mask)
                    warning('No neurons from subpopulation id %d found', iS);
                    subpop_pc_data{iS} = zeros([tool.numSubfactorsBySubpopulation(iS), size(subpop_data_centered, 2), size(subpop_data_centered, 3)]);
                    continue;
                end
                subpop_data_centered = all_data_centered(mask, :, :);
                
                % check there are enough neurons as subfactors needed across all datasets
                if tool.numSubfactorsBySubpopulation(iS) > nnz(mask)
                    % not enough, do pca and pad with small random values
                    warning('More subfactors for subpopulation %d than neurons across all datasets', iS);
                    nFound = nnz(mask);
                    szCat = [tool.numSubfactorsBySubpopulation(iS) - nFound, size(subpop_data_centered, 2), size(subpop_data_centered, 3)];
                    [partial, tool.subpop_sigma(iS)] = do_pca_tensor_normalize(subpop_data_centered, nFound);
                    subpop_pc_data{iS} = cat(1, partial, 0.1*randn(szCat));
                else
                    [subpop_pc_data{iS}, tool.subpop_sigma(iS)] = do_pca_tensor_normalize(subpop_data_centered, tool.numSubfactorsBySubpopulation(iS));
                end
                
                tool.target_in_subfac{iS} = subpop_pc_data{iS};
                tool.target_out_subfac{iS} = subpop_pc_data{iS};
            end
            prog.finish();
            
            tool.pcAvg_bySubpopulation = subpop_pc_data;
            
            %% 3: SESSION-SPECIFIC RATES TO SUBFACTORS READIN MATRICES, BIASES
            %     SESSION-SPECIFIC ALL-SUBFACTORS TO RATES READOUT MATRICES
            %     SESSION-SPECIFIC SUBFACTORS TO FACTORS READIN MATRICES
            
            prog = ProgressBar(tool.nDatasets, 'Computing readin / readout matrices by dataset');
            
            [readinMatrices_cxs, readinMatrices_cxs_mask, biases_c, readin_subfactor_subpopulation_ids, ...
                readinMatrices_sxf, readoutMatrices_sxc, readoutMatrices_sxc_mask, out_biases_c] = deal(cell(tool.nDatasets, 1));
            [tool.recon_in_subfac] = deal(cell(tool.nDatasets, nSubpop));
            tool.recon_in_fac = cell(tool.nDatasets, 1);
            
            extremeValueSingleTrial = nan(tool.nDatasets, 1);
            for iDS = 1:tool.nDatasets
                prog.update(iDS);
                % center the rates and use the means as the bias
                rates = tool.conditionAvgsByDataset{iDS};
                biases_c{iDS} = nanmean(rates(:, :), 2); %#ok<*PROPLC>
                rates_centered = rates - biases_c{iDS};
                
                subPopIds = tool.subpopulationIdsByDataset{iDS};
                subpopsPresent = ismember(1:tool.nSubpopulations, unique(subPopIds));
                
                nRates = size(tool.conditionAvgsByDataset{iDS}, 1);
                                
                [readin_each_subpop, readin_each_subpop_mask, which_subpop] = deal(cell(tool.nSubpopulations, 1));
                for iS = 1:tool.nSubpopulations
                    if ~subpopsPresent(iS)
                        tool.recon_in_subfac{iDS, iS} = zeros(0, T, C);
                       	continue;
                    end
                    % do the regression of subpop pcs onto the rates
                    [temp_only_neurons_from_subpop, tool.recon_in_subfac{iDS, iS}] = ...
                        do_regression(rates_centered(subPopIds == iS, :, :), subpop_pc_data{iS});
                    
                    % assign into the appropriate rows (neurons from that subpopulation) in the matrix
                    readin_each_subpop{iS} = zeros(nRates, tool.numSubfactorsBySubpopulation(iS));
                    readin_each_subpop{iS}(subPopIds == iS, :) = temp_only_neurons_from_subpop;
                    
                    % assign true the appropriate rows (neurons from that subpopulation) in the mask
                    readin_each_subpop_mask{iS} = false(nRates, tool.numSubfactorsBySubpopulation(iS));
                    readin_each_subpop_mask{iS}(subPopIds == iS, :) = true;
                    
                    % identify which subpopulation is associated with each subfacotr
                    which_subpop{iS} = repmat(iS, tool.numSubfactorsBySubpopulation(iS), 1);
                end
                
                % concatenate over subpops to yield one c x s matrix and mask
                readinMatrices_cxs{iDS} = cat(2, readin_each_subpop{subpopsPresent});
                readinMatrices_cxs_mask{iDS} = cat(2, readin_each_subpop_mask{subpopsPresent});
                readin_subfactor_subpopulation_ids{iDS} = cat(1, which_subpop{subpopsPresent});
                               
                % compute matrix to predict in factors (global pcs) from
                % subfactors present
                
                % old version: use target subpop pcs as the input
                % subpop_present_pc_data = cat(1, subpop_pc_data{subpopsPresent});

                % new version: use the actual reconstruction (rates * W_cxs)                
                subpop_present_pc_data = cat(1, tool.recon_in_subfac{iDS, subpopsPresent});
                
                [readinMatrices_sxf{iDS}, tool.recon_in_fac{iDS}] = do_regression(subpop_present_pc_data, global_pc_data);

                % SCALING:
                % Important: the condition averaged rates projected into subfactors are reasonably scaled now.
                % but the single trials need not be, and this breaks things. We need to scale down this cxs matrices
                % to get the maximum single trial deviation between -1 and 1
                singleTrial = tool.computeCountsRebinned(iDS); % nCh x T x tRial
                singleTrial = singleTrial - biases_c{iDS};
                Wcxf = readinMatrices_cxs{iDS} * readinMatrices_sxf{iDS};
                singleTrial = TensorUtils.linearCombinationAlongDimension(singleTrial, 1, Wcxf');
                extremeValueSingleTrial(iDS) = quantile(abs(singleTrial(:)), 0.999);
                
                tool.target_in_fac = global_pc_data;
                
                % compute readout matrices to predict log rates from
                % subfactors included on that dataset
                [readout_each_subpop, readout_each_subpop_mask, which_subpop] = deal(cell(tool.nSubpopulations, 1)); 
                
                % for readout:
                % rates = exp(subfac * W_sxc + biases) where the bias is inside the exponentiation
                % so we want to regress:
                % log_rates = subfac * W_sxc + b for W and b
                % but there are many rates == 0, so instead of taking log,
                % we linearize log(rates) ~= rates - 1 which is fine
                % because rates is in spikes / bin, not spikes / second, so
                % it is generally small
                % rates - 1 = subfac * W_sxc + b
                log_rates = rates - 1; % this is approximate, linearization around rates = 1
                log_biases_c = nanmean(log_rates(:, :), 2);
                log_rates_centered = log_rates - log_biases_c;
                
                for iS = 1:tool.nSubpopulations
                    if ~subpopsPresent(iS)
                        % readout matrix will multiply an s x 1 vector of subfactors, so 
                        % so we pad missing subpops with zeros (so they
                        % don't contribute)
                        readout_each_subpop{iS} = zeros(tool.numSubfactorsBySubpopulation(iS), nRates);
                        readout_each_subpop_mask{iS} = false(tool.numSubfactorsBySubpopulation(iS), nRates);
                        tool.recon_out_rates{iDS, iS} = zeros(0, T, C);
                        tool.target_out_rates{iDS, iS} = zeros(0, T, C);
                    else
                        % do the regression of log rates onto subpop pcs=
                        [temp_only_neurons_from_subpop, tool.recon_out_rates{iDS, iS}] = ...
                            do_regression(subpop_pc_data{iS}, log_rates_centered(subPopIds == iS, :, :));
                        tool.target_out_rates{iDS, iS} = log_rates_centered(subPopIds == iS, :, :);
                        
                        % assign into the appropriate rows (neurons from that subpopulation) in the matrix
                        readout_each_subpop{iS} = zeros(tool.numSubfactorsBySubpopulation(iS), nRates);
                        readout_each_subpop{iS}(:, subPopIds == iS) = temp_only_neurons_from_subpop;
                    
                        % assign true the appropriate columns (neurons from that subpopulation) in the mask
                        readout_each_subpop_mask{iS} = false(tool.numSubfactorsBySubpopulation(iS), nRates);
                        readout_each_subpop_mask{iS}(:, subPopIds == iS) = true;
                    
                        % identify which subpopulation is associated with each subfacotr
                        which_subpop{iS} = repmat(iS, tool.numSubfactorsBySubpopulation(iS), 1);
                    end
                end
                         
                % concatenate over subpops to yield one nSubpops(iDS) x
                % num_channels 
                readoutMatrices_sxc{iDS} = cat(1, readout_each_subpop{:});
                readoutMatrices_sxc_mask{iDS} = cat(1, readout_each_subpop_mask{:});
                out_biases_c{iDS} = cat(1, log_biases_c);
            end
            prog.finish();
            
            % SCALING: we have computed the extreme values for each dataset on single trials,
            % now we pick the max and scale down all matrices uniformly
            % we don't bother to rescale the recon and target fields here, though we could
            scaleDownBy = max(extremeValueSingleTrial);
            for iDS = 1:tool.nDatasets
                readinMatrices_cxs{iDS} = readinMatrices_cxs{iDS} ./ sqrt(scaleDownBy);
                readinMatrices_sxf{iDS} = readinMatrices_sxf{iDS} ./ sqrt(scaleDownBy);
            end

            %% 4: SHARED READOUT FROM FACTORS TO SUBFACTORS
            % map from factors to subfactor pcs
            all_subpop_pc_data = cat(1, subpop_pc_data{:});
            [shared_readoutMatrix_fxs, recon_out_subfac_cat] = do_regression(global_pc_data, all_subpop_pc_data);
            
            % split into subpopulations
            tool.recon_out_subfac = TensorUtils.splitAlongDimension(recon_out_subfac_cat, 1, tool.numSubfactorsBySubpopulation);
            tool.target_out_subfac = subpop_pc_data;
            
            %% 5: ASSIGN EVERYTHING COMPUTED
            tool.readinMatrices_cxs = readinMatrices_cxs;
            tool.readinMatrices_cxs_mask = readinMatrices_cxs_mask;
            tool.readin_subfactor_subpopulation_ids = readin_subfactor_subpopulation_ids;
            tool.biases_c = biases_c;
            tool.readinMatrices_sxf = readinMatrices_sxf;

            tool.shared_readoutMatrix_fxs = shared_readoutMatrix_fxs;
            tool.readoutMatrices_sxc = readoutMatrices_sxc;
            tool.readoutMatrices_sxc_mask = readoutMatrices_sxc_mask;
            tool.out_biases_c = out_biases_c;
        end
        
        function plot_recon_in_subfac(tool, iSub, varargin)
            reconByDataset = tool.recon_in_subfac(:, iSub);
            target = tool.target_in_subfac{iSub};
            tool.internal_plotReconVsTarget(reconByDataset, target, varargin{:}, 'basisName', 'In Subfac')
        end

        function plot_recon_in_fac(tool, varargin)
            reconByDataset = tool.recon_in_fac;
            target = tool.target_in_fac;
            tool.internal_plotReconVsTarget(reconByDataset, target, varargin{:}, 'basisName', 'In Factor')
        end

        function plot_recon_out_subfac(tool, iSub, varargin)
            recon = tool.recon_out_subfac{iSub};
            target = tool.target_out_subfac{iSub};
            tool.internal_plotReconVsTarget(recon, target, varargin{:}, 'basisName', 'Out Subfac')
        end
        
        function plot_recon_out_rates(tool, iDataset, iSub, varargin)
            recon = tool.recon_out_rates{iDataset, iSub};
            target = tool.target_out_rates{iDataset, iSub};
            tool.internal_plotReconVsTarget(recon, target, varargin{:}, 'basisName', 'Rate')
        end
        
        function internal_plotReconVsTarget(tool, reconByDataset, target, varargin)
            % plotAlignmentReconstruction(tool, basisIdx, conditionIdx, varargin)
            
            if iscell(reconByDataset)
                nonEmpty = ~cellfun(@isempty, reconByDataset);
                if ~any(nonEmpty)
                    error('No datasets have reconstruction data');
                end
                reconTensor = reconByDataset{find(nonEmpty, 1)};
            else
                reconTensor = reconByDataset;
            end
            nBases = size(reconTensor, 1);
            nConditions = size(reconTensor, 3);

            p = inputParser();
            p.addOptional('basisIdx', 1:min(6, nBases), @(x) isnumeric(x) && isvector(x));
            p.addOptional('conditionIdx', 1:min(6, nConditions), @(x) isnumeric(x) && isvector(x));
            p.addParameter('basisName', 'Basis', @ischar);
            p.addParameter('flip', false, @islogical); % stack conditions vertically and factors horizontally instead of vice versa
            p.addParameter('datasetIdx', 1:tool.nDatasets, @isvector);
            p.parse(varargin{:});
            flip = p.Results.flip;

            clf;
            
            basisIdx = p.Results.basisIdx;
            conditionIdx = p.Results.conditionIdx;

            datasetIdx = p.Results.datasetIdx;
            if islogical(datasetIdx)
                datasetIdx = find(datasetIdx);
            end

            nBasesPlot = numel(basisIdx);
            nConditionsPlot = numel(conditionIdx);

            cmap = LFADS.Utils.hslmap(numel(datasetIdx));

            for iB = 1:nBasesPlot
                b = basisIdx(iB);
                for iC = 1:nConditionsPlot
                    c = conditionIdx(iC);
                    if ~iscell(tool.conditions)
                        this_cond_desc = sprintf('Condition %d', tool.conditions(c));
                    else
                        this_cond_desc = tool.conditions{c};
                        assert(ischar(this_cond_desc));
                    end

                    if flip
                        iSub = sub2ind([nBasesPlot, nConditionsPlot], iB, iC);
                        LFADS.Utils.subtightplot(nConditionsPlot, nBasesPlot, iSub, 0.01, [0.01 0.05], [0.05 0.01]);
                    else
                        iSub = sub2ind([nConditionsPlot, nBasesPlot], iC, iB);
                        LFADS.Utils.subtightplot(nBasesPlot, nConditionsPlot, iSub, 0.01, [0.01 0.05], [0.05 0.01]);
                    end

                    % plot single dataset reconstructions
                    target_trace = squeeze(target(b, :, c));
                    
                    if iscell(reconByDataset)
                        for iDS = 1:numel(datasetIdx)
                            if isempty(reconByDataset{iDS}), continue; end
                            recon_ds = squeeze(reconByDataset{iDS}(b, :, c));
                            h = plot(recon_ds, 'LineWidth', 0.5, 'Color', cmap(iDS, :));
                            %LFADS.Utils.setLineOpacity(h, 0.5);
                            LFADS.Utils.showInLegend(h, tool.datasetNames{datasetIdx(iDS)});
                            hold on;
                        end
                    else
                        recon_this = squeeze(reconTensor(b, :, c));
                        h = plot(recon_this, 'LineWidth', 0.5, 'Color', cmap(1, :));
                        %LFADS.Utils.setLineOpacity(h, 0.5);
                        LFADS.Utils.showInLegend(h, 'recon');
                        hold on;
                    end

                    hTarget = plot(target_trace, 'k--', 'LineWidth', 0.5);
                    LFADS.Utils.showInLegend(hTarget, 'Target');

                    hold on;

                    axis tight;
                    box off;
                    set(gca, 'YTick', [], 'XTick', []);
                    axis off;

                    if (flip && iC == 1)
                        h = title(sprintf('%s %d', p.Results.basisName, b));
                        h.Visible = 'on';
                        h.FontWeight = 'normal';
                    elseif (~flip && iB == 1)
                        h = title(this_cond_desc);
                        h.Visible = 'on';
                        h.FontWeight = 'normal';
                    end

                    if (flip && iB == 1)
                        h = ylabel(this_cond_desc);
                        h.Visible = 'on';
                    elseif (~flip && iC == 1)
                        h = ylabel(sprintf('%s %d', p.Results.basisName, b));
                        h.Visible = 'on';
                    end
                    if (~flip && iB == 1 && iC == nConditionsPlot) || (flip && iC ==1 && iB == nBasesPlot)
                        legend(gca, 'show', 'Location', 'Best');
                        legend boxoff;
                    end
                end
            end
        end
    end
end
