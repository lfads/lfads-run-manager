classdef Run < handle & matlab.mixin.CustomDisplay
    % Represents a single LFADS experiment on a specific set of datasets. Runs are grouped into :ref:`LFADS_RunCollection`
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.
    
    methods(Abstract)
        % These methods will need to be implemented in a subclass that provides the custom behavior for your application
        
        seq = convertDatasetToSequenceStruct(r, dataset, mode, varargin);
        % Converts the loaded data within a dataset into a sequence struct. The sequence data returned will be a
        % struct array where each element corresponds to a trial. You can include any metadata or information fields
        % that you like for future reference or analysis. At a minimum, you must include field `.y`, `.y_time`, and `.params.dtMS`
        %
        % For each trial:
        %   - `.y` will contain spike data as an `nNeurons` x `nTimeBins` array of spike counts. The spike binning is
        %       for this is determined by you, and can be left at 1 ms. Later,
        %       this data will be rebinned according to RunParams .spikeBinMs field.
        %   - `.y_time` provides a time vector corresponding to the time
        %       bins of `.y`, which should be identical on each trial
        %   - `.binWidthMs` specifies the time bin width for `.y`
        %   - optionally: `.conditionId` specifies the condition to which each trial 
        %       belongs. This information isn't passed to LFADS. It is used
        %       only when building the alignment matrices for multi-session
        %       stitching, if trial-averaging is employed.
        %
        % Parameters
        % ------------
        % dataset : :ref:`LFADS_Dataset`
        %   The :ref:`LFADS_Dataset` instance from which data were loaded
        % mode (string) : typically 'export' indicating sequence struct
        %   will be exported for LFADS, or 'alignment' indicating that this
        %   struct will be used to generate alignment matrices
        %
        % Returns
        % ----------
        % seq : struct Array
        %   sequence formatted data. A struct array where each elemnt corresponds to a specific trial.
        
    end
    
    properties
        name char = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk
        
        comment char = '' % Textual comment for convenience
        
        version uint32 = 4; % Internal versioning allowing for graceful evolution of path settings
    end
    
    properties
        runCollection % :ref:`LFADS_RunCollection` instance to which this run belongs
        params % :ref:`LFADS_RunParams` instance shared by all runs in the collection, contains parameter settings
        paramIndexInRunCollection % which RunParams index in the RunCollection this params corresponds to
        datasets % Array of :ref:`LFADS_Dataset` instances which this particular Run will utilize
        
        sequenceData cell % nDatasets cell array of sequence struct data
        posteriorMeans % nDatasets array of :ref:`LFADS_PosteriorMeans` when loaded

        inputInfo % parameters used to train model
    end
    
    properties(Dependent)
        nDatasets % Number of datasets used by this run
        datasetCollection % Dataset collection used by this run (and all runs in the same RunCollection)
        path % Unique folder within rootPath including paramStr/name
        
        datasetIndsInCollection % indices of each dataset into datasetCollection.datasets
        
        paramsString % string representation of params generated using .params.generateString()
        
        pathCommonData % Path on disk where original data files are saved, shared by all Runs in this collection
        
        pathSequenceFiles % Path on disk where sequence files may be saved
        sequenceFileNames % List of sequence file names (sans path)
        
        pathLFADSInput % Path on disk where LFADS input hd5 files for this run will be symlinked into, from their true location in pathCommonData
        lfadsInputFileNames % List of LFADS input hd5 files (sans path)

        lfadsInputInfoFileNames % List of LFADS input info .mat files (sans path)
        
        pathLFADSOutput % Path on disk where LFADS output will be written
        
        fileShellScriptLFADSTrain % Location on disk where shell script for LFADS training will be written
        fileShellScriptLFADSPosteriorMeanSample  % Location on disk where shell script for LFADS posterior mean sampling will be written
        fileShellScriptLFADSWriteModelParams % Location on disk where shell script for LFADS write model params will be written
        
        fileLFADSOutput % output from training and sampling can be tee'd here
        
        fileModelParams % Location on disk where model params will be written
        
        sessionNameTrain % name of tmux session that will be created if useSession = true is passed to writeShellScriptLFADSTrain
        sessionNamePosteriorMean % name of tmux session that will be created if useTmuxSession = true is passed to writeShellScriptLFADSPosteriorMean
        sessionNameWriteModelParams % name of tmux session that will be created if useTmuxSession = true is passed to writeShellScriptLFADSWriteModelParams
    end
    
    methods
        function r = Run(varargin)

        end
        
        function tf = usesDifferentSequenceDataForAlignment(r) 
            % tf = usesDifferentSequenceDataForAlignment() 
            %
            % Returns true if you would like the Run to call your
            % convertDatasetToSequenceStruct with a mode == 'alignment'
            % argument when constructing the alignment matrices. This
            % allows you to specify a different set of trials used for
            % constructing alignment matrices, e.g. only correct trials.
            
            tf = false;
        end
        
        function alignmentMatrices = prepareAlignmentMatrices(r, seqData)
            % Prepares alignment matrices to seed the stitching process when 
            % using multiple days of sequence data for LFADS input file generation. 
            % Generate alignment matrices which specify the initial guess at the 
            % encoder matrix that converts neural activity from each dataset 
            % to a common set of factors (for stitching). Each alignment matrix
            % should be nNeurons (that session) x nFactors.
            % 
            % The default implementation computes trial-averages (averaging
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
            
            % tally up the total number of channels across all datasets
            num_channels = 0;
            
            % gather global unique conditions
            condField = 'conditionId';
            c = cell(r.nDatasets, 1);
            conditionsEachTrial = cell(r.nDatasets, 1);
            for nd = 1:r.nDatasets
               conditionsEachTrial{nd} = {seqData{nd}.(condField)}';
               if isscalar(conditionsEachTrial{nd}{1})
                   conditionsEachTrial{nd} = cell2mat(conditionsEachTrial{nd});
               end
                   
               c{nd} = unique(conditionsEachTrial{nd});
            end
            conditions = unique(cat(1, c{:}));
            if isnumeric(conditions)
                conditions = conditions(~isnan(conditions));
            elseif iscellstr(conditions)
                conditions = conditions(~cellfun(@isempty, conditions));
            else
                error('conditionId field must contain numeric or string values');
            end
            
            % split each struct by condition
            for nd = 1:r.nDatasets
                datasetInfo(nd).seq = seqData{nd}; %#ok<AGROW>
                
                % we are going to make a big array that spans all days.
                % store down the indices for this day
                this_day_num_channels = size(datasetInfo(nd).seq(1).y,1);
                datasetInfo(nd).this_day_inds = num_channels + (1 : this_day_num_channels); %#ok<AGROW>
                num_channels = num_channels + this_day_num_channels; 
            end
            
            % we are going to make a matrix that is
            % num_channels x (time_per_trial x total_num_trials)
            bin_size = r.params.spikeBinMs;
            time_per_trial = floor(size(datasetInfo(1).seq(1).y, 2) / ...
                bin_size);
            all_data = nan(num_channels, time_per_trial * numel(conditions), ...
                'single');
            
            % fill up the data matrix with binned data
            for nd = 1:numel(datasetInfo)
                [~, ic] = ismember(conditionsEachTrial{nd}, conditions);
                
                % start at the zero time point for each day
                data_time_ind = 0;
                for ncond = 1:numel(conditions)
                    trials_to_use_this_condition = find(ic==ncond);
                    n_trials_this_condition = numel(trials_to_use_this_condition);
                    binned_data = nan(numel(datasetInfo(nd).this_day_inds),time_per_trial, n_trials_this_condition);
                    for nt = 1:n_trials_this_condition
                        trial = datasetInfo(nd).seq(trials_to_use_this_condition(nt));
                        tmp = reshape(trial.y(:, 1:(bin_size*time_per_trial))', ...
                            [], time_per_trial, size(trial.y,1));
                        binned_data(:, :, nt) = squeeze(sum(tmp))';                  
                    end
                    
                    % average over trials in the condition
                    binned_data = nanmean(binned_data, 3);
                    if any(isnan(binned_data(:))) && n_trials_this_condition > 0
                        binned_data(isnan(binned_data)) = 0;
                        warning('NaN data found on condition %d dataset %d', ncond, nd);
                    end
                    all_data(datasetInfo(nd).this_day_inds, data_time_ind ...
                        + (1:time_per_trial)) = binned_data;
                    data_time_ind = data_time_ind + time_per_trial;
                end
            end
            
            % apply PCA
            try
                keep_pcs = pca(all_data', 'Rows', 'pairwise', 'NumComponents', r.params.c_in_factors_dim);
            catch
                keep_pcs = pca(all_data', 'Rows', 'complete', 'NumComponents', r.params.c_in_factors_dim);
            end
            
            % project all data into pca space
            all_data_centered = bsxfun(@minus, all_data, nanmean(all_data, 2));
            dim_reduced_data = keep_pcs' * all_data_centered;
            
            % get a mapping from each day to the lowD space
            [this_day_data, this_data_predicted] = cellvec(numel(datasetInfo)); % each is nChannels x time x 
            for nd = 1:numel(datasetInfo)
                this_day_data{nd} = all_data(datasetInfo(nd).this_day_inds, :);
                
                % figure out which timepoints are valid
                tMask = ~any(isnan(this_day_data{nd}), 1) & ~any(isnan(dim_reduced_data), 1);
                dim_reduced_data_this = bsxfun(@minus, dim_reduced_data(:, tMask), ...
                    nanmean(dim_reduced_data(:, tMask), 2));
                this_data_centered = bsxfun(@minus, this_day_data{nd}(:, tMask), ...
                    nanmean(this_day_data{nd}(:, tMask), 2));
                
                % TEMPORARY to deal with bias, TODO remove this when LFADS
                % has bias term
                if false
                    % regress this day's data against the global PCs
                    datasetInfo(nd).alignment_matrix_cxf = (this_data_centered' \ dim_reduced_data_this');
                else
                    % regress uncentered data onto the centered PCs,
                    % assumes last row of this_day_data is all ones so it
                    % can provide the bias that cancel's the mean of the
                    % data on this day
                    this_data_notCentered = this_day_data{nd}(:, tMask);
                    datasetInfo(nd).alignment_matrix_cxf = (this_data_notCentered' \ dim_reduced_data_this');
                end
                if any(isnan(datasetInfo(nd).alignment_matrix_cxf(:)))
                    error('NaNs in the the alignment matrix');
                end
                    
                this_data_predicted{nd} = datasetInfo(nd).alignment_matrix_cxf' * this_data_centered;
            end
            
            alignmentMatrices = {datasetInfo.alignment_matrix_cxf};
        end
        
        function tf = eq(a, b)
            % Overloaded == operator to enable equality if name, params,
            % datasets, and runCollection fields all match.
            
            tf = false(size(a));
            assert(isequal(size(b), size(a)), 'Sizes must match');
            for i = 1:numel(a)
                tf(i) = strcmp(a(i).name, b(i).name) && isequal(a(i).params, b(i).params) ...
                    && isequal(a(i).datasets, b(i).datasets) && isequal(a(i).runCollection, b(i).runCollection);
            end
        end
        
        function dc = get.datasetCollection(r)
            dc = r.runCollection.datasetCollection;
        end
        
        function p = get.path(r)
            if isempty(r.runCollection)
                p = '';
            elseif r.version < 3
                % collectionPath_name
                p = fullfile(r.runCollection.path, r.name);
            else
                % collectionPath / paramString / name
                p = fullfile(r.runCollection.path, r.paramsString, r.name);
            end
        end
        
        function str = get.paramsString(r)
            if ~isempty(r.params)
                str = r.params.generateHashName();
            else
                str = '';
            end
        end
        
        function p = get.pathCommonData(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.runCollection.path, r.params.generateInputDataHashName());
            end
        end
        
        function p = get.pathSequenceFiles(r)
            if isempty(r.runCollection)
                p = '';
            elseif r.version < 4
                p = fullfile(r.path, 'seq');
            else
                p = fullfile(r.pathCommonData, 'seq');
            end
        end
        
        function names = get.lfadsInputFileNames(r)
            if r.version < 2
                names = arrayfun(@(ds) sprintf('lfads_%s_spikes.h5',  [r.nameWithParams '_' ds.name]), ...
                    r.datasets, 'UniformOutput', false);
            else
                % remove redundant info from name already in the path
                names = arrayfun(@(ds) sprintf('lfads_%s.h5',  ds.name), ...
                    r.datasets, 'UniformOutput', false);
            end
        end
        
        function names = get.lfadsInputInfoFileNames(r)
            names = arrayfun(@(ds) sprintf('inputInfo_%s.mat',  ds.name), ...
                    r.datasets, 'UniformOutput', false);
        end

        function p = get.pathLFADSInput(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.path, 'lfadsInput');
            end
        end
        
        function p = get.pathLFADSOutput(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.path, 'lfadsOutput');
            end
        end
        
        function f = get.fileShellScriptLFADSTrain(r)
            f = fullfile(r.path, 'lfads_train.sh');
        end
        
        function f = get.fileShellScriptLFADSPosteriorMeanSample(r)
            f = fullfile(r.path, 'lfads_posterior_mean_sample.sh');
        end
        
        function f = get.fileShellScriptLFADSWriteModelParams(r)
            f = fullfile(r.path, 'lfads_write_model_params.sh');
        end
        
        function f = get.fileModelParams(r)
            f = fullfile(r.pathLFADSOutput, 'model_params');
        end
        
        function f = get.fileLFADSOutput(r)
            f = fullfile(r.path, 'lfads.out');
        end
        
        function n = get.nDatasets(r)
            n = numel(r.datasets);
        end
        
        function names = get.sequenceFileNames(r)
            if r.version < 2
                names = arrayfun(@(ds) [r.nameWithParams '_' ds.name '_seq.mat'], r.datasets, 'UniformOutput', false);
            else
                % simpler file names without extra redundant info already
                % in path
                names = arrayfun(@(ds) ['seq_' ds.name '.mat'], r.datasets, 'UniformOutput', false);
            end
        end
        
        function [trainList, validList] = getLFADSPosteriorSampleMeanFiles(r)
            % Generates the list of training and validation LFADS posterior mean files for loading, without path
            %
            % Returns
            % ---------
            % trainList : cellstr
            %   list of file names for training posterior samples
            % validList : cellstr
            %   list of file names for validation posterior samples
            
            if r.version < 2
                trainList = arrayfun(@(ds) sprintf('model_runs_%s_spikes.h5_train_posterior_sample',  [r.nameWithParams '_' ds.name]), ...
                    r.datasets, 'UniformOutput', false);
                validList = arrayfun(@(ds) sprintf('model_runs_%s_spikes.h5_valid_posterior_sample',  [r.nameWithParams '_' ds.name]), ...
                    r.datasets, 'UniformOutput', false);
            else
                trainList = arrayfun(@(ds) sprintf('model_runs_%s.h5_train_posterior_sample', ds.name), ...
                    r.datasets, 'UniformOutput', false);
                validList = arrayfun(@(ds) sprintf('model_runs_%s.h5_valid_posterior_sample', ds.name), ...
                    r.datasets, 'UniformOutput', false);
            end
        end
        
        function sess = get.sessionNameTrain(r)
            sess = sprintf('train_%s_%s', r.name, r.paramsString);
            sess = strrep(sess, '.', '_');
        end
        
        function sess = get.sessionNamePosteriorMean(r)
            sess = sprintf('pm_%s_%s', r.name, r.paramsString);
            sess = strrep(sess, '.', '_');
        end
        
        function idx = get.datasetIndsInCollection(r)
            [~, idx] = ismember(r.datasets, r.datasetCollection.datasets);
        end
    end
    
    methods(Hidden)
        function h = getFirstLineHeader(r)
            className = class(r);
            h = sprintf('%s "%s" (%d datasets)\n', className, r.name, r.nDatasets);
        end
    end
    
    methods (Access = protected)
        function header = getHeader(r)
            if ~isscalar(r)
                header = getHeader@matlab.mixin.CustomDisplay(r);
            else
                rc = r.runCollection;
                header = sprintf('%s\n  Path: %s\n  Data: %s\n  %s "%s" : %s\n\n  %d datasets in "%s"\n', ...
                    r.getFirstLineHeader(), r.path, r.pathCommonData, ...
                    class(r.params), r.paramsString, r.params.generateShortDifferencesString, ...
                    r.nDatasets, r.datasetCollection.name);
                for s = 1:r.nDatasets
                    header = cat(2, header, sprintf('    [%2d] %s', s, r.datasets(s).getHeader()));
                end
            end
        end
    end
    
    methods
        function prepareForLFADS(r, regenerate)
            % Generate all files needed to run LFADS.
            if nargin < 2
                regenerate = false;
            end
            r.makeLFADSInput(regenerate);
        end
        
        function deleteLFADSOutput(r, varargin)
            p = inputParser();
            p.addParameter('confirmed', false, @islogical);
            p.parse(varargin{:});
            
            if ~p.Results.confirmed
                resp = input(sprintf('Are you sure you want to delete %s:\n', r.pathLFADSOutput), 's');
                if lower(resp(1)) ~= 'y'
                    return;
                end
            end
            
            fprintf('Deleting %s\n', r.pathLFADSOutput);
            cmd = sprintf('rm --preserve-root -rf "%s"', r.pathLFADSOutput);
            [s, res] = system(cmd);
            if s
                error('Error deleting output: %s', res);
            end
            
            if exist(r.fileLFADSOutput, 'file')
                delete(r.fileLFADSOutput);
            end
            
            donefile = fullfile(r.path, 'lfads.done');
            if exist(donefile, 'file')
                delete(donefile);
            end
        end
        
        function makeSequenceFiles(r)
            % Generate the seqence files and save them to disk
            
            if isempty(r.datasets)
                fprintf('No datasets added to Run\n');
                return;
            end
             
            prog = LFADS.Utils.ProgressBar(numel(r.datasets), 'Generating sequence files');
            for iDS = 1:numel(r.datasets)
                prog.update(iDS, 'Generating sequence file for %s', ds.name);
                r.generateSequenceStructForDataset(iDS, true);
            end
            prog.finish();
        end
        
        function deleteSequenceFiles(r)
            % Delete the seqence files saved to disk
            
            if isempty(r.datasets)
                fprintf('No datasets added to Run\n');
            end
             
            for iDS = 1:numel(r.datasets)
                seqFile = fullfile(r.pathSequenceFiles, r.sequenceFileNames{iDS});
                if exist(seqFile, 'file')
                    delete(seqFile); 
                end
            end
        end
        
        function seq = generateSequenceStructForDataset(r, datasetIndex, saveToDisk, mode)
            if nargin < 3
                saveToDisk = false;
            end 
            if nargin < 4
                mode = 'export';
            end
            ds = r.datasets(datasetIndex);

            % call user function on dataset
            seq = r.convertDatasetToSequenceStruct(ds, mode);

            % check the sequence struct returned
            seq = r.checkSequenceStruct(seq);
                
            if saveToDisk
                LFADS.Utils.mkdirRecursive(r.pathSequenceFiles);
                seqFile = fullfile(r.pathSequenceFiles, r.sequenceFileNames{datasetIndex});
                save(seqFile, 'seq');
            end
        end

        function out = loadInputInfo(r)
            fnames = cellfun(@(x) fullfile(r.pathLFADSInput, x), r.lfadsInputInfoFileNames, 'UniformOutput', false);
            for iDS = 1:numel(fnames)
                out(iDS) = load(fnames{iDS}); %#ok<AGROW>
            end
        end

        function seqCell = loadSequenceData(r, reload, mode)
            % seq = loadSequenceData([reload = True])
            % Load the sequence files from disk if they exist, or generates
            % them if not. Caches them in .sequenceData.
            %
            % Args:
            %   reload (bool) : Reload sequence data from disk even if already found in
            %     .sequenceData. Default = false
            %   mode (string) : Either 'export' if data is destined for
            %     LFADS as the raw spike data it will operate on. Or
            %     'alignment' if destined for building alignment matrices,
            %     which may be a subset of the trials (e.g., if only correct
            %     trials should be used for this purpose).
            %
            % Returns:
            %   seqData (cell of struct arrays) : nDatasets cell array of sequence structures loaded from sequence files on disk

            if nargin < 2
                reload = false;
            end
            if nargin < 3
                mode = 'export';
            end
            if strcmp(mode, 'export') && ~reload && ~isempty(r.sequenceData)
                seqCell = r.sequenceData;
                return;
            end
            
            % we will try to load the files from disk if already generated
            % (only on mode == 'export')
            seqFiles = cellfun(@(file) fullfile(r.pathSequenceFiles, file), r.sequenceFileNames, 'UniformOutput', false);
            
            if r.nDatasets > 1
                prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading/generating sequence data for %s', mode);            else
                prog = [];
            end
            seqCell = cell(r.nDatasets, 1);
            for iDS = 1:r.nDatasets
                if ~isempty(prog), prog.update(iDS); end
                if ~exist(seqFiles{iDS}, 'file') || ~strcmp(mode, 'export')
                    % generate on the fly if file not found or mode !=
                    % export
                    seqCell{iDS} = r.generateSequenceStructForDataset(iDS, false, mode);
                else
                    % load from disk if exists and mode == export
                    tmp = load(seqFiles{iDS});
                    seqCell{iDS} = r.checkSequenceStruct(tmp.seq);
                end
            end
            if ~isempty(prog), prog.finish(); end
            
            r.sequenceData = r.modifySequenceDataPostLoading(seqCell);
        end
        
        function seq = checkSequenceStruct(r, seq)
            assert(isfield(seq, 'y'), 'Sequence struct missing y field');
            assert(isfield(seq, 'y_time'), 'Sequence struct missing y_time field');
                
            % convert old params.dtMS --> binWidthMs 
            if ~isfield(seq, 'binWidthMs') 
                if isfield(seq, 'params') % convert params.dtMS to spikeBinWidthMs
                    for iS = 1:numel(seq)
                        seq(iS).binWidthMs = seq(iS).params.dtMS;
                    end 
                else
                    error('Sequence struct missing binWidthMs field');
                end
            end
            assert(numel(unique([seq.binWidthMs])) == 1, 'binWidthMs mismatch');
        end
            
        function seq = modifySequenceDataPostLoading(r, seq)
            % Optionally make any changes or do any post-processing of sequence data upon loading
            
        end
        
        function deleteLFADSInputFiles(r)
            % Delete the seqence files saved to disk
            
            if isempty(r.datasets)
                fprintf('No datasets added to Run\n');
            end
             
            fnames = r.lfadsInputFileNames;
            fnamesInfo = r.lfadsInputInfoFileNames;
            for iDS = 1:numel(r.datasets)
                file = fullfile(r.pathCommonData, fnames{iDS});
                if exist(file, 'file')
                    delete(file); 
                end
                
                file = fullfile(r.pathCommonData, fnamesInfo{iDS});
                if exist(file, 'file')
                    delete(file); 
                end
                
                file = fullfile(r.pathLFADSInput, fnames{iDS});
                if exist(file, 'file')
                    delete(file); 
                end
                
                file = fullfile(r.pathLFADSInput, fnamesInfo{iDS});
                if exist(file, 'file')
                    delete(file); 
                end
            end
        end
        
        function makeLFADSInput(r, regenerate)
            % Generate the LFADS input HD5 files and save them to disk in the pathCommonData folder. 
            % If a file already exists, keep the existing file unless
            % regenerate is true. Then symlink the HD5 files used by this
            % run into pathLFADSInput.
            %
            % Args:
            %   regenerate (bool) : Regenerate HD5 files on disk. If false,
            %     the existing files will be left alone.
            
            if nargin < 2
                regenerate = false;
            end
            
            seqs = {};
            validInds = {};
            trainInds = {};
            
            par = r.params;
            if regenerate
                r.deleteSequenceFiles();
            end
            
            % check which files need to be regenerate
            maskGenerate = false(r.nDatasets, 1);
            fnames = r.lfadsInputFileNames;
            inputInfoNames = r.lfadsInputInfoFileNames;
            if ~regenerate
                for iDS = 1:r.nDatasets
                    fname = fullfile(r.pathCommonData, fnames{iDS});
                    if ~exist(fname, 'file')
                        maskGenerate(iDS) = true;
                    end
                    
                    fname = fullfile(r.pathCommonData, inputInfoNames{iDS});
                    if ~exist(fname, 'file')
                        maskGenerate(iDS) = true;
                    end
                end
            end
                    
            if any(maskGenerate)
                seqData = r.loadSequenceData(regenerate);
                
                % if there are multiple datasets, we need an alignment matrix
                if r.nDatasets > 1 && r.params.useAlignmentMatrix
                    % call out to abstract dataset specific method
                    useAlignMatrices = true;
                    
                    % ask for specific dataset for building the alignment
                    % matrices, which may be a subset of all trials, e.g.
                    % correct trials only. If return is empty, use the full
                    % seqData
                    if r.usesDifferentSequenceDataForAlignment()
                        seqDataForAlignmentMatrices = r.loadSequenceData(regenerate, 'alignment');
                    else
                        seqDataForAlignmentMatrices = seqData;
                    end
                    alignmentMatrices= r.prepareAlignmentMatrices(seqDataForAlignmentMatrices);
                else
                    useAlignMatrices = false;
                end

                % choose validation and training trial indices
                [validIndsCell, trainIndsCell] = deal(cell(r.nDatasets, 1));
                for nd = 1:r.nDatasets
                    allInds = 1:numel(seqData{nd});
                    validIndsCell{nd} = 1 : (r.params.trainToTestRatio+1) : numel(seqData{nd});
                    trainIndsCell{nd} = setdiff(allInds, validIndsCell{nd});
                end

                % support old .params.dtMS field
                if isfield(seqData{1}(1), 'params') && isfield(seqData{1}(1).params, 'dtMS')
                    inputBinSizeMs = seqData{1}(1).params.dtMS;
                elseif isfield(seqData{1}(1), 'binWidthMs')
                    inputBinSizeMs = seqData{1}(1).binWidthMs;
                else
                    error('Sequence data lacks binWidthMs field');
                end

                % arguments for the 'seq_to_lfads' call below
                seqToLFADSArgs = {'binSizeMs', par.spikeBinMs,  ...
                    'inputBinSizeMs', inputBinSizeMs, ...
                    'trainInds', trainIndsCell(maskGenerate), 'testInds', validIndsCell(maskGenerate)};

                if useAlignMatrices
                    seqToLFADSArgs{end+1} = 'alignment_matrix_cxf';
                    seqToLFADSArgs{end+1} = alignmentMatrices(maskGenerate);
                end

                % write the actual lfads input file
                LFADS.Utils.mkdirRecursive(r.pathCommonData);
                LFADS.seq_to_lfads(seqData(maskGenerate), r.pathCommonData, r.lfadsInputFileNames, ...
                    seqToLFADSArgs{:});
                
                % save input info file for each dataset generated
                inputInfoNames = r.lfadsInputInfoFileNames;
                for iDS = 1:r.nDatasets
                    paramInputDataHash = r.params.generateInputDataHash(); %#ok<*NASGU>
                    if maskGenerate(iDS)
                        trainInds = trainIndsCell{iDS};
                        validInds = validIndsCell{iDS};
                        fname = fullfile(r.pathCommonData, inputInfoNames{iDS});
                        save(fname, 'trainInds', 'validInds', 'paramInputDataHash');
                    end
                end
            end
            
            % check which files need to be symlinked from pathCommonData
            % into pathLFADSInput
            LFADS.Utils.mkdirRecursive(r.pathLFADSInput);
            maskLink = true(r.nDatasets, 1);
            fnames = r.lfadsInputFileNames;
            fnamesInputInfo = r.lfadsInputInfoFileNames;
            for iDS = 1:r.nDatasets
                % make link relative (from
                % runCollection/param_HASH/runName/lfadsInput/ to
                % runCollection/data_HASH/)
                origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnames{iDS});
                % origName = fullfile(r.pathCommonData, fnames{iDS});

                linkName = fullfile(r.pathLFADSInput, fnames{iDS});
                if ~exist(linkName, 'file') || regenerate
                    LFADS.Utils.makeSymLink(origName, linkName, false);
                end
                
                % make link relative
                origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnamesInputInfo{iDS});
                % origName = fullfile(r.pathCommonData, fnamesInputInfo{iDS});
                linkName = fullfile(r.pathLFADSInput, fnamesInputInfo{iDS});
                if ~exist(linkName, 'file') || regenerate
                    LFADS.Utils.makeSymLink(origName, linkName, false);
                end
            end
        end
        
        function f = writeShellScriptLFADSTrain(r, varargin)
            % f = writeShellScriptLFADSTrain(cuda_visible_device, display, varargin)
            % Write a shell script used for running the LFADS python code
            %
            % Args:
            %   cuda_visible_devices : int
            %     which GPUs to make visible to CUDA, e.g. 0 or 1
            %   display : int
            %     which display to use for internal LFADS plot generation, e.g. 500 -->
            %     DISPLAY=:500
            %   useTmuxSession : bool
            %     if true, will prefix the command so that it runs within a
            %     new tmux session
            %   appendPosteriorMeanSample : bool (false)
            %     if true, will append the command to run posterior mean
            %     sampling after training is finished
            %
            % Returns:
            %   file : string
            %     Full path to shell script which can be used to begin or resume LFADS training
            
            p = inputParser();
            p.addOptional('cuda_visible_devices', [], @isscalar);
            p.addOptional('display', '', @(x) isnumeric(x) && mod(x,1)==0);
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', true, @islogical);
            p.addParameter('header', '#!/bin/bash', @ischar);
            p.addParameter('appendPosteriorMeanSample', false, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});
            
            f = r.fileShellScriptLFADSTrain;
            fid = fopen(f, 'w');
            trainString = r.buildLFADSTrainingCommand(...
                'cuda_visible_devices', p.Results.cuda_visible_devices, ...
                'display', p.Results.display, ...
                'useTmuxSession', p.Results.useTmuxSession, ...
                'keepSessionAlive', p.Results.keepSessionAlive, ...
                'teeOutput', false); % teeify later
            
            if p.Results.appendPosteriorMeanSample
                % run only if train succeeds
                pmString = r.buildCommandLFADSPosteriorMeanSample(...
                    'cuda_visible_devices', p.Results.cuda_visible_devices, ...
                    'useTmuxSession', p.Results.useTmuxSession, ...
                    'keepSessionAlive', p.Results.keepSessionAlive, ...
                    'teeOutput', false); % teeify later
                if p.Results.teeOutput
                    trainString = sprintf('(%s && %s)', trainString, pmString);
                else
                    trainString = sprintf('%s && %s', trainString, pmString);
                end
            end
            
            if p.Results.teeOutput
                trainString = LFADS.Utils.teeify_string(trainString, r.fileLFADSOutput, false);
            end
            
            fprintf(fid, '%s\n%s\n', p.Results.header, trainString);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end
        
        function runLFADSTrainingCommand(r)
            %   function runLFADSTrainingCommand(r)
            system( sprintf('sh %s', r.fileShellScriptLFADSTrain) );
        end
        
        function outputString = buildLFADSTrainingCommand(r, varargin)
            p = inputParser();
            p.addOptional('cuda_visible_devices', [], @(x) isempty(x) || isscalar(x));
            p.addOptional('display', '', @(x) isempty(x) || (isnumeric(x) && mod(x,1)==0));
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', true, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});
            
            outputString = sprintf(['python $(which run_lfads.py) --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s'], ...
                r.pathLFADSInput, r.pathLFADSOutput);
            
            % use the method from +LFADS/RunParams.m
            optionsString = r.params.generateCommandLineOptionsString(r);
            outputString = sprintf('%s%s', outputString, optionsString);
            
            % set cuda visible devices
            if ~isempty(p.Results.cuda_visible_devices)
                outputString = sprintf('CUDA_VISIBLE_DEVICES=%i %s', ...
                    p.Results.cuda_visible_devices, outputString);
            end
            % set the display variable
            if ~isempty(p.Results.display)
                outputString = sprintf('DISPLAY=:%i %s', ...
                    p.Results.display, outputString);
            end
            
            if p.Results.teeOutput
                outputString = LFADS.Utils.teeify_string(outputString, r.fileLFADSOutput, false);
            end
            
            % if requested, tmux-ify the command
            if p.Results.useTmuxSession
                outputString = LFADS.Utils.tmuxify_string( outputString, r.sessionNameTrain, 'keepSessionAlive', p.Results.keepSessionAlive);
            end
        end
        
        function cmd = buildCommandLFADSPosteriorMeanSample(r, varargin)
            % Generates the command string for LFADS posterior mean sampling
            %
            % Returns
            % --------
            % cmd : string
            %   Shell command for running LFADS posterior mean sampling
            
            p = inputParser();
            p.addOptional('inputParams', @iscell)
            p.addParameter('loadHyperparametersFromFile', false, @islogical);
            p.addParameter('batchSize', 512, @isscalar);
            p.addParameter('cuda_visible_devices', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', false, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});
            batchSize = p.Results.batchSize;
            
            if p.Results.loadHyperparametersFromFile
                % this is the old way of doing it that isn't necessary now
                % that the LFADS run_lfads.py code has stabilized
                
                lfdir = r.pathLFADSOutput;

                % load the params that were used for training            
                params = lfadsi_read_parameters(lfdir); %#ok<*PROPLC>

                % make sure these are up to date
                params.data_dir = r.pathLFADSInput;
                params.lfads_save_dir = r.pathLFADSOutput;

                params.batch_size = batchSize; % this is the number of samples used to calculate the posterior mean
                params.checkpoint_pb_load_name = 'checkpoint_lve';

                % add in allow growth field
                params.allow_gpu_growth = r.params.c_allow_gpu_growth;

                % need to remove "dataset_names" and "dataset_dims" and
                % "temporal_spike_jitter_width"
                params = rmfield(params, {'dataset_names', 'dataset_dims', 'temporal_spike_jitter_width'});
                use_controller = boolean(params.ci_enc_dim);
            
                execstr = 'python';
                if ~isempty(p.Results.inputParams)
                    inputParams = p.Results.inputParams;
                    % take the arguments passed, add new params, or overwrite existing ones
                    f = fields(inputParams);
                    for nn = 1:numel(f)
                        params.(f{nn}) = inputParams.(f{nn});
                    end
                end
                
                f = fields(params);
                optstr = '';
                for nf = 1:numel(f)
                    fval = params.(f{nf});
                    %convert any numbers to strings
                    if islogical(fval)
                        if fval
                            fval = 'True';
                        else
                            fval = 'False';
                        end
                    elseif isnumeric(fval)
                        fval = num2str(fval);
                    end
                    optstr = strcat(optstr, sprintf(' --%s=%s',f{nf}, fval));
                end

                optstr = strcat(optstr, sprintf(' --kind=posterior_sample'));

                % put the command together
                cmd = sprintf('%s $(which run_lfads.py) %s', execstr, optstr);
            else
                % use the RunParams to generate the params
                paramsString = r.params.generateCommandLineOptionsString(r, 'omitFields', {'c_temporal_spike_jitter_width', 'batch_size'});
                
                cmd = sprintf(['python $(which run_lfads.py) --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s --kind=posterior_sample --batch_size=%d --checkpoint_pb_load_name=checkpoint_lve %s'], ...
                r.pathLFADSInput, r.pathLFADSOutput, batchSize, paramsString);
            end        
            
            % set cuda visible devices
            if ~isempty(p.Results.cuda_visible_devices)
                cmd = sprintf('CUDA_VISIBLE_DEVICES=%i %s', ...
                    p.Results.cuda_visible_devices, cmd);
            end
            
            if p.Results.teeOutput
                cmd = LFADS.Utils.teeify_string(cmd, r.fileLFADSOutput, true);
            end
            
            % if requested, tmux-ify the command
            if p.Results.useTmuxSession
                cmd = LFADS.Utils.tmuxify_string(cmd, r.sessionNamePosteriorMean, 'keepSessionAlive', p.Results.keepSessionAlive );
                fprintf('Tmux Session is %s\n  tmux a -t %s\n\n', r.sessionNamePosteriorMean, r.sessionNamePosteriorMean);
            end
        end
        
        function cmd = buildCommandLFADSWriteModelParams(r)
            % Generates the command string for LFADS write model params
            %
            % Returns
            % --------
            % cmd : string
            %   Shell command for running LFADS write model params

            % use the RunParams to generate the params
            parstr = r.params.generateCommandLineOptionsString(r, 'omitFields', {'c_temporal_spike_jitter_width'});

            cmd = sprintf(['python $(which run_lfads.py) --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s --kind=write_model_params --checkpoint_pb_load_name=checkpoint_lve %s'], ...
                r.pathLFADSInput, r.pathLFADSOutput, parstr);
        end
        
        function runLFADSPosteriorMeanCommand(r)
            %   function runLFADSPosteriorMeanCommand(r)
            system( sprintf('sh %s', r.fileShellScriptLFADSPosteriorMeanSample) );
        end
        
        function f = writeShellScriptLFADSPosteriorMeanSample(r, varargin)
            % Write a shell script used for running the LFADS posterior mean sampling. 
            %
            % Args:
            %   cuda_visible_devices : int
            %     which GPUs to make visible to CUDA, e.g. 0 or 1
            %   useTmuxSession : bool
            %     if true, will prefix the command so that it runs within a
            %     new tmux session
            %
            % Returns
            % --------
            % file : string
            %   Full path to shell script which can be used to perform LFADS posterior mean sampling on the lowest
            %   validation error checkpoint
           
            p = inputParser();
            p.addOptional('cuda_visible_devices', [], @isscalar);
            p.addParameter('header', '#!/bin/bash\n', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            f = r.fileShellScriptLFADSPosteriorMeanSample;
            fid = fopen(f, 'w');
            
            outputString = r.buildCommandLFADSPosteriorMeanSample('cuda_visible_devices', p.Results.cuda_visible_devices, p.Unmatched);
            
            fprintf(fid, [p.Results.header, outputString]);
            fclose(fid);
            LFADS.Utils.chmod('ug+rx', f);
        end
        
        function runLFADSWriteModelParamsCommand(r)
            system(sprintf('sh %s', r.fileShellScriptLFADSWriteModelParams));
        end
        
        function f = writeShellScriptLFADSWriteModelParams(r, varargin)
            % Write a shell script used for running the LFADS write model params.
            %
            % Args:
            %   cuda_visible_devices : int
            %     which GPUs to make visible to CUDA, e.g. 0 or 1
            %   useTmuxSession : bool
            %     if true, will prefix the command so that it runs within a
            %     new tmux session
            %
            % Returns
            % --------
            % file : string
            %   Full path to shell script which can be used to perform LFADS posterior mean sampling on the lowest
            %   validation error checkpoint
           
            p = inputParser();
            p.addOptional('cuda_visible_devices', [], @isscalar);
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', false, @islogical);
            p.addParameter('header', '#!/bin/bash\n', @ischar);
            p.parse(varargin{:});
            
            f = r.fileShellScriptLFADSWriteModelParams;
            fid = fopen(f, 'w');
            
            outputString = r.buildCommandLFADSWriteModelParams();

            % set cuda visible devices
            if ~isempty(p.Results.cuda_visible_devices)
                outputString = sprintf('CUDA_VISIBLE_DEVICES=%i %s', ...
                    p.Results.cuda_visible_devices, outputString);
            end
            
            % if requested, tmux-ify the command
            if p.Results.useTmuxSession
                outputString = LFADS.Utils.tmuxify_string( outputString, r.sessionNameWriteModelParams, 'keepSessionAlive', p.Results.keepSessionAlive );
                fprintf('Tmux Session is %s\n  tmux a -t %s\n\n', r.sessionNameWriteModelParams, r.sessionNameWriteModelParams);
            end
            
            fprintf(fid, [p.Results.header outputString]);
            fclose(fid);
            LFADS.Utils.chmod('ug+rx', f);
        end
        
        function pms = loadPosteriorMeans(r, reload)
            % pmData = loadPosteriorMeans(r, reload)
            % After the posterior mean shell script has been run, this will load the posterior mean samples from disk
            % and convert them into :ref:`LFADS_PosteriorMeans` instances. These will also be cached in r.posteriorMeans
            %
            % Parameters
            % ------------
            % reload : bool
            %   if false, data stored in r.posteriorMeans will be returned
            %   if all datasets are loaded. if true, new data will be
            %   loaded from disk always.
            %
            % Returns
            % --------
            % pmData : string
            %   nDatasets cell of :ref:`LFADS_PosteriorMeans` data loaded from disk
            
            if nargin < 2
                reload = false;
            end
            if ~isempty(r.posteriorMeans) && all([r.posteriorMeans.isValid]) && ~reload
                pms = LFADS.Utils.makecol(r.posteriorMeans);
                return;
            end
            
            seq = r.loadSequenceData();
            
            info = r.loadInputInfo();
            % check hashes actually match
            thisHash = r.params.generateInputDataHash();
            for iDS = 1:r.nDatasets
                if ~isequal(info(iDS).paramInputDataHash, thisHash)
                    error('Input data param hash saved for run %d in %s does not match', iDS, r.lfadsInputInfoFileNames{1});
                end
            end
            
            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
            prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading posterior means for each dataset');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                if ~exist(fullfile(r.pathLFADSOutput, trainList{iDS}), 'file')
                    warning('LFADS Posterior Mean train file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, trainList{iDS}));
                    pms = [];
                    break;
                end
                if ~exist(fullfile(r.pathLFADSOutput, validList{iDS}), 'file')
                    warning('LFADS Posterior Mean valid file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, validList{iDS}));
                    pms = [];
                    break;
                end
                
                pmData = LFADS.Utils.loadPosteriorMeans(fullfile(r.pathLFADSOutput, validList{iDS}), ....
                    fullfile(r.pathLFADSOutput, trainList{iDS}), ...
                    info(iDS).validInds, info(iDS).trainInds);
                
                if isfield(seq{iDS}, 'binWidthMs')
                    dt_y = seq{iDS}(1).binWidthMs;
                else
                    dt_y = seq{iDS}(1).params.dtMS;
                end
                dt_pm = r.params.spikeBinMs;
                rebin = dt_pm / dt_y;
                time = seq{iDS}(1).y_time(1:rebin:end);
                time = time(1:size(pmData.rates, 2));
                
                pms(iDS) = LFADS.PosteriorMeans(pmData, r.params, time); %#ok<AGROW>
            end
            prog.finish();
            
            pms = LFADS.Utils.makecol(pms);
            r.posteriorMeans = pms;
        end
        
        function seqs = addPosteriorMeansToSeq(r)
            % function seqs = addPosteriorMeansToSeq(r)
            % returns a sequence that has posterior mean
            % values integrated
            
            if isempty(r.posteriorMeans) || ~all([r.posteriorMeans.isValid])
                r.loadPosteriorMeans;
            end
            
            if isempty(r.sequenceData) || numel(r.sequenceData) == 0
                r.sequenceData = r.loadSequenceData();
            end
            seqs = r.sequenceData;
            
            % iterate over datasets
            for iDS = 1:numel(seqs)
                pm = r.posteriorMeans(iDS);
                
                for ntr = 1:numel(seqs{iDS})
                    seqs{iDS}(ntr).rates = squeeze(pm.rates(:,:,ntr));
                    seqs{iDS}(ntr).factors = squeeze(pm.factors(:,:,ntr));
                    seqs{iDS}(ntr).generator_states = squeeze(pm.generator_states(:,:,ntr));
                    if ~isempty(pm.controller_outputs)
                        seqs{iDS}(ntr).controller_outputs = ...
                            squeeze(pm.controller_outputs(:,:,ntr));
                    else
                        seqs{iDS}(ntr).controller_outputs = [];
                    end
                end
            end
            
            r.sequenceData = seqs;
        end
        
        function readouts = loadReadoutMatricesByDataset(r)
            fname = r.fileModelParams;
            assert(exist(fname, 'file') > 1, 'model_params file not found. Ensure that runCommandLFADSWriteModelParams has been run');
            
            biasStr = '/LFADS_glm_fac_2_logrates_%s.h5_b:0';
            weightsStr = '/LFADS_glm_fac_2_logrates_%s.h5_W:0';
            
            for iD = 1:r.nDatasets
                readouts(iD).rates_b = h5read(fname, sprintf(biasStr, r.datasets(iD).name)); %#ok<AGROW>
                readouts (iD).rates_W = h5read(fname, sprintf(weightsStr, r.datasets(iD).name)); %#ok<AGROW>
            end
            
            readouts = readouts';
            
        end
            
    end
end
