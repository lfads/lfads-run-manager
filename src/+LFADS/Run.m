classdef Run < handle & matlab.mixin.CustomDisplay
    % Represents a single LFADS experiment on a specific set of datasets. Runs are grouped into :ref:`LFADS_RunCollection`
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.

    methods
        % These methods will need to be implemented in a subclass that provides the custom behavior for your application

        function out = generateCountsForDataset(r, dataset, mode, varargin) %#ok<STOUT,INUSD>
            % Generate binned spike count tensor for a single dataset.
            %% Generate binned spike count tensor for a single dataset.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            %
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices. You can
            %   include a different subset of the data (or different time
            %   windows) for the alignment process separately from the actual
            %   data exported to LFADS, or return the same for both. Alignment
            %   is only relevant for multi-dataset models. If you wish to use
            %   separate data for alignment, override the method usesDifferentDataForAlignment
            %   to return true as well.
            %
            % Returns
            % ----------
            % out: a scalar struct with the following fields:
            % .counts : nTrials x nChannels x nTime tensor
            %   spike counts in time bins in trials x channels x time. These
            %   should be total counts, not normalized rates, as they will be
            %   added during rebinning.
            %
            % .timeVecMs: nTime x 1 vector
            %   of timepoints in milliseconds associated with each time bin. You can start this
            %   wherever you like, but timeVecMs(2) - timeVecMs(1) will be
            %   treated as the spike bin width used when the data are later
            %   rebinned to match run.params.spikeBinMs
            %
            % .conditionId: nTrials x 1 vector
            %   of unique conditionIds. Can be cell array of strings or
            %   vector of unique integers.
            %
            % .truth: nTrials x nChannels x nTime tensor
            %   for synthetic datasets, provides the ground-truth counts
            %   for each trial. Same size as .counts. 
            %
            % .externalInputs: nTrials x nExternalInputs x nTime tensor
            %   specifies the observed, external inputs which will be
            %   passed either to the generator directly or to the encoder.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            %
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices. You can
            %   include a different subset of the data (or different time
            %   windows) for the alignment process separately from the actual
            %   data exported to LFADS, or return the same for both. Alignment
            %   is only relevant for multi-dataset models. If you wish to use
            %   separate data for alignment, override the method usesDifferentDataForAlignment
            %   to return true as well.
            %
            % Returns
            % ----------
            % counts : nTrials x nChannels x nTime tensor
            %   spike counts in time bins in trials x channels x time. These
            %   should be total counts, not normalized rates, as they will be
            %   added during rebinning.
            %
            % timeVecMs: nTime x 1 vector
            %   of timepoints in milliseconds associated with each time bin. You can start this
            %   wherever you like, but timeVecMs(2) - timeVecMs(1) will be
            %   treated as the spike bin width used when the data are later
            %   rebinned to match run.params.spikeBinMs
            %
            % conditionId: nTrials x 1 vector
            %   of unique conditionIds. Can be cell array of strings or
            %   vector of unique integers.

            error('Implement generateCountsForDataset in your subclass');
        end

    end

    properties
        name char = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk

        comment char = '' % Textual comment for convenience

        version uint32 = 20171107; % Internal versioning allowing for graceful evolution of path settings, set by RunCollection automatically
    end

    properties
        runCollection % :ref:`LFADS_RunCollection` instance to which this run belongs
        params % :ref:`LFADS_RunParams` instance shared by all runs in the collection, contains parameter settings
        paramIndexInRunCollection % which RunParams index in the RunCollection this params corresponds to
        datasets % Array of :ref:`LFADS_Dataset` instances which this particular Run will utilize

        sequenceData cell % nDatasets cell array of sequence struct data
        posteriorMeans % nDatasets array of :ref:`LFADS_PosteriorMeans` when loaded

        inputInfo % parameters used to train model

        multisessionAlignmentTool % LFADS.MultisessionAlignmentTool instance used to generate alignment matrices and visualize results
        modelTrainedParams % LFADS.ModelTrainedParams instance holding the learned LFADS network parameters
    end

    properties(Dependent)
        nDatasets % Number of datasets used by this run
        datasetNames % nDatasets x 1 cell array of names
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
        
        fileFitLog % location on disk where fit log will be written

        sessionNameTrain % name of tmux session that will be created if useSession = true is passed to writeShellScriptLFADSTrain
        sessionNamePosteriorMean % name of tmux session that will be created if useTmuxSession = true is passed to writeShellScriptLFADSPosteriorMean
        sessionNameWriteModelParams % name of tmux session that will be created if useTmuxSession = true is passed to writeShellScriptLFADSWriteModelParams
    end

    methods
        function r = Run(varargin)

        end

        function tf = usesDifferentDataForAlignment(r)  %#ok<MANU>
            % tf = usesDifferentDataForAlignment()
            %
            % Returns true if you would like the Run to call your
            % convertDatasetToSequenceStruct with a mode == 'alignment'
            % argument when constructing the alignment matrices. This
            % allows you to specify a different set of trials used for
            % constructing alignment matrices, e.g. only correct trials.

            tf = false;
        end

        function seq = convertDatasetToSequenceStruct(r, dataset, mode, varargin)
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
            %       `y_true` has the same size as `y` and provides
            %       ground-truth counts.  `externalInputs` provides
            %       observed inputs.
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

            out = r.generateCountsForDataset(dataset, mode, varargin{:});
            
            assert(isstruct(out) && isscalar(out) && isfield(out, 'counts'));
            counts = out.counts;
            if isfield(out, 'timeVecMs')
                timeVecMs = out.timeVecMs;
            else
                timeVecMs = 1:size(counts, 3);
            end
                
            if isfield(out, 'conditionId') 
                conditionId = out.conditionId;
            else
                conditionId = [];
            end
            
            if isfield(out, 'truth')
                truth = out.truth;
            else
                truth = [];
            end
            
            if isfield(out, 'externalInputs')
                externalInputs = out.externalInputs;
            else
                externalInputs = [];
            end
            
            assert(isnumeric(counts) && ndims(counts) == 3, 'counts must be 3 dim tensor');
            assert(isnumeric(timeVecMs) && isvector(timeVecMs), 'timeVecMs must be numeric vector');
            assert(size(counts, 3) == numel(timeVecMs), 'size(counts, 3) must match numel(timeVecMs)');
            assert(isempty(conditionId) || numel(conditionId) == size(counts, 1), 'numel(conditionId) must be size(counts, 1)');
            assert(isempty(truth) || isequal(size(counts), size(truth)), 'size(truth) must be same as size(counts)');
            assert(isempty(externalInputs) || (size(counts, 1) == size(externalInputs, 1) && size(counts, 2) == size(counts, 2)), 'size(externalInputs) must match size(counts) along dims 1, 2');

            nTrials = size(counts, 1);

            binWidthMs = timeVecMs(2) - timeVecMs(1);

            for iTrial = nTrials:-1:1
                seq(iTrial).y = squeeze(counts(iTrial, :, :)); % nChannels x nTime
                seq(iTrial).y_time = timeVecMs;
                seq(iTrial).binWidthMs = binWidthMs;
                if ~isempty(conditionId)
                    if iscell(conditionId)
                        seq(iTrial).conditionId = conditionId{iTrial};
                    else
                        seq(iTrial).conditionId = conditionId(iTrial);
                    end
                end
                if ~isempty(truth)
                     seq(iTrial).y_true = squeeze(truth(iTrial, :, :)); % nChannels x nTime
                end
                if ~isempty(externalInputs)
                     seq(iTrial).externalInputs = squeeze(externalInputs(iTrial, :, :));  % nExternalInputs x nTime
                end
            end
        end

        function [alignmentMatrices, alignmentBiases] = doMultisessionAlignment(r, regenerate)
            if nargin < 2
                regenerate = false;
            end
            assert((r.nDatasets > 1 && r.params.useAlignmentMatrix) || (r.nDatasets == 1 && r.params.useSingleDatasetAlignmentMatrix), ...
                'Alignment matrices can only be built for multi-dataset runs where useAlignmentMatrix is true or single-dataset runs where useSingleDatasetAlignmentMatrix is true');

            % ask for specific dataset for building the alignment
            % matrices, which may be a subset of all trials, e.g.
            % correct trials only.
            seqDataForAlignmentMatrices = r.loadSequenceData(regenerate, 'alignment');
            [alignmentMatrices, alignmentBiases] = r.prepareAlignmentMatrices(seqDataForAlignmentMatrices);
        end

        function [alignmentMatrices, alignmentBiases] = prepareAlignmentMatrices(r, seqData)
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

            r.multisessionAlignmentTool = LFADS.MultisessionAlignmentTool(r, seqData, {r.datasets.name}');
            [alignmentMatrices, alignmentBiases] = r.multisessionAlignmentTool.computeAlignmentMatricesUsingTrialAveragedPCR();
        end

        function visualizeAlignmentMatrixReconstructions(r, nFactorsOrFactorIdx, nConditionsOrConditionIdx)
            if isempty(r.multisessionAlignmentTool)
                r.doMultisessionAlignment();
            end

            r.multisessionAlignmentTool.plotAlignmentReconstruction(nFactorsOrFactorIdx, nConditionsOrConditionIdx);
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
            elseif r.version < 20171107
                p = fullfile(r.runCollection.path, r.params.generateInputDataHashName());
            else
                % append run name under data_HASH folder so that different
                % run specs get different h5 files (with potentially
                % different alignment matrices)
                p = fullfile(r.runCollection.path, r.params.generateInputDataHashName(), r.name);
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
        
        function f = get.fileFitLog(r)
            f = fullfile(r.pathLFADSOutput, 'fitlog.csv');
        end


        function f = get.fileLFADSOutput(r)
            f = fullfile(r.path, 'lfads.out');
        end

        function n = get.nDatasets(r)
            n = numel(r.datasets);
        end

        function names = get.datasetNames(r)
            names = {r.datasets.name}';
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
                trainList = arrayfun(@(ds) sprintf('model_runs_%s.h5_train_posterior_sample_and_average', ds.name), ...
                    r.datasets, 'UniformOutput', false);
                validList = arrayfun(@(ds) sprintf('model_runs_%s.h5_valid_posterior_sample_and_average', ds.name), ...
                    r.datasets, 'UniformOutput', false);
            end
        end

        function sess = get.sessionNameTrain(r)
            sess = sprintf('train_%s__%s', r.name, r.paramsString);
            sess = strrep(sess, '.', '_');
        end

        function sess = get.sessionNamePosteriorMean(r)
            sess = sprintf('pm_%s__%s', r.name, r.paramsString);
            sess = strrep(sess, '.', '_');
        end

        function sess = get.sessionNameWriteModelParams(r)
            sess = sprintf('writeParams_%s_%s', r.name, r.paramsString);
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
            % loads fields saved into inputInfo.mat which are essentially
            % cached things needed for post-processing, e.g. training vs.
            % validation trial inds, posterior mean time vectors
            %
            % Returns:
            %   out (nDatasets x 1 struct array)

            fnames = cellfun(@(x) fullfile(r.pathLFADSInput, x), r.lfadsInputInfoFileNames, 'UniformOutput', false);
            for iDS = 1:numel(fnames)
                out(iDS) = load(fnames{iDS}); %#ok<AGROW>
            end
            out = out';
            r.inputInfo = out;
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
            
            if strcmp(mode, 'alignment') && ~r.usesDifferentDataForAlignment()
                % subsequent steps should go identically for alignment if
                % usesDifferentDataForAlignment returns false
                mode = 'export';
            end
            
            if strcmp(mode, 'export') && ~reload && ~isempty(r.sequenceData)
                % return cached sequence data if already loaded
                seqCell = r.sequenceData;
                return;
            end

            % we will try to load the files from disk if already generated
            % (only on mode == 'export')
            seqFiles = cellfun(@(file) fullfile(r.pathSequenceFiles, file), r.sequenceFileNames, 'UniformOutput', false);

            if r.nDatasets > 1
                prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading/generating sequence data for %s', mode); 
            else
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

            if strcmp(mode, 'export')
                r.sequenceData = r.modifySequenceDataPostLoading(seqCell);
            end
        end

        function seq = checkSequenceStruct(r, seq) %#ok<INUSL>
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

        function seq = modifySequenceDataPostLoading(r, seq) %#ok<INUSL>
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
                % delete the original files
                file = fullfile(r.pathCommonData, fnames{iDS});
                if exist(file, 'file')
                    delete(file);
                end

                file = fullfile(r.pathCommonData, fnamesInfo{iDS});
                if exist(file, 'file')
                    delete(file);
                end

                % and delete the symlinks that point to them
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

        function assertParamsOkayForSequenceData(r, seqData)
            % Checks parameter settings against the sequence data.
            % 1. are there enough trials given the batch size and
            % trainTestRatio

            nTrials = cellfun(@numel, seqData);
            nRequired = r.params.c_batch_size * (1+r.params.trainToTestRatio);

            idxTooFew = find(nTrials <= nRequired); % lfads code uses a > rather than >=, so we must do the same

            assert(isempty(idxTooFew), 'Issue with Run %s: %d trials are required for c_batch_size=%d and trainToTestRatio=%d. Datasets %s have too few trials', ...
                r.name, nRequired, r.params.c_batch_size, r.params.trainToTestRatio, vec2str(idxTooFew));
            
            function str = vec2str(vec)
                % returns a string representation of a vector

                str = ['[' num2str(LFADS.Utils.makerow(vec)) ']'];

                if size(vec, 1) > size(vec, 2)
                    % include a transpose tick if its a column vector
                    str = [str ''''];
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
            else
                maskGenerate = true(r.nDatasets, 1);
            end

            if any(maskGenerate)
                regenerate = false;
                seqData = r.loadSequenceData(regenerate); % this will set r.sequenceData

                r.assertParamsOkayForSequenceData(seqData);

                % no need to regenerate if alignment and export use the
                % same data, since we would have just regenerated them via
                % loadSequenceData above
                regenerateAlignmentData = regenerate && r.usesDifferentDataForAlignment();
                
                if r.nDatasets > 1 && r.params.useAlignmentMatrix
                    % generate alignment matrices for stitching run 
                    useAlignMatrices = true;
                    [alignmentMatrices, alignmentBiases] = r.doMultisessionAlignment(regenerateAlignmentData);
                    
                elseif r.version >= 20171107 && r.nDatasets == 1 && r.params.useSingleDatasetAlignmentMatrix
                    % generate alignment matrix for single run (just PCA down to c_factors_dim)
                    useAlignMatrices = true;
                    [alignmentMatrices, alignmentBiases] = r.doMultisessionAlignment(regenerateAlignmentData);
                    
                else
                    % no alignment matrices
                    useAlignMatrices = false;
                end

                % choose validation and training trial indices
                [validIndsCell, trainIndsCell] = deal(cell(r.nDatasets, 1));
                for iDS = 1:r.nDatasets
                    allInds = 1:numel(seqData{iDS});
                    validIndsCell{iDS} = 1 : (r.params.trainToTestRatio+1) : numel(seqData{iDS});
                    trainIndsCell{iDS} = setdiff(allInds, validIndsCell{iDS});
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
                    seqToLFADSArgs{end+1} = 'alignment_bias_c';
                    seqToLFADSArgs{end+1} = alignmentBiases(maskGenerate);
                end

                % write the actual lfads input file
                LFADS.Utils.mkdirRecursive(r.pathCommonData);
                
                LFADS.Interface.seq_to_lfads(seqData(maskGenerate), r.pathCommonData, r.lfadsInputFileNames, ...
                    seqToLFADSArgs{:});

                % save input info file for each dataset generated
                inputInfoNames = r.lfadsInputInfoFileNames;
                for iDS = 1:r.nDatasets
                    paramInputDataHash = r.params.generateInputDataHash(); %#ok<*NASGU>
                    if maskGenerate(iDS)
                        trainInds = trainIndsCell{iDS};
                        validInds = validIndsCell{iDS};

                        % save time vectors used in the sequence files to
                        % facilitate fast loading of posterior means sampling
                        seq_timeVector = seqData{iDS}.y_time;
                        seq_binSizeMs = inputBinSizeMs;

                        % save the rebinned spike counts and condition ids
                        % too
                        counts = cat(3, seqData{iDS}.y); % nNeurons x nTime x nChannels
                        if isnumeric(seqData{iDS}(1).conditionId)
                            conditionId = cat(1, seqData{iDS}.conditionId);
                        else
                            conditionId = LFADS.Utils.makecol({seqData{iDS}.conditionId});
                        end
                        
                        % optionally include ground truth
                        extra = {};
                        if isfield(seqData{iDS}, 'y_true')
                            truth = seqData{iDS}.y_true;
                            extra = union(extra, 'truth');
                        end
                        if isfield(seqData{iDS}, 'externalInputs')
                            externalInputs = seqData{iDS}.externalInputs;
                            extra = union(extra, 'externalInputs');
                        end

                        fname = fullfile(r.pathCommonData, inputInfoNames{iDS});
                        save(fname, 'trainInds', 'validInds', 'paramInputDataHash', 'seq_timeVector', 'seq_binSizeMs', 'conditionId', 'counts', extra{:});
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
                % make link relative (from link location in
                % runCollection/param_HASH/runName/lfadsInput/)
                if r.version < 20171107
                    % to % runCollection/data_HASH/file.h5
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnames{iDS});
                else
                    % % runCollection/data_HASH/run_name/file.h5
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnames{iDS});
                end

                linkName = fullfile(r.pathLFADSInput, fnames{iDS});
                if ~exist(linkName, 'file') || regenerate
                    LFADS.Utils.makeSymLink(origName, linkName, false);
                end

                % make link relative
                if r.version < 20171107
                    % runCollection/data_HASH/inputInfo.mat
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), fnamesInputInfo{iDS});
                else
                    % runCollection/data_HASH/run_name/inputInfo.mat
                    origName = fullfile('..', '..', '..', r.params.generateInputDataHashName(), r.name, fnamesInputInfo{iDS});
                end
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
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addParameter('appendPosteriorMeanSample', false, @islogical);
            p.addParameter('appendWriteModelParams', false, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.addParameter('prependPathToLFADS', false, @islogical); % prepend an export path to run_lfads.py
            p.addParameter('virtualenv', '', @ischar); % prepend source activate environment name
            p.parse(varargin{:});

            runLFADSPath = LFADS.Utils.find_run_lfads_py();

            f = r.fileShellScriptLFADSTrain;
            fid = fopen(f, 'w');
            trainString = r.buildLFADSTrainingCommand(...
                'path_run_lfads_py', p.Results.path_run_lfads_py, ...
                'cuda_visible_devices', p.Results.cuda_visible_devices, ...
                'display', p.Results.display, ...
                'useTmuxSession', false, ...
                'teeOutput', false); % teeify later

            if p.Results.appendPosteriorMeanSample
                % run only if train succeeds
                pmString = r.buildCommandLFADSPosteriorMeanSample(...
                    'path_run_lfads_py', p.Results.path_run_lfads_py, ...
                    'cuda_visible_devices', p.Results.cuda_visible_devices, ...
                    'useTmuxSession', false, ...
                    'teeOutput', false); % teeify later

                if p.Results.appendWriteModelParams
                    % run only if train succeeds
                    writeParamsString = r.buildCommandLFADSWriteModelParams(...
                        'path_run_lfads_py', p.Results.path_run_lfads_py, ...
                        'cuda_visible_devices', p.Results.cuda_visible_devices, ...
                        'useTmuxSession', false, ...
                        'teeOutput', false); % teeify later
                    if p.Results.teeOutput
                        trainString = sprintf('(%s && %s && %s)', trainString, pmString, writeParamsString);
                    else
                        trainString = sprintf('%s && %s && %s', trainString, pmString, writeParamsString);
                    end
                else
                    if p.Results.teeOutput
                        trainString = sprintf('(%s && %s)', trainString, pmString);
                    else
                        trainString = sprintf('%s && %s', trainString, pmString);
                    end
                end
            end

            if p.Results.teeOutput
                trainString = LFADS.Utils.teeify_string(trainString, r.fileLFADSOutput, false);
            end

            % we do all of the tmux at once on the combined commands
             if p.Results.useTmuxSession
                trainString = LFADS.Utils.tmuxify_string(trainString, r.sessionNameTrain, 'keepSessionAlive', p.Results.keepSessionAlive);
            end

            fprintf(fid, '%s\n', p.Results.header);

            if ~isempty(p.Results.virtualenv)
                fprintf(fid, 'source activate %s\n', p.Results.virtualenv);
            end

            if p.Results.prependPathToLFADS
                folder = LFADS.Utils.find_run_lfads_py(false);
                if ~isempty(folder)
                    fprintf(fid, 'export PATH="%s:$PATH"\n', folder);
                end
            end
            
            checkLFADSFoundString = LFADS.Run.generateCheckRunLFADSPyFoundString(p.Results.path_run_lfads_py);
            fprintf(fid, '\n%s\n', checkLFADSFoundString);

            fprintf(fid, '%s\n', trainString);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end

        function runLFADSTrainingCommand(r)
            %   function runLFADSTrainingCommand(r)
            system( sprintf('sh %s', r.fileShellScriptLFADSTrain) );
        end

        function outputString = buildLFADSTrainingCommand(r, varargin)
            p = inputParser();
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addOptional('cuda_visible_devices', [], @(x) isempty(x) || isscalar(x));
            p.addOptional('display', '', @(x) isempty(x) || (isnumeric(x) && mod(x,1)==0));
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', true, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});

            outputString = sprintf(['python %s --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s'], ...
                p.Results.path_run_lfads_py, ...
                LFADS.Utils.GetFullPath(r.pathLFADSInput), LFADS.Utils.GetFullPath(r.pathLFADSOutput));

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
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addParameter('loadHyperparametersFromFile', false, @islogical);
            p.addParameter('num_samples_posterior', r.params.num_samples_posterior, @isscalar); % can be used to manually overwrite
            p.addParameter('cuda_visible_devices', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', false, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});

            if p.Results.loadHyperparametersFromFile
                % this is the old way of doing it that isn't necessary now
                % that the LFADS run_lfads.py code has stabilized

                lfdir = r.pathLFADSOutput;

                % load the params that were used for training
                params = LFADS.Interface.read_parameters(lfdir); %#ok<*PROPLC>

                % make sure these are up to date
                params.data_dir = LFADS.Utils.GetFullPath(r.pathLFADSInput);
                params.lfads_save_dir = LFADS.Utils.GetFullPath(r.pathLFADSOutput);

                params.checkpoint_pb_load_name = 'checkpoint_lve';
                params.batch_size = p.Results.num_samples_posterior;

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

                optstr = strcat(optstr, sprintf(' --kind=posterior_sample_and_average'));

                % put the command together
                cmd = sprintf('%s %s %s', execstr, p.Results.path_run_lfads_py, optstr);
            else
                % use the RunParams to generate the params
                paramsString = r.params.generateCommandLineOptionsString(r, 'omitFields', {'c_temporal_spike_jitter_width', 'c_batch_size'});
                
                cmd = sprintf(['python %s --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s --kind=posterior_sample_and_average --batch_size=%d --checkpoint_pb_load_name=checkpoint_lve %s'], ...
                p.Results.path_run_lfads_py, ...
                LFADS.Utils.GetFullPath(r.pathLFADSInput), LFADS.Utils.GetFullPath(r.pathLFADSOutput), p.Results.num_samples_posterior, paramsString);
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

        function cmd = buildCommandLFADSWriteModelParams(r, varargin)
            p = inputParser();
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addParameter('cuda_visible_devices', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', false, @islogical);
            p.addParameter('teeOutput', false, @islogical);
            p.parse(varargin{:});

            % Generates the command string for LFADS write model params
            %
            % Returns
            % --------
            % cmd : string
            %   Shell command for running LFADS write model params


            % use the RunParams to generate the params
            paramsString = r.params.generateCommandLineOptionsString(r, 'omitFields', {'c_temporal_spike_jitter_width'});

            cmd = sprintf(['python %s --data_dir=%s --data_filename_stem=lfads ' ...
            '--lfads_save_dir=%s --kind=write_model_params --checkpoint_pb_load_name=checkpoint_lve %s'], ...
                p.Results.path_run_lfads_py, ...
                LFADS.Utils.GetFullPath(r.pathLFADSInput), LFADS.Utils.GetFullPath(r.pathLFADSOutput), paramsString);

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
                fprintf('Tmux Session is %s\n  tmux a -t %s\n\n', r.sessionNameWriteModelParams, r.sessionNameWriteModelParams);
            end
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
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addOptional('cuda_visible_devices', [], @isscalar);
            p.addParameter('header', '#!/bin/bash\n', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            f = r.fileShellScriptLFADSPosteriorMeanSample;
            fid = fopen(f, 'w');

            outputString = r.buildCommandLFADSPosteriorMeanSample(...
                'path_run_lfads_py', p.Results.path_run_lfads_py, ...
                'cuda_visible_devices', p.Results.cuda_visible_devices, p.Unmatched);

            fprintf(fid, '%s\n', p.Results.header);
            checkLFADSFoundString = LFADS.Run.generateCheckRunLFADSPyFoundString(p.Results.path_run_lfads_py);
            fprintf(fid, '\n%s\n', checkLFADSFoundString);
            
            fprintf(fid, '%s\n', outputString);
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
            p.addParameter('path_run_lfads_py', '$(which run_lfads.py)', @ischar);
            p.addParameter('header', '#!/bin/bash\n', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            f = r.fileShellScriptLFADSWriteModelParams;
            fid = fopen(f, 'w');

            outputString = r.buildCommandLFADSWriteModelParams(...
                'path_run_lfads_py', p.Results.path_run_lfads_py, ...
                p.Unmatched);

            fprintf(fid, '%s\n', p.Results.header);
            checkLFADSFoundString = LFADS.Run.generateCheckRunLFADSPyFoundString(p.Results.path_run_lfads_py);
            fprintf(fid, '\n%s\n', checkLFADSFoundString);
            
            fprintf(fid, '%s\n', outputString);
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

            info = r.loadInputInfo();
            % check hashes actually match
            thisHash = r.params.generateInputDataHash();
            for iDS = 1:r.nDatasets
                if ~isequal(info(iDS).paramInputDataHash, thisHash)
                    warning('Input data param hash saved for run %d in %s does not match', iDS, r.lfadsInputInfoFileNames{1});
                end
            end

            % determine whether we need to load the sequence data. if the
            % cached field pm_timeVector is prsent in info, we don't need
            % them
            if ~isfield(info, 'seq_timeVector') || ~isfield(info, 'seq_binSizeMs')
                seq = r.loadSequenceData();
                for iDS = 1:r.nDatasets
                    if isfield(seq{iDS}, 'binWidthMs')
                        info(iDS).seq_binSizeMs = seq{iDS}(1).binWidthMs;
                    else
                        info(iDS).seq_binSizeMs = seq{iDS}(1).params.dtMS;
                    end
                    info(iDS).seq_timeVector = seq{iDS}(1).y_time;
                end
            end

            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
            prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading posterior means for each dataset');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                if ~exist(fullfile(r.pathLFADSOutput, trainList{iDS}), 'file')
                    oldFile = strrep(trainList{iDS}, 'posterior_sample_and_average', 'posterior_sample');
                    if exist(fullfile(r.pathLFADSOutput, oldFile), 'file')
                        trainList{iDS} = oldFile;
                    else
    %                     warning('LFADS Posterior Mean train file not found for dataset %d: %s', ...
    %                         iDS, fullfile(r.pathLFADSOutput, trainList{iDS}));
                        pms = [];
                        return;
                    end
                end
                if ~exist(fullfile(r.pathLFADSOutput, validList{iDS}), 'file')
                    oldFile = strrep(validList{iDS}, 'posterior_sample_and_average', 'posterior_sample');
                    if exist(fullfile(r.pathLFADSOutput, oldFile), 'file')
                        validList{iDS} = oldFile;
                    else
    %                     warning('LFADS Posterior Mean valid file not found for dataset %d: %s', ...
    %                         iDS, fullfile(r.pathLFADSOutput, validList{iDS}));
                        pms = [];
                        return;
                    end
                end

                if isfield(info, 'conditionId')
                    conditionIds = info(iDS).conditionId;
                else
                    conditionIds = [];
                end
                if isfield(info, 'counts')
                    rawCounts = info(iDS).counts;
                else
                    rawCounts = [];
                end

                pmData = LFADS.Utils.loadPosteriorMeans(fullfile(r.pathLFADSOutput, validList{iDS}), ....
                    fullfile(r.pathLFADSOutput, trainList{iDS}), ...
                    info(iDS).validInds, info(iDS).trainInds);

                dt_pm = r.params.spikeBinMs;
                dt_y = info(iDS).seq_binSizeMs;
                rebin = dt_pm / dt_y;
                time = info(iDS).seq_timeVector(1:rebin:end);
                time = time(1:size(pmData.rates, 2));

                pms(iDS) = LFADS.PosteriorMeans(pmData, r.params, time, conditionIds, rawCounts); %#ok<AGROW>
            end
            prog.finish();

            pms = LFADS.Utils.makecol(pms);
            r.posteriorMeans = pms;
        end

        function tf = checkPosteriorMeansExist(r, verbose)
            if nargin < 2
                verbose = false;
            end
            tf = true;

            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
             for iDS = 1:r.nDatasets
                if ~exist(fullfile(r.pathLFADSOutput, trainList{iDS}), 'file')
                    tf = false;
                    if verbose
                        fprintf('LFADS train file not found %s\n', trainList{iDS});
                        continue;
                    else
                        return;
                    end
                end
                if ~exist(fullfile(r.pathLFADSOutput, validList{iDS}), 'file')
                    tf = false;
                    if verbose
                        fprintf('LFADS valid file not found %s\n', validList{iDS});
                        continue;
                    else
                        return;
                    end
                end
             end
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
                    seqs{iDS}(ntr).generator_ics = squeeze(pm.generator_ics(:,ntr));
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

        function mtp = loadModelTrainedParams(r, varargin)
            fname = r.fileModelParams;
            assert(exist(fname, 'file') > 1, 'model_params file not found. Ensure that runCommandLFADSWriteModelParams has been run');

            r.modelTrainedParams = LFADS.ModelTrainedParams(fname, r.datasetNames, varargin{:});
            mtp = r.modelTrainedParams;
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
        
        function fitLog = loadFitLog(r, varargin)
            fname = r.fileFitLog;
            assert(exist(fname, 'file') > 1, 'fitlog.csv file not found in lfadsOutput directory. Ensure that model has been trained');

            r.fitLog = LFADS.FitLog(fname, r.name);
            fitLog = r.fitLog;
        end
    end
    
    methods(Static)
        function str = generateCheckRunLFADSPyFoundString(run_lfads_path)
            if nargin < 1
                run_lfads_path = '$(which run_lfads.py)';
            end
            str = ['path_to_run_lfads=' run_lfads_path newline ...
                   'if [ ! -n "$path_to_run_lfads" ]; then' newline ...
                   '    echo "Error: run_lfads.py not found on PATH. Ensure you add LFADS to your system PATH."' newline ...
                   '    exit 1' newline ...
                   'fi' newline];
        end
    end
end
