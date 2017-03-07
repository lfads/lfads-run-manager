classdef Run < handle & matlab.mixin.CustomDisplay
    % Represents a single LFADS experiment on a specific set of datasets. Runs are grouped into :ref:`LFADS_RunCollection`
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.
    
    methods(Abstract)
        % These methods will need to be implemented in a subclass that provides the custom behavior for your application
        
        seq = convertDatasetToSequenceStruct(r, dataset, data);
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
        %   - `.params.dtMS` specifies the time bin width for `.y`
        %
        % Parameters
        % ------------
        % dataset : :ref:`LFADS_Dataset`
        %   The :ref:`LFADS_Dataset` instance from which data were loaded
        % data : typically struct array
        %   data contents of Dataset as returned by `loadData`
        %
        % Returns
        % ----------
        % seq : struct Array
        %   sequence formatted data. A struct array where each elemnt corresponds to a specific trial.
        
    end
    
    properties
        name char = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk
        
        comment char = '' % Textual comment for convenience
        
        version uint32 = 3; % Internal versioning allowing for graceful evolution of path settings
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
        
        pathSequenceFiles % Path on disk where sequence files will be saved
        sequenceFileNames % List of sequence file names (sans path)
        
        pathLFADSInput % Path on disk where LFADS input hd5 files will be saved
        lfadsInputFileNames % List of LFADS input hd5 files (sans path)

        lfadsInputInfoFileName % List of LFADS input info .mat files (sans path)
        
        pathLFADSOutput % Path on disk where LFADS output will be written
        
        fileShellScriptLFADSTrain % Location on disk where shell script for LFADS training will be written
        fileShellScriptLFADSPosteriorMeanSample  % Location on disk where shell script for LFADS posterior mean sampling will be written
        
        sessionNameTrain % name of tmux session that will be created if useSession = true is passed to writeShellScriptLFADSTrain
        sessionNamePosteriorMean % name of tmux session that will be created if useTmuxSession = true is passed to writeShellScriptLFADSPosteriorMean
    end
    
    methods
        function alignmentMatrices = prepareAlignmentMatrices(seqData)
            % Prepares alignment matrices to seed the stitching process when using multiple days of sequence data for
            % LFADS input file generation. Generate
            % alignment matrices which specify the initial guess at the encoder matrix that converts neural activity
            % from each dataset to a common set of factors (for stitching). Specify training and validation indices
            % (subsets of trials) on each day.
            %
            % Parameters
            % ------------
            % seqData : `nDatasets` cell of struct arrays of sequence data
            %   Sequence data for each dataset as returned by `convertDatasetToSequenceStruct`
            %
            % Returns
            % ----------
            % alignmentMatrices : `nDatasets` cell of `nNeurons` x `nFactors` matrices
            %   For each dataset, an initial guess at the encoder matrices which maps `nNeurons` (for that dataset) to a
            %   common set of `nFactors` (up to you to pick this). Seeding this well helps the stitching process. Typically,
            %   PC regression can provide a reasonable set of guesses.
            
            alignmentMatrices = {};
            error('You must override prepareAlignmentMatrices in your Run class if params.useAlignmentMatrices is set to true');
        end
        
        function r = Run(varargin)
            
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
        
        function p = get.pathSequenceFiles(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.path, 'seq');
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
        
        function name = get.lfadsInputInfoFileName(r)
            name = 'lfadsInputInfo.mat';
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
                header = sprintf('%s\n  Path: %s\n\n  %s "%s" : %s\n\n  %d datasets in "%s"\n', r.getFirstLineHeader(), r.path, ...
                    class(r.params), r.paramsString, r.params.generateShortDifferencesString, ...
                    r.nDatasets, r.datasetCollection.name);
                for s = 1:r.nDatasets
                    header = cat(2, header, sprintf('    [%2d] %s', s, r.datasets(s).getHeader()));
                end
            end
        end
    end
    
    methods
        function prepareForLFADS(r)
            % Generate all files needed to run LFADS.
            % Shortcut for generating and saving sequence files, LFADS input HD5 input files, and the shell script for running LFADS python code.
            r.makeSequenceFiles();
            r.makeLFADSInput();
            f = r.writeShellScriptLFADSTrain();
            fprintf('Shell script for training "%s": \n  %s\n', r.name, f);
        end
        
        function makeSequenceFiles(r)
            % Generate the seqence files and save them to disk
            if isempty(r.datasets)
                fprintf('No datasets added to Run\n');
                return;
            end
            
            LFADS.Utils.mkdirRecursive(r.pathSequenceFiles);
            
            sequenceFileNames = r.sequenceFileNames; %#ok<PROP>
            
            prog = LFADS.Utils.ProgressBar(numel(r.datasets), 'Generating sequence files');
            for iDS = 1:numel(r.datasets)
                ds = r.datasets(iDS);
                prog.update(iDS, 'Generating sequence files for %s', ds.name);
                
                % call user function
                seq = r.convertDatasetToSequenceStruct(ds);
                
                seqFile = fullfile(r.pathSequenceFiles, sequenceFileNames{iDS}); %#ok<PROP>
                save(seqFile, 'seq');
            end
            prog.finish();
        end
        
        function out = loadInputInfo(r)
            fname = fullfile(r.path, r.lfadsInputInfoFileName);
            r.inputInfo = load(fname);
            out = r.inputInfo;
        end

        function seq = loadSequenceData(r, reload)
            % seq = loadSequenceData([reload = True])
            % Load the sequence files from disk, caches them in
            % .sequenceData.
            %
            % Args:
            %   reload (bool) : Reload sequence data from disk even if already found in
            %     .sequenceData. Default = false
            %
            % Returns:
            %   seqData (cell of struct arrays) : nDatasets cell array of sequence structures loaded from sequence files on disk
            
            
            if nargin < 2
                reload = false;
            end
            if ~reload && ~isempty(r.sequenceData)
                seq = r.sequenceData;
                return;
            end
            
            seqFiles = cellfun(@(file) fullfile(r.pathSequenceFiles, file), r.sequenceFileNames, 'UniformOutput', false);
            
            prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading sequence files');
            seq = cell(r.nDatasets, 1);
            for nd = 1:r.nDatasets
                prog.update(nd);
                if ~exist(seqFiles{nd}, 'file')
                    error('Could not located sequence file for dataset %d: %s', nd, seqFiles{nd});
                end
                tmp = load(seqFiles{nd});
                seq{nd} = tmp.seq;
            end
            prog.finish();
            
            r.sequenceData = r.modifySequenceDataPostLoading(seq);
        end
        
        function seq = modifySequenceDataPostLoading(r, seq)
            % Optionally make any changes or do any post-processing of sequence data upon loading
            
        end
        
        function makeLFADSInput(r)
            % Generate the LFADS input HD5 files and save them to disk
            
            seqs = {};
            validInds = {};
            trainInds = {};
            
            par = r.params;
            
            % load sequence data
            seqFiles = cellfun(@(file) fullfile(r.pathSequenceFiles, file), r.sequenceFileNames, 'UniformOutput', false);
            
            % load seq files
            seqData = cell(r.nDatasets, 1);
            prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading sequence files');
            for nd = 1:r.nDatasets
                prog.update(nd);
                tmp = load(seqFiles{nd});
                seqData{nd} = tmp.seq;
            end
            prog.finish();
            
            % if there are multiple datasets, we need an alignment matrix
            if r.nDatasets > 1 && r.params.useAlignmentMatrix
                % call out to abstract dataset specific method
                useAlignMatrices = true;
                alignmentMatrices= r.prepareAlignmentMatrices(seqData);
            else
                useAlignMatrices = false;
            end
            
            % choose validation and training trial indices
            [validInds, trainInds] = deal(cell(r.nDatasets, 1));
            for nd = 1:r.nDatasets
                allInds = 1:numel(seqData{nd});
                validInds{nd} = 1 : (r.params.trainToTestRatio+1) : numel(seqData{nd});
                trainInds{nd} = setdiff(allInds, validInds{nd});
            end
            % arguments for the 'seq_to_lfads' call below
            seqToLFADSArgs = {'binSizeMs', par.spikeBinMs,  ...
                'inputBinSizeMs', seqData{1}(1).params.dtMS, ...
                'trainInds', trainInds, 'testInds', validInds};
            
            if useAlignMatrices
                seqToLFADSArgs{end+1} = 'alignment_matrix_cxf';
                seqToLFADSArgs{end+1} = alignmentMatrices;
            end
            
            % write the actual lfads input file
            LFADS.seq_to_lfads(seqData, r.pathLFADSInput, r.lfadsInputFileNames, ...
                seqToLFADSArgs{:});
            
            fname = fullfile(r.path, r.lfadsInputInfoFileName);
            params = r.params; %#ok<*NASGU,PROP>
            save(fname, 'trainInds', 'validInds', 'params');
        end
        
        function f = writeShellScriptLFADSTrain(r, varargin)
            % f = writeShellScriptLFADSTrain(cuda_visible_device, display, varargin)
            % Write a shell script used for running the LFADS python code
            %
            % Args:
            %   cuda_visible_devices : int
            %     which GPUs to make visible to CUDA, e.g. 0 or 1
            %   display : int
            %     which display to use for plot generation, e.g. 500 -->
            %     DISPLAY=:500
            %   useTmuxSession : bool
            %     if true, will prefix the command so that it runs within a
            %     new tmux session
            %
            % Returns:
            %   file : string
            %     Full path to shell script which can be used to begin or resume LFADS training
            
            p = inputParser();
            p.addOptional('cuda_visible_devices', [], @isscalar);
            p.addOptional('display', '', @(x) isnumeric(x) && mod(x,1)==0);
            p.addParameter('useTmuxSession', false, @islogical);
            p.addParameter('keepSessionAlive', true, @islogical);
            p.parse(varargin{:});
            
            f = r.fileShellScriptLFADSTrain;
            fid = fopen(f, 'w');
            outputString = r.buildLFADSTrainingCommand();
            
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
            
            % if requested, tmux-ify the command
            if p.Results.useTmuxSession
                outputString = LFADS.Utils.tmuxify_string( outputString, r.sessionNameTrain, 'keepSessionAlive', p.Results.keepSessionAlive);
                fprintf('Tmux Session is %s\n  tmux a -t %s\n\n', r.sessionNameTrain, r.sessionNameTrain);
            end
            
            fprintf(fid, outputString);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end
        
        function runLFADSTrainingCommand(r)
            %   function runLFADSTrainingCommand(r)
            system( sprintf('sh %s', r.fileShellScriptLFADSTrain) );
        end
        
        function outputString = buildLFADSTrainingCommand(r)
            outputString = sprintf(['python $(which run_lfads.py) --data_dir=%s --data_filename_stem=lfads ' ...
                '--lfads_save_dir=%s'], ...
                r.pathLFADSInput, r.pathLFADSOutput);
            
            % use the method from +LFADS/RunParams.m
            optionsString = r.params.generateCommandLineOptionsString(r);
            outputString = sprintf('%s%s', outputString, optionsString);
        end
        
        function cmd = buildCommandLFADSPosteriorMeanSample(r, inputParams, varargin)
            % Generates the command string for LFADS posterior mean sampling
            %
            % Returns
            % --------
            % cmd : string
            %   Shell command for running LFADS posterior mean sampling
            
            p = inputParser();
            p.addParameter('loadHyperparametersFromFile', false, @islogical);
            p.addParameter('batchSize', 512, @isscalar);
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
                if exist('inputParams','var')
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
        end
        
        function runLFADSPosteriorMeanCommand(r)
            %   function runLFADSPosteriorMeanCommand(r)
            system( sprintf('sh %s', r.fileShellScriptLFADSPosteriorMeanSample) );
        end
        
        function f = writeShellScriptLFADSPosteriorMeanSample(r, varargin)
            % Write a shell script used for running the LFADS posterior mean sampling. This must be run after the LFADS
            % training has been started.
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
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            f = r.fileShellScriptLFADSPosteriorMeanSample;
            fid = fopen(f, 'w');
            
            outputString = r.buildCommandLFADSPosteriorMeanSample(p.Unmatched);

            % set cuda visible devices
            if ~isempty(p.Results.cuda_visible_devices)
                outputString = sprintf('CUDA_VISIBLE_DEVICES=%i %s', ...
                    p.Results.cuda_visible_devices, outputString);
            end
            
            % if requested, tmux-ify the command
            if p.Results.useTmuxSession
                outputString = LFADS.Utils.tmuxify_string( outputString, r.sessionNamePosteriorMean, 'keepSessionAlive', p.Results.keepSessionAlive );
                fprintf('Tmux Session is %s\n  tmux a -t %s\n\n', r.sessionNamePosteriorMean, r.sessionNamePosteriorMean);
            end
            
            fprintf(fid, outputString);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
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
            
            info = load(fullfile(r.path, r.lfadsInputInfoFileName));
            
            if ~isequal(info.params, r.params)
                warning('Params saved for run in lfadsInputInfo.mat do not match. See params inside %s.', fullfile(r.path, r.lfadsInputInfoFileName));
            end
                
            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
            prog = LFADS.Utils.ProgressBar(r.nDatasets, 'Loading posterior means for each dataset');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                if ~exist(fullfile(r.pathLFADSOutput, trainList{iDS}), 'file')
                    warning('LFADS Posterior Mean train file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, trainList{iDS}));
                    continue;
                end
                if ~exist(fullfile(r.pathLFADSOutput, validList{iDS}), 'file')
                    warning('LFADS Posterior Mean valid file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, validList{iDS}));
                    continue;
                end
                
                pmData = LFADS.Utils.loadPosteriorMeans(fullfile(r.pathLFADSOutput, validList{iDS}), ....
                    fullfile(r.pathLFADSOutput, trainList{iDS}), ...
                    info.validInds{iDS}, info.trainInds{iDS});
                
                dt_y = seq{iDS}(1).params.dtMS;
                dt_pm = r.params.spikeBinMs;
                rebin = dt_pm / dt_y;
                time = seq{iDS}(1).y_time(1:rebin:end);
                
                pms(iDS) = LFADS.PosteriorMeans(pmData, info.params, time); %#ok<AGROW>
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
                    seqs{iDS}(ntr).controller_outputs = ...
                        squeeze(pm.controller_outputs(:,:,ntr));
                end
            end
            
            r.sequenceData = seqs;
        end
    end
end
