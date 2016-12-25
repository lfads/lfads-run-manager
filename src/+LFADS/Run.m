classdef Run < handle & matlab.mixin.CustomDisplay
    % Represents a single LFADS experiment on a specific set of datasets. Runs are grouped into :ref:`LFADS_RunCollection`
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.

    methods(Abstract)
        % These methods will need to be implemented in a subclass that provides the custom behavior for your application

        seq = convertDatasetToSequenceStruct(r, dataset, data);
        % Converts the loaded data within a dataset into a sequence struct. The sequence data returned will be a
        % struct array where each element corresponds to a trial. You can include any metadata or information fields
        % that you like for future reference or analysis. At a minimum, you must include field `.y`. For each trial,
        % `.y` will contain spike data as an `nNeurons` x `nTimeBins` array of spike counts. The spike binning is
        % typically specified in this Run's `.params`  in the field `spikeBinMs`.
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

        [alignmentMatrices, trainInds, validInds] = prepareSequenceDataForLFADS(seqData);
        % Prepares multiple days of sequence data for LFADS input file generation. Generate
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
        %   PC regression can provide a reasonable set of guesses
        %
        % trainInds : `nDatasets` cell of indices into each dataset's trial array
        %   Trial indices for each datatset to use for training
        % validInds : `nDatasets` cell of indices into each dataset's trial array
        %   Trial indices for each datatset to use for validation
    end

    properties
        name = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk

        comment = '' % Textual comment for convenience

        version = 2; % Internal versioning allowing for graceful evolution of path settings
    end

    properties
        runCollection % :ref:`LFADS_RunCollection` instance to which this run belongs
        datasets % Array of :ref:`LFADS_Dataset` instances which this particular Run will utilize
        
        sequenceData % nDatasets cell array of sequence struct data
        posteriorMeans % nDatasets array of :ref:`LFADS_PosteriorMeans` when loaded
    end

    properties(Dependent)
        nDatasets % Number of datasets used by this run
        datasetCollection % Dataset collection used by this run (and all runs in the same RunCollection)
        path % Unique folder within rootPath including name_paramSuffix
        params % :ref:`LFADS_RunParams` instance shared by all runs in the collection, contains parameter settings

        pathSequenceFiles % Path on disk where sequence files will be saved
        sequenceFileNames % List of sequence file names (sans path)

        pathLFADSInput % Path on disk where LFADS input hd5 files will be saved
        lfadsInputFileNames % List of LFADS input hd5 files (sans path)

        pathLFADSOutput % Path on disk where LFADS output will be written

        fileShellScriptLFADSTrain % Location on disk where shell script for LFADS training will be written
        fileShellScriptLFADSPosteriorMeanSample  % Location on disk where shell script for LFADS posterior mean sampling will be written

        nameWithParams % Concatenation of this Run's name and the serialized representation of the RunParams
    end

    methods
        function r = Run(name, runCollection)
            % run = Run(name, runCollection)
            %
            % Parameters
            % ------------
            % name : string
            %   Unique name for this run within the collection
            %
            % runCollection : :ref:`LFADS_RunCollection` instance
            %   RunCollection to which this run should be added

            r.name = name;
            runCollection.addRun(r);
        end

        function p = get.params(r)
            p = r.runCollection.params;
        end

        function dc = get.datasetCollection(r)
            dc = r.runCollection.datasetCollection;
        end

        function p = get.path(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.runCollection.path, r.name);
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

        function n = get.nameWithParams(r)
            if isempty(r.params)
                paramStr = '';
            else
                paramStr = ['_', r.params.generateSuffix()];
            end
            n = [r.name, paramStr];
        end

        function selectDatasetsByIndex(r, idx)
            % Specify the datasets that this Run should use from the set of datasets in its DatasetCollection by indices
            %
            % Parameters
            % ------------
            % idx : logical mask or indices
            %   Selection applied to this run's DatasetCollection's array of datasets

            r.datasets = r.datasetCollection.datasets(idx);
        end

        function selectDatasetsByName(r, names)
            % Specify the datasets that this Run should use from the set of datasets in its DatasetCollection using name matching
            %
            % Parameters
            % ------------
            % names : string or cellstr
            %   Name or names to search for within this run's DatasetCollection's array of datasets
            r.datasets = r.datasetCollection.matchDatasetsByName(names);
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
    end

    methods(Hidden)
        function h = getFirstLineHeader(r)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(r);
            h = sprintf('%s %s (%d datasets)\n', className, r.name, r.nDatasets);
        end
    end

    methods (Access = protected)
       function header = getHeader(r)
          if ~isscalar(r)
             header = getHeader@matlab.mixin.CustomDisplay(r);
          else
             rc = r.runCollection;
             header = sprintf('%s\n  %d datasets in %s\n', r.getFirstLineHeader(), r.nDatasets, rc.path);
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
            fprintf('Shell script for training\n%s\n', f);
        end

        function makeSequenceFiles(r)
            % Generate the seqence files and save them to disk
            if isempty(r.datasets)
                fprintf('No datasets added to Run\n');
                return;
            end

            fprintf('Saving sequence files in %s\n', r.pathSequenceFiles);
            LFADS.Utils.mkdirRecursive(r.pathSequenceFiles);

            sequenceFileNames = r.sequenceFileNames; %#ok<PROP>

%             prog = LFADS.Utils.ProgressBar(numel(r.datasets), 'Generating sequence files');
            for iDS = 1:numel(r.datasets)
                ds = r.datasets(iDS);
%                 prog.update(iDS, 'Generating sequence files for %s', ds.name);

                % call user
                seq = r.convertDatasetToSequenceStruct(ds);

                seqFile = fullfile(r.pathSequenceFiles, sequenceFileNames{iDS}); %#ok<PROP>
                fprintf('Saving seq file to %s\n', seqFile);
                save(seqFile, 'seq');
            end
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
            
            r.sequenceData = seq;
        end

        function makeLFADSInput(r)
            % Generate the LFADS input HD5 files and save them to disk

            seqs = {};
            validInds = {};
            trainInds = {};

            par = r.params;

            dtMS = par.spikeBinMs; % this is bin size into LFADS, want to make larger for speedy training times (was 2 orig)
            pcs_to_keep = par.pcsKeep;
            do_trial_averaging = par.pcTrialAvg;

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

            % call out to abstract dataset specific method
            [alignmentMatrices, trainInds, validInds] = r.prepareSequenceDataForLFADS(seqData);

            LFADS.seq_to_lfads(seqData, r.pathLFADSInput, r.lfadsInputFileNames, ...
                         'trainInds',trainInds, 'testInds', validInds, ...
                         'binSizeMs', dtMS, ...
                         'inputBinSizeMs', seqData{1}(1).params.dtMS, ...
                         'alignment_matrix_cxf', alignmentMatrices);

            fname = fullfile(r.path, 'lfadsInputInfo.mat');
            params = r.params; %#ok<*NASGU,PROP>
            save(fname, 'trainInds', 'validInds', 'params');
        end

        function f = writeShellScriptLFADSTrain(r)
            % file = writeShellScriptLFADSTrain()
            % Write a shell script used for running the LFADS python code
            %
            % Returns
            % --------
            % file : string
            %   Full path to shell script which can be used to begin or resume LFADS training

            f = r.fileShellScriptLFADSTrain;
            fid = fopen(f, 'w');
            fprintf(fid, 'SESS_NAME="%s"\n', r.nameWithParams);
            fprintf(fid, 'L2_GEN_SCALE="500"\n');
            fprintf(fid, 'L2_CON_SCALE="500"\n');
            fprintf(fid, 'DATADIR="%s"\n', r.pathLFADSInput);
            fprintf(fid, 'OUTDIR="%s"\n', r.pathLFADSOutput);
            fprintf(fid, 'python $(which run_lfads.py) --data_dir=$DATADIR --data_filename_stem=lfads --lfads_save_dir=$OUTDIR --cell_clip_value=5 --factors_dim=8 --in_factors_dim=8 --ic_enc_dim=100 --ci_enc_dim=100 --gen_dim=100 --keep_prob=%g --learning_rate_decay_factor=%g --device=/gpu:0 --co_dim=4 --do_causal_controller=false --l2_gen_scale=$L2_GEN_SCALE --l2_con_scale=$L2_CON_SCALE --batch_size=%d --kl_increase_steps=%d --l2_increase_steps=%d --controller_input_lag=1 \n', ...
                r.params.keepProb, r.params.learningRateDecayFactor, r.params.batchSize, r.params.regularizerIncreaseSteps, r.params.regularizerIncreaseSteps);
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end

        function cmd = buildCommandLFADSPosteriorMeanSample(r, varargin)
            % Generates the command string for LFADS posterior mean sampling
            %
            % Returns
            % --------
            % cmd : string
            %   Shell command for running LFADS posterior mean sampling

            lfdir = r.pathLFADSOutput;

            % set some default parameters
            %if ~exist('varargin','var'), varargin ={}; end
            %keys = varargin(1:2:end);

            params = lfadsi_read_parameters(lfdir); %#ok<*PROPLC>

            % make sure these are up to date
            params.data_dir = r.pathLFADSInput;
            params.lfads_save_dir = r.pathLFADSOutput;

            params.batch_size = r.params.batchSize;
            params.checkpoint_pb_load_name = 'checkpoint_lve';
%             default_keys = {'batch_size', 'checkpoint_pb_load_name'};
%             default_values = {r.params.batchSize, 'checkpoint_lve'}; % was 512 not r.params.batchSize
%             for nkey = 1:numel(default_keys)
%                 % if this key is not defined, set it to the default value
%                 if ~any(strcmp(keys, default_keys{nkey}))
%                     keys{end+1} = default_keys{nkey}; %#ok<*AGROW>
%                     keys{end+1} = default_values{nkey};
%                 end
%             end
%             
            % need to remove "dataset_names" and "dataset_dims"
            params = rmfield(params, {'dataset_names', 'dataset_dims'});
            use_controller = boolean(params.ci_enc_dim);

%             tmp = lfads_path();
%             if ~use_controller
%                 lfpath = fullfile(tmp.path, tmp.lfnc);
%             else
%                 lfpath = fullfile(tmp.path, tmp.lf);
%             end

            execstr = 'python';

            % take the arguments passed in via varargin, add new params, or overwrite existing ones
            for nn = 1:2:numel(varargin)
                % have to specially handle CUDA_VISIBLE_DEVICES as an environment ...
                %   variable rather than a command line param
                if strcmpi(varargin{nn},'device')
                    execstr = strcat(sprintf('CUDA_VISIBLE_DEVICES=%i', varargin{nn+1}),  ...
                                     execstr);
                else
                    params.(varargin{nn}) = varargin{nn+1};
                end
            end

            f = fields(params);
            optstr = '';
            for nf = 1:numel(f)
                fval = params.(f{nf});
                %convert any numbers to strings
                if isnumeric(fval), fval = num2str(fval); end
                optstr = strcat(optstr, sprintf(' --%s=%s',f{nf}, fval));
            end

            optstr = strcat(optstr, sprintf(' --kind=posterior_sample'));

            % put the command together
            cmd = sprintf('%s $(which run_lfads.py) %s', execstr, optstr);
        end

        function f = writeShellScriptLFADSPosteriorMeanSample(r)
            % file = writeShellScriptLFADSPosteriorMeanSample()
            % Write a shell script used for running the LFADS posterior mean sampling. This must be run after the LFADS
            % training has been started.
            %
            % Returns
            % --------
            % file : string
            %   Full path to shell script which can be used to perform LFADS posterior mean sampling on the lowest
            %   validation error checkpoint

            f = r.fileShellScriptLFADSPosteriorMeanSample;
            fid = fopen(f, 'w');
            fprintf(fid, r.buildCommandLFADSPosteriorMeanSample());
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end

        function pms = loadPosteriorMeans(r, reload)
            % pmData = loadPosteriorMeans(r)
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
                pms = r.posteriorMeans;
                return;
            end
                
            info = load(fullfile(r.path, 'lfadsInputInfo.mat'));
            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
            prog = ProgressBar(r.nDatasets, 'Loading posterior means for each dataset');
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
                pms(iDS) = LFADS.PosteriorMeans(pmData, info.params);
            end
            prog.finish();
            
            r.posteriorMeans = pms;
        end
    end
end
