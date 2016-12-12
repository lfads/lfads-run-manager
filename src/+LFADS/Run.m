classdef Run < handle & matlab.mixin.CustomDisplay
% Wraps a set of runs often sharing common parameter settings but utilizing
% different sets of datasets

    properties
        name = ''
        comment = ''
        
        version = 2;
    end
    
    properties(SetAccess=?LFADS.RunCollection)
        runCollection
        datasets
    end
    
    properties(Dependent)
        nDatasets
        datasetCollection
        path % unique folder within rootPath including name_paramSuffix
        params
        
        pathSequenceFiles
        sequenceFileNames
        
        pathLFADSInput
        lfadsInputFileNames
        
        pathLFADSOutput
        
        fileShellScriptLFADSTrain
        fileShellScriptLFADSPosteriorMeanSample
        
        nameWithParams
    end
    
    methods(Abstract)
        seq = convertDatasetToSequenceStruct(r, dataset, data);
        
        [seqData, alignmentMatrices, trainInds, validInds] = prepareSequenceDataForLFADS(seqData);
    end
    
    methods
        function r = Run(name, runCollection)
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
            r.datasets = r.datasetCollection.datasets(idx);
        end
        
        function selectDatasetsByName(r, names)
            r.datasets = r.datasetCollection.matchDatasetsByName(names);
        end
        
        function [trainList, validList] = getLFADSPosteriorSampleMeanFiles(r)
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
    
    methods
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
            r.makeSequenceFiles();
            r.makeLFADSInput();
            f = r.writeShellScriptLFADSTrain();
            fprintf('Shell script for training\n%s\n', f);
        end
        
        function makeSequenceFiles(r)
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
        
        function seq = loadSequenceFiles(r)
            % load seq files
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
        end
        
        function makeLFADSInput(r)
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
            [seqData, alignmentMatrices, trainInds, validInds] = r.prepareSequenceDataForLFADS(seqData);

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
           % function lfadsi_run_posterior_mean_sample(lfdir, varargin)
            %  calls the python script to sample the posterior and generate parameter estimates
            %  arguments
            %    lfdir: path to the trained network
            %
            %    varargin: overloads the parameters used to fit the model, e.g.
            %     (..., 'checkpoint_pb_load_name', 'checkpoint_lve',
            %     'batch_size', 512, 'device', NDEVICE,
            %     'ps_nexamples_to_process', num_trials_to_process)
            %
            %   checkpoint_pb_load_name selects which file to use for
            %   checkpoints
            %
            %   batch_size selects how many samples to use for computing
            %   posterior average
            %
            %   device specifies which gpu to run this on
            %
            %   defaults (overwrites hyperparams file with these):
            %       batch_size: 512
            %       checkpoint_pb_load_name: 'checkpoint_lve'

            lfdir = r.pathLFADSOutput;
           
            % set some default parameters
            if ~exist('varargin','var'), varargin ={}; end
            keys = varargin(1:2:end);

            default_keys = {'batch_size', 'checkpoint_pb_load_name'};
            default_values = {r.params.batchSize, 'checkpoint_lve'}; % was 512 not r.params.batchSize
            for nkey = 1:numel(default_keys)
                % if this key is not defined, set it to the default value
                if ~any(strcmp(keys, default_keys{nkey}))
                    keys{end+1} = default_keys{nkey}; %#ok<*AGROW>
                    keys{end+1} = default_values{nkey};
                end
            end

            params = lfadsi_read_parameters(lfdir); %#ok<*PROPLC>

            % make sure these are up to date
            params.data_dir = r.pathLFADSInput;
            params.lfads_save_dir = r.pathLFADSOutput;
            
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
            f = r.fileShellScriptLFADSPosteriorMeanSample;
            fid = fopen(f, 'w');
            fprintf(fid, r.buildCommandLFADSPosteriorMeanSample());
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end
            
        function pms = loadPosteriorMeanSamples(r)
            info = load(fullfile(r.path, 'lfadsInputInfo.mat'));
            [trainList, validList] = r.getLFADSPosteriorSampleMeanFiles();
            pms = cell(r.nDatasets, 1);
            
            for iDS = 1:r.nDatasets
                if ~exist(fullfile(r.pathLFADSOutput, trainList{iDS}), 'file')
                    error('LFADS Posterior Mean train file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, trainList{iDS}));
                end
                if ~exist(fullfile(r.pathLFADSOutput, validList{iDS}), 'file')
                    error('LFADS Posterior Mean valid file not found for dataset %d: %s', ...
                        iDS, fullfile(r.pathLFADSOutput, validList{iDS}));
                end
                pms{iDS} = LFADS.Utils.loadPosteriorMeans(fullfile(r.pathLFADSOutput, validList{iDS}), ....
                    fullfile(r.pathLFADSOutput, trainList{iDS}), ...
                    info.validInds{iDS}, info.trainInds{iDS});
                pms{iDS}.params = info.params;
            end
        end
    end
end
