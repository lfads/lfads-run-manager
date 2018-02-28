classdef RunCollection < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
    % A set of :ref:`LFADS_Run` instances utilizing
    % different subsets of :ref:`LFADS_Dataset` insstances in a :ref:`LFADS_DatasetCollection`.
    % A RunCollection is a logical grouping of LFADS runs, where multiple
    % runs specifications are each run with multiple parameter
    % settings. The run specifications are provided as arrays of RunSpec instances which
    % indicate the name and which datasets are included for each run specification. The
    % parameter sweeps are specified using arrays of RunParams instances.
    % From these, a `nRunSpecs` x `nParams` matrix of :ref:`LFADS_Run`
    % instances will be generated. These in turn can be used to train
    % individual LFADS models on a particular set of datasets (or single
    % dataset) and particular parameter settings.

    properties
        name char = '' % Name of this RunCollection, determines its relative path on disk
        comment char = '' % Textual comment for convenience
        rootPath char = ''; % Root path on disk under which individual Runs will be stored

        version uint32 = 20171107; % version used for graceful evolution with backwards compatibility

        datasetCollection % DatasetCollection instance
    end

    properties(SetAccess=protected)
        runs % `nRunSpecs` x `nParams` matrix of :ref:`LFADS_Run` instances
        params % array of RunParams instances
        runSpecs % array of :ref:`LFADS_RunSpec` instances used to specify the runs for each param setting
    end

    properties(Dependent)
        nParams % number of parameter settings in `params`
        nRunSpecs % number of run specifications
        nRunsTotal % number of runs total (equal to `nParams` * `nRunsEachParam`

        nDatasets % number of datasets within the datasetCollection
        datasetNames % nDatasets x 1 cell array of dataset names
        
        path % unique folder given by rootPath/name
        pathsCommonDataForParams % cellstr of folders given by rootPath/name/{data_HASH} for each params
        pathsForParams % cellstr of folders given by rootPath/name/{paramStr}
        fileShellScriptTensorboard % path the location where a shell script to launch TensorBoard for all runs will be written
        fileSummaryText % path where summary text info will be written

        fileShellScriptRunQueue % path of shell script to launch lfadsqueue to train and sample all runs within
    end

    methods
        function rc = RunCollection(varargin)
            % rc = RunCollection(rootPath, name, datasetCollection[, runParams, runSpecs])
            %
            % Args:
            %   rootPath : string
            %     location on disk under which run collections will be saved
            %   name : string
            %     unique identifier for this RunCollection, also used as subfolder within `rootPath`
            %   datasetCollection : :ref:`LFADS_DatasetCollection`
            %     DatasetCollection to be used by runs within this RunCollection
            %   runParams : array of :ref:`LFADS_RunParams` instances
            %     Parameter settings which will be used uniformly by all runs, will also be serialized in the folder path
            %   runSpecs : array of :ref:`LFADS_RunSpec` instances
            %     RunSpec instances describing the names and lists of
            %     datasets for each run that will be run with each set of
            %     params

            p = inputParser();
            p.addOptional('rootPath', '', @ischar);
            p.addOptional('name', '', @ischar);
            p.addOptional('datasetCollection', [], @(x) isa(x, 'LFADS.DatasetCollection'));
            p.addOptional('runParams', [], @(x) isa(x, 'LFADS.RunParams'));
            p.addOptional('runSpecs', [], @(x) isa(x, 'LFADS.DatasetCollection'));
            p.parse(varargin{:});

            rc.rootPath = p.Results.rootPath;
            rc.name = p.Results.name;
            rc.datasetCollection = p.Results.datasetCollection;
            rc.params = p.Results.runParams;
            rc.runSpecs = p.Results.runSpecs;
        end

        function addRunSpec(rc, runSpecs)
            % Adds new RunSpec instance(s) to this RunCollection.
            % Automatically appends new runs to `.runs`.
            %
            % Args:
            %   runSpec : LFADS.RunSpec
            %     new RunSpec instances describing the names and datasets
            %     included. Each of these new RunSpecs will
            %     be run on each of the parameter settings in `.params`

            assert(isa(runSpecs, 'LFADS.RunSpec'), 'Must be LFADS.RunSpec instance(s)');
            runSpecs = LFADS.Utils.makecol(runSpecs(:));

            for i = 1:numel(runSpecs)
                assert(isequal(runSpecs(i).datasetCollection, rc.datasetCollection), 'DatasetCollection of added RunSpecs must match RunCollection');
            end

            if isempty(rc.runSpecs)
                rc.runSpecs = runSpecs;
            else
                % check for existing runs by name and replace them
                [tf, idx] = rc.ismemberRunSpecs(runSpecs);
                if any(tf)
                    warning('Replacing existing runs with matching name');
                    rc.runSpecs(idx(tf)) = runSpecs(tf);
                end
                rc.runSpecs = cat(1, rc.runSpecs, runSpecs(~tf));
            end

            rc.generateRuns();
        end

        function clearRunSpecs(rc)
            % Flush list of run specs

            rc.runSpecs = [];
            rc.generateRuns();
        end

        function filterRunSpecs(rc, mask)
            % Apply selection mask to list of runs specs
            %
            % Parameters:
            %   mask : logical or indices
            %     selection mask applied to `.runSpecs`

            rc.runSpecs = rc.runSpecs(mask);
            rc.generateRuns();
        end

        function addParams(rc, params)
            % Adds new LFADS.RunParams instance(s) to this RunCollection.
            % Automatically appends new runs to `.runs`.
            %
            % Args:
            %   params : :ref:`LFADS_RunParams`
            %     array of parameter settings to run each RunSpec on

            assert(isa(params, 'LFADS.RunParams'), 'Must be LFADS.RunParams instance(s)');
            params = LFADS.Utils.makecol(params(:));
            if isempty(rc.params)
                rc.params = params;
            else
                for iP = 1:numel(params)
                    params.version = rc.version;
                end
                
                % check for existing runs by name and replace them
                [tf, idx] = rc.ismemberParams(params);
                if any(tf)
                    warning('Replacing existing runs with matching name');
                    rc.params(idx(tf)) = params(tf);
                end
                rc.params = cat(1, rc.params, params(~tf));
            end

            if rc.version < 3 && numel(rc.params) > 1
                error('Multiple params not supported for version < 3');
            end

            rc.generateRuns();
        end

        function clearParams(rc)
            % Flush list of run params

            rc.params = [];
            rc.generateRuns();
        end

        function filterParams(rc, mask)
            % Apply selection mask to .params
            %
            % Parameters:
            %   mask : logical or indices
            %     selection mask applied to `.runs`

            rc.params = rc.params(mask);
            rc.generateRuns();
        end

        function clearAll(rc)
            rc.runSpecs = [];
            rc.params = [];
            rc.generateRuns();
        end

        function generateRuns(rc)
            assert(~isempty(rc.datasetCollection), 'DatasetCollection is empty');
            assert(~isempty(rc.rootPath), 'rootPath is empty');
            assert(~isempty(rc.name), 'name is empty');

            if rc.nRunsTotal > 0
                % start with existing matrix
                matOrig = rc.runs;

                for iS = rc.nRunSpecs:-1:1
                    spec = rc.runSpecs(iS);
                    if ischar(spec.runClassName)
                        clsFn = str2func(spec.runClassName);
                    elseif isa(clsFn, 'function_handle')
                        % okay as is
                    else
                        error('Unknown runClassName, must be string or function_handle to constructor');
                    end

                    for iP = rc.nParams:-1:1
                        % create the new run by manually assigning
                        % properties. This avoids the need for the user to
                        % preserve the constructor
                        new = clsFn();
                        new.name = spec.name;
                        new.runCollection = rc;
                        new.params = rc.params(iP);
                        new.paramIndexInRunCollection = iP;
                        new.datasets = spec.datasets;
                        
                        % set the Run version to match RunCollectionVersion
                        new.version = rc.version;

                        % check whether the old run matches, keep it if so,
                        % so that we don't break references unless
                        % something has changed
                        if size(matOrig, 1) >= iS && size(matOrig, 2) >= iP
                            orig = matOrig(iS, iP);
                            if orig == new
                                new = orig;
                            end
                        end

                        mat(iS, iP) = new;
                    end
                end

                rc.runs = mat;
            else
                rc.runs = [];
            end
        end

        function prepareForLFADS(rc, regenerate)
            if nargin < 2
                regenerate = false;
            end
            
            if rc.nRunsTotal == 0
                if rc.nParams == 0
                    error('No RunParams have been added. Use rc.addParams');
                else
                    error('No RunSpecs have been added. Use rc.addRunSpec');
                end
            end
            
            prog = LFADS.Utils.ProgressBar(rc.nRunsTotal, 'Generating input data for each run');
            for iR = 1:rc.nRunsTotal
                prog.update(iR);
                rc.runs(iR).prepareForLFADS(regenerate);
            end
            prog.finish();

            rc.writeSummaryText();
        end

        function deleteLFADSOutput(rc, varargin)
            resp = input('Are you sure you want to delete EVERYTHING (type yes): ', 's');
            if ~strcmpi(resp, 'yes')
                return;
            end

            for iR = 1:rc.nRunsTotal
                rc.runs(iR).deleteLFADSOutput('confirmed', true);
            end
        end

        function out_file = writeShellScriptRunQueue(rc, varargin)
            p = inputParser();
            p.addParameter('rerun', false, @islogical);
            p.addParameter('gpuList', [], @(x) isempty(x) || isvector(x));
            p.addParameter('display', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('oneTaskPerGPU', true, @islogical); % false allows GPUs to multi-task, which seems to be slower than just running tasks sequentially
            p.addParameter('maxTasksSimultaneously', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('gpuMemoryRequired', 2000, @isscalar); % in MB
            
            p.addParameter('runIdx', 1:rc.nRunsTotal, @isvector); % subsets and orders the runs to include in the script
            
            p.addParameter('prependPathToLFADS', false, @islogical); % prepend an export path to run_lfads.py
            p.addParameter('virtualenv', '', @ischar); % prepend source activate environment name
            
            p.parse(varargin{:});

            rc.writeTensorboardShellScript();

            if isempty(p.Results.display)
                display = LFADS.Utils.getDisplay();
                if isempty(display)
                    error('Please specify a display');
                end
            else
                display = p.Results.display;
            end

            out_file = rc.fileShellScriptRunQueue;

            fid = fopen(out_file, 'w');
            
            if p.Results.prependPathToLFADS
                folder = LFADS.Utils.find_lfadsqueue_py();
                if ~isempty(folder)
                    fprintf(fid, 'import sys\n');
                    fprintf(fid, 'sys.path.append("%s")\n', folder);
                end
            end
            
            fprintf(fid, 'import lfadsqueue as lq\n\n');
            fprintf(fid, 'queue_name = "%s"\n', rc.name);
            fprintf(fid, 'tensorboard_script = "%s"\n', LFADS.Utils.GetFullPath(rc.fileShellScriptTensorboard));

            if ~isempty(p.Results.gpuList)
                gpuListStr = ['[', LFADS.Utils.strjoin(p.Results.gpuList, ', '), ']'];
                fprintf(fid, 'gpu_list = %s\n\n', gpuListStr);
            else
                fprintf(fid, 'gpu_list = []\n\n');
            end

            fprintf(fid, 'task_specs = [');
            
            % determine full ordering
            runIdx = p.Results.runIdx;
            assert(all(ismember(runIdx, 1:rc.nRunsTotal)), 'runIdx must be in 1:nRuns');
            
            nRunsIncluded = numel(runIdx);
            prog = LFADS.Utils.ProgressBar(nRunsIncluded, 'Writing shell scripts for each run');
            
            for iiR = 1:nRunsIncluded
                prog.update(iiR);
                iR = runIdx(iiR);
                rc.runs(iR).writeShellScriptLFADSTrain('display', display, 'useTmuxSession', false, ...
                    'appendPosteriorMeanSample', true, ...
                    'appendWriteModelParams', true, ...
                    'teeOutput', true, ...
                    'prependPathToLFADS', p.Results.prependPathToLFADS, ...
                    'virtualenv', p.Results.virtualenv);
                
                outfile = LFADS.Utils.GetFullPath(rc.runs(iR).fileLFADSOutput);
                donefile = LFADS.Utils.GetFullPath(fullfile(rc.runs(iR).path, 'lfads.done'));

                name = sprintf('lfads_%s_run%03d_%s', rc.runs(iR).paramsString, iR, rc.runs(iR).name); %#ok<*PROPLC>
                name = strrep(name, '.', '_');
                fprintf(fid, '{"name": "%s", "command": "bash %s", "memory_req": %d, "outfile": "%s", "donefile": "%s"}, \n', ...
                    name, LFADS.Utils.GetFullPath(rc.runs(iR).fileShellScriptLFADSTrain), ...
                    p.Results.gpuMemoryRequired, outfile, donefile);
            end
            prog.finish();
            fprintf(fid, ']\n\n');

            if p.Results.rerun
                donefileStr = ', ignore_donefile=True';
            else
                donefileStr = '';
            end

            if ~isempty(p.Results.maxTasksSimultaneously)
                maxStr = sprintf(', max_tasks_simultaneously=%d', p.Results.maxTasksSimultaneously);
            else
                maxStr = '';
            end
            
            if p.Results.oneTaskPerGPU
                oneTaskStr = ', one_task_per_gpu=True';
            else
                oneTaskStr = ', one_task_per_gpu=False';
            end
            fprintf(fid, 'tasks = lq.run_lfads_queue(queue_name, tensorboard_script, task_specs, gpu_list=gpu_list%s%s%s)\n\n', maxStr, donefileStr, oneTaskStr);
            fclose(fid);
        end
    end

    methods
        function p = get.path(rc)
            if rc.version >= 3
                p = fullfile(rc.rootPath, rc.name);
            else
                if isempty(rc.params)
                    paramSuffix = '';
                else
                    paramSuffix = ['_' rc.params(1).generateString()];
                end
                p = fullfile(rc.rootPath, [rc.name paramSuffix]);
            end
        end

        function pcell = get.pathsForParams(rc)
            p = rc.path;
            pcell = arrayfun(@(par) fullfile(p, par.generateHashName()), rc.params, 'UniformOutput', false);
        end

        function pcell = get.pathsCommonDataForParams(rc)
            p = rc.path;
            pcell = arrayfun(@(par) fullfile(p, par.generateInputDataHashName()), rc.params, 'UniformOutput', false);
        end

        function f = get.fileShellScriptTensorboard(r)
            f = fullfile(r.path, 'launch_tensorboard.sh');
        end

        function runTensorboard(r)
            system( sprintf('sh %s', r.fileShellScriptTensorboard) );
        end

        function f = get.fileSummaryText(r)
            f = fullfile(r.path, 'summary.txt');
        end

        function f = get.fileShellScriptRunQueue(rc)
            f = fullfile(rc.path, 'run_lfadsqueue.py');
        end

        function n = get.nParams(rc)
            n = numel(rc.params);
        end

        function n = get.nRunsTotal(rc)
            n = rc.nRunSpecs * rc.nParams;
        end

        function n = get.nRunSpecs(rc)
            n = numel(rc.runSpecs);
        end

        function n = get.nDatasets(rc)
            n = rc.datasetCollection.nDatasets;
        end
        
        function names = get.datasetNames(r)
            names = {r.datasetCollection.datasets.name}';
        end

        function [tf, idx] = ismemberRunSpecs(rc, runSpecSearch)
            % Determine if any run specs that match runSpecSearch are found in .runSpecs.
            % If runSpecSearch is string or cellstr, matches by
            % name. If runSpecSearch is an array of RunSpec instances,
            % finds them using isequal. If runSpecSearch is a vector of
            % indices, selects from runSpecs directly
            %
            % Args:
            %   runSpecSearch : array of LFADS.RunSpec, strings, or indices into .runSpecs
            %
            % Returns:
            %   tf : logical
            %     does each runSpec exist within .runSpecs
            %   idx : indices
            %     which index in .runSpecs
            %

            if ischar(runSpecSearch)
                runSpecSearch = {runSpecSearch};
            end
            if iscellstr(runSpecSearch)
                [tf, idx] = ismember(runSpecSearch, {rc.runSpecs.name});

            elseif isa(runSpecSearch, 'LFADS.RunSpec')
                [tf, idx] = ismember(runSpecSearch, rc.runSpecs);

            else
                % assume is selection
                idx = runSpecSearch;
                tf = true(size(idx));
            end
            idx = LFADS.Utils.makecol(idx(:));
        end

        function [runSpecs, idx] = findRunSpecs(rc, runSpecSearch)
            % Find run specs that match runSpecSearch. Throws an error if any run is not
            % found. If runSpecSearch is string or cellstr, matches by
            % name. If runSpecSearch is an array of RunSpec instances,
            % finds them using isequal. If runSpecSearch is a vector of
            % indices, selects from runSpecs directly.
            %
            % Args:
            %   runSpecSearch : array of LFADS.RunSpec, strings, or indices into .runSpecs
            %
            % Returns:
            %   runs : LFADS.Run
            %     matching runs
            %   idx : vector of indices
            %     vector of indices into `.runSpecs` of matching run
            %     instances

            [tf, idx] = rc.ismemberRunSpecs(runSpecSearch);
            assert(all(tf), 'Some run spec names could not be found in this RunCollection');
            runSpecs = rc.runSpecs(idx);
        end

        function [tf, idx] = ismemberParams(rc, paramSearch)
            % Determine if any run params that match paramSearch are found in .params.
            % If paramSearch is an array of RunParams instances,
            % finds them using isequal. If paramSearch is a vector of
            % indices, selects from runParams directly.
            %
            % Args:
            %   paramSearch : array of LFADS.RunParams, or indices into
            %   .params, or string or cell of strings like param_HASH or
            %   HASH
            %
            % Returns:
            %   params : LFADS.RunParams
            %     matching params
            %   idx : vector of indices
            %     selection into `.params` of matches

            if isa(paramSearch, 'LFADS.RunParams')
                if ~strcmp(class(paramSearch), class(rc.params))
                    % if classes don't match, cannot be member
                    tf = false(size(rc.params));
                    idx = [];
                else
                    [tf, idx] = ismember(paramSearch, rc.params);
                end

            elseif ischar(paramSearch) || iscellstr(paramSearch)
                if ischar(paramSearch)
                    paramSearch = {paramSearch}; 
                end
                
                for i = 1:numel(paramSearch)
                    if ~strncmp(paramSearch{i}, 'param_', 6)
                        paramSearch{i} = ['param_' paramSearch{i}];
                    end
                end
                [tf, idx] = ismember(paramSearch, {rc.params.paramHashString});
                
            else
                % assume is selection
                idx = paramSearch;
                tf = true(size(idx));
            end
            idx = LFADS.Utils.makecol(idx(:));
        end

        function [params, idx] = findParams(rc, paramSearch)
            % Find run params that match paramSearch. Throws an error if any RunParams is not
            % found. If paramSearch is an array of RunParams instances,
            % finds them using isequal. If paramSearch is a vector of
            % indices, selects from runParams directly
            %
            % Args:
            %   paramSearch : array of LFADS.RunParams or indices into .params
            %
            % Returns:
            %   params : LFADS.RunParams
            %     matching params
            %   idx : vector of indices
            %     selection into `.params` of matches

            [tf, idx] = rc.ismemberParams(paramSearch);
            assert(all(tf), 'Some run params could not be found in this RunCollection');
            params = rc.params(idx);
        end

        function [runs, rowIdx, colIdx] = findRuns(rc, runSpecSearch, paramSearch)
            % Returns a subset of the .runs matrix for a certain subset of RunSpecs and RunParams.
            % Throws an error if any run is not found.
            %
            % Args:
            %   runSpecSearch : array of LFADS.RunSpec, strings, or indices into .runSpecs
            %   paramSearch : array of LFADS.RunParams, or indices into .params
            %
            % Returns:
            %   runs : LFADS.Run matrix
            %     matching runs as matrix
            %   rowIdx : vector of indices
            %     vector of indices matching into `.runSpecs`
            %   colIdx : vector of indices
            %     vector of indices matching into .params

            if isempty(runSpecSearch)
                rowIdx = 1:rc.nRunSpecs;
            else
                [~, rowIdx] = rc.findRunSpecs(runSpecSearch);
            end
            if nargin > 2 && ~isempty(paramSearch)
                [~, colIdx] = rc.findParams(paramSearch);
            else
                colIdx = 1:rc.nParams;
            end

            runs = rc.runs(rowIdx, colIdx);
        end

        function str = getTensorboardCommand(rc, varargin)
            % Generates the shell command text to launch a TensorBoard displaying all runs within this collection
            %
            % Returns:
            %   cmd : string
            %     Command which can luanch TensorBoard from command line

            ip = inputParser();
            ip.addOptional('port', [], @(x) isscalar(x) || isempty(x));
            ip.addParameter('useTmuxSession', false, @islogical);
            ip.addParameter('tmuxName', rc.name, @ischar);
            ip.parse(varargin{:});

            runEntry = cell(rc.nRunSpecs, rc.nParams);
            for s = 1:rc.nRunSpecs
                for p = 1:rc.nParams
                    runEntry{s,p} = sprintf('%s/%s:"%s"', ...
                        rc.runs(s, p).params.generateHashName(), ...
                        rc.runs(s,p).name, ...
                        LFADS.Utils.GetFullPath(rc.runs(s,p).pathLFADSOutput));
                end
            end

            if ~isempty(ip.Results.port)
                portStr = sprintf('--port=%i', ip.Results.port);
            else
                portStr = '';
            end

            str = sprintf('tensorboard --logdir=%s %s', LFADS.Utils.strjoin(runEntry, ','), portStr);

            if ip.Results.useTmuxSession
                str = LFADS.Utils.tmuxify_string(str, ip.Results.tmuxName);
            end
        end

        function f = writeTensorboardShellScript(rc, varargin)
            % Generates the shell command text to launch a TensorBoard displaying all runs within this collection
            % and saves it to a file inside the RunCollection's path
            %
            % Returns:
            %   file : string
            %     Path to the shell script, will match `.fileShellScriptTensorboard`

            f = rc.fileShellScriptTensorboard;
            fid = fopen(f, 'w');
            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '%s "$@"\n', rc.getTensorboardCommand(varargin{:}));
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end

        function text = generateSummaryText(rc)
            newline = sprintf('\n'); %#ok<SPRINTFN>
            text = sprintf('%s "%s" (%d runs total)\n', ...
                 class(rc), rc.name, rc.nRunsTotal);
            text = cat(2, text, sprintf('  Path: %s\n', rc.path));
            text = cat(2, text, sprintf('  Dataset Collection "%s" (%d datasets) in %s\n\n', ...
                 rc.datasetCollection.name, rc.nDatasets, rc.datasetCollection.path));

            sep = sprintf('------------------------\n');
            text = cat(2, text, sprintf('  %s\n  %d Run Specifications:\n\n', sep, rc.nRunSpecs));
            for s = 1:rc.nRunSpecs
                text = cat(2, text, rc.runSpecs(s).generateSummaryText(4, s), newline);
            end

            text = cat(2, text, sprintf('  %s\n  %d Parameter Settings:\n\n', sep, rc.nParams));
            for p = 1:rc.nParams
               text = cat(2, text, rc.params(p).generateSummaryText(4, p), newline);
            end
        end

        function f = writeSummaryText(rc)
            % Generates a text file sumamrizing the details of the RunSpecs and RunParams
            % within this run collection.
            %
            % Returns:
            %   file : string
            %     Path to the text file, which will match .fileSummaryText
            f = rc.fileSummaryText;
            fid = fopen(f, 'w');
            fprintf(fid, rc.generateSummaryText());
            fclose(fid);
            LFADS.Utils.chmod('uga+rx', f);
        end

        function loadSequenceData(rc, reload)
            % Call `loadSequenceData` on each run in this collection

            if nargin < 2
                reload = false;
            end
            prog = LFADS.Utils.ProgressBar(rc.nRunsTotal, 'Loading sequence data');
            for i = 1:rc.nRunsTotal
                prog.update(i);
                rc.runs(i).loadSequenceData(reload);
            end
            prog.finish();
        end
        
        function tf = checkPosteriorMeansExist(rc, verbose)
            if nargin < 2
                verbose = false;
            end
            tf = false(rc.nRunSpecs, rc.nParams);
            
            for iR = 1:rc.nRunsTotal
                tf(iR) = rc.runs(iR).checkPosteriorMeansExist(verbose);
            end
        end

        function loadPosteriorMeans(rc, reload)
            % Call `loadPosteriorMeans` on each run in this collection

            if nargin < 2
                reload = false;
            end
            prog = LFADS.Utils.ProgressBar(rc.nRunsTotal, 'Loading posterior mean samples');
            for i = 1:rc.nRunsTotal
                prog.update(i);
                rc.runs(i).loadPosteriorMeans(reload);
            end
            prog.finish();
        end
    end

    methods (Access = protected)
       function header = getHeader(rc)
          if ~isscalar(rc)
             header = getHeader@matlab.mixin.CustomDisplay(rc);
          else
             className = class(rc);
             header = sprintf('%s "%s" (%d runs total)\n  Dataset Collection "%s" (%d datasets) in %s\n', ...
                 className, rc.name, rc.nRunsTotal, rc.datasetCollection.name, rc.nDatasets, rc.datasetCollection.path);
             header = cat(2, header, sprintf('  Path: %s\n\n', rc.path));
             header = cat(2, header, sprintf('  %d parameter settings\n', rc.nParams));
             for p = 1:rc.nParams
                 header = cat(2, header, sprintf('  [%d %s %s] %s "%s" %s\n', ...
                     p, rc.params(p).generateHashName, rc.params(p).generateInputDataHashName, ...
                     class(rc.params(p)), rc.params(p).name, rc.params(p).generateShortDifferencesString()));
             end

             header = cat(2, header, sprintf('\n  %d run specifications\n', rc.nRunSpecs));
             for s = 1:rc.nRunSpecs
                 header = cat(2, header, sprintf('  [%2d] %s\n', s, rc.runSpecs(s).getFirstLineHeader()));
             end
          end
       end
    end
end
