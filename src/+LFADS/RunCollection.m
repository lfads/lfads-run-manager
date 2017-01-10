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

        version uint32 = 3; % version used for backwards compatibility
        
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
        path % unique folder given by rootPath/name
        fileShellScriptTensorboard % path the location where a shell script to launch TensorBoard for all runs will be written
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
            rc.runParams = p.Results.runParams;
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
            if isempty(rc.runSpecs)
                rc.runSpecs = runSpecs;
            else
                % check for existing runs by name and replace them
                names = arrayfun(@(oldR) oldR.name, rc.runSpecs, 'UniformOutput', false);
                [tf, idx] = ismember(r.name, names);
                if any(tf)
                    warning('Replacing existing runs with matching name');
                    rc.runs(idx(tf)) = runSpecs(tf);
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

        function filterRuns(rc, mask)
            % Apply selection mask to list of runs
            %
            % Parameters:
            %   mask : logical or indices
            %     selection mask applied to `.runs`
            
            rc.runs = rc.runs(mask);
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
            if isempty(rc.runParams)
                rc.runParams = params;
            else
                rc.runParams = cat(1, rc.runSpecs, runParams);
            end
            
            if rc.version < 3 && numel(rc.runParams) > 1 
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
            
            rc.runs = [];
            for iS = rc.nRunSpecs:-1:1
                spec = rc.runSpecs(iS);
                clsFn = str2func(spec.getRunClassName());
                
                for iP = rc.nParams:-1:1
                    % create the run
                    rc.runs(iS, iP) = clsFn(spec.name, rc, rc.params(iP), spec.datasets);
                end
            end
        end
    end
        
    methods
        function p = get.path(rc)
            if rc.version < 3
                p = fullfile(rc.rootPath, rc.name);
            else
                paramSuffix = rc.params(1).generateSuffix();
                p = fullfile(rc.rootPath, [rc.name '_' paramSuffix]);
            end
        end

        function f = get.fileShellScriptTensorboard(r)
            f = fullfile(r.path, 'launch_tensorboard.sh');
        end
        
        function n = get.nParams(rc)
            n = numel(rc.params);
        end

        function n = get.nRunsTotal(rc)
            n = numel(rc.runs);
        end
        
        function n = get.nRunSpecs(rc)
            n = numel(rc.runSpecs);
        end

        function n = get.nDatasets(rc)
            n = rc.datasetCollection.nDatasets;
        end
        
        function [runSpecs, idx] = findRunSpecs(rc, runSpecSearch)
            % Find run specs that match runSpecSearch. Throws an error if any run is not
            % found. If runSpecSearch is string or cellstr, matches by
            % name. If runSpecSearch is an array of RunSpec instances,
            % finds them using isequal. If runSpecSearch is a vector of
            % indices, selects from runSpecs directly
            %
            % Args:
            %   runSpecSearch : array of LFADS.RunSpec, strings, or indices into .runSpecs
            %
            % Returns:
            %   runs (LFADS.Run) : matching runs
            %   idx : vector of indices into `.runSpecs` of matching run
            %     instances
            %
            
            if ischar(runSpecSearch)
                runSpecSearch = {runSpecSearch};
            end
            if iscellstr(runSpecSearch)
                [tf, idx] = ismember(runSpecSearch, {rc.runSpecs.name});
                assert(all(tf), 'Some run spec names could not be found in this RunCollection');
            elseif isa(runSpecSearch, 'LFADS.RunSpec')
                [tf, idx] = ismember(runSpecSearch, rc.runSpecs);
                assert(all(tf), 'Some run spec instances could not be found in this RunCollection');
            else
                % assume is selection
                idx = runSpecSearch;
            end
            idx = LFADS.Utils.makecol(idx(:));
            runSpecs = rc.runsSpecs(idx);
        end
        
        function [params, idx] = findRunParams(rc, paramSearch)
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
            %   idx : vector of indices into `.params` of matches
            %
            
            if isa(paramSearch, 'LFADS.RunParams')
                [tf, idx] = ismember(paramSearch, rc.params);
                assert(all(tf), 'Some run params instances could not be found in this RunCollection');
            else
                % assume is selection
                idx = paramSearch;
            end
            idx = LFADS.Utils.makecol(idx(:));
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
            
            [~, rowIdx] = rc.findRunSpecs(runSpecSearch);
            [~, colIdx] = rc.findParams(paramSearch);
            runs = rc.runs(rowIdx, colIdx);
        end

        function str = getTensorboardCommand(rc)
            % Generates the shell command text to launch a TensorBoard displaying all runs within this collection
            %
            % Returns:
            %   cmd : string
            %     Command which can luanch TensorBoard from command line

            runEntry = cellvec(rc.nRuns);
            for r = 1:rc.nRuns
                runEntry{r} = sprintf('%s:%s', rc.runs(r).name, rc.runs(r).pathLFADSOutput);
            end
            str = sprintf('tensorboard --logdir=%s', strjoin(runEntry, ','));
        end

        function f = writeTensorboardShellScript(rc)
            % Generates the shell command text to launch a TensorBoard displaying all runs within this collection
            % and saves it to a file inside the RunCollection's path
            %
            % Returns:
            %   file : string
            %     Path to the shell script, will match `.fileShellScriptTensorboard`

            f = rc.fileShellScriptTensorboard;
            fid = fopen(f, 'w');
            fprintf(fid, rc.getTensorboardCommand());
            fclose(fid);
            chmod('uga+rx', f);
        end
        
        function loadSequenceData(rc, reload)
            % Call `loadSequenceData` on each run in this collection

            if nargin < 2
                reload = false;
            end
            prog = LFADS.Utils.ProgressBar(rc.nRuns, 'Loading sequence data');
            for i = 1:rc.nRuns
                prog.update(i);
                rc.runs(i).loadSequenceData(reload);
            end
            prog.finish();
        end
        
        function loadPosteriorMeans(rc, reload)
            % Call `loadPosteriorMeans` on each run in this collection

            if nargin < 2
                reload = false;
            end
            prog = LFADS.Utils.ProgressBar(rc.nRuns, 'Loading posterior mean samples');
            for i = 1:rc.nRuns
                prog.update(i);
                rc.runs(i).loadPosteriorMeans(reload);
            end
            prog.finish();
        end
    end

%     methods
%         function rc2 = copyClearAll(rc)
%             % Make a copy of this run collection, but without the run specs and params inside
%             %
%             % Returns
%             % ---------
%             % rc : :ref:`LFADS_RunCollection`
%             %   Copy of RunCollection sans runs
%             rc2 = copy(rc);
%             rc2.clearAll();
%         end
%     end

    methods (Access = protected)
       function header = getHeader(rc)
          if ~isscalar(rc)
             header = getHeader@matlab.mixin.CustomDisplay(rc);
          else
             className = matlab.mixin.CustomDisplay.getClassNameForHeader(rc);
             header = sprintf('%s %s in %s\n', className, rc.name, rc.path);
             header = cat(2, header, sprintf('  %d RunParams settings\n', rc.nParams));
             for p = 1:rc.nParams
                 header = cat(2, header, sprintf('  [%2d] %s\n', rc.params(p).getSuffix));
             end
             
             header = cat(2, header, sprintf('  %d runs per param setting, %d runs total\n', rc.nRunsEachParam, rc.nRunsTotal));
             for s = 1:rc.nRunsEachParam
                 header = cat(2, header, sprintf('  [%2d] %s', s, rc.runs(s).getFirstLineHeader()));
             end
          end
       end
    end
end
