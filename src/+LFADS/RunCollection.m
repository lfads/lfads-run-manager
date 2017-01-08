classdef RunCollection < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
    % A set of :ref:`LFADS_Run` instances sharing common parameter settings but utilizing
    % different subsets of :ref:`LFADS_Dataset` insstances in a :ref:`LFADS_DatasetCollection`.
    % A RunCollection is a logical grouping of otherwise independent LFADS runs, which are constrained
    % to share a common Dataset collection and run parameters.

    properties
        name = '' % Name of this RunCollection, determines its relative path on disk
        comment = '' % Textual comment for convenience
        rootPath = ''; % Root path on disk under which individual Runs will be stored

        datasetCollection % DatasetCollection instance

        params % RunParams instance
    end

    properties(SetAccess=protected)
        runs
    end

    properties(Dependent)
        nRuns
        nDatasets
        path % unique folder within rootPath including name_paramSuffix
        fileShellScriptTensorboard % path the location where a shell script to launch TensorBoard for all runs will be written
    end

    methods
        function rc = RunCollection(rootPath, name, datasetCollection, runParams)
            %  rc = RunCollection(rootPath, name, datasetCollection, runParams)
            %
            % Parameters
            % ------------
            % rootPath : string
            %   location on disk under which run collections will be saved
            % name : string
            %   unique identifier for this RunCollection, also used as subfolder within `rootPath`
            % datasetCollection : :ref:`LFADS_DatasetCollection`
            %   DatasetCollection to be used by runs within this RunCollection
            % runParams : :ref:`LFADS_RunParams`
            %   Parameter settings which will be used uniformly by all runs, will also be serialized in the folder path

            rc.name = name;
            rc.rootPath = rootPath;
            rc.datasetCollection = datasetCollection;
            rc.params = runParams;
        end

        function p = get.path(rc)
            if isempty(rc.params)
                paramStr = '';
            else
                paramStr = ['_', rc.params.generateSuffix()];
            end
            p = fullfile(rc.rootPath, [rc.name, paramStr]);
        end

        function f = get.fileShellScriptTensorboard(r)
            f = fullfile(r.path, 'launch_tensorboard.sh');
        end
        
        function runs = getRunsByName(rc, names)
            % Find runs by their names. Throws an error if any run is not
            % found.
            %
            % Args:
            %   names (string or cellstr) : single name or cell array of
            %     names to find
            % Returns:
            %   runs (LFADS.Run) : matching runs
            %   idx (uint) : list of indices into `.runs` of matching run
            %     instances
            %
            
            if ischar(names)
                names = {names};
            end
            [tf, idx] = ismember(names, {rc.runs.name});
            assert(all(tf), 'Some run names could not be found in this RunCollection');
            
            runs = rc.runs(idx);
        end

        function clearRuns(rc)
            % Flush list of runs

            rc.runs = [];
        end

        function filterRuns(rc, mask)
            % Apply selection mask to list of runs
            %
            % Parameters
            % ------------
            % mask : logical or indices
            %   selection mask applied to `.runs`
            rc.runs = rc.runs(mask);
        end

        function n = get.nRuns(rc)
            n = numel(rc.runs);
        end

        function n = get.nDatasets(rc)
            n = rc.datasetCollection.nDatasets;
        end

        function str = getTensorboardCommand(rc)
            % Generates the shell command text to launch a TensorBoard displaying all runs within this collection
            %
            % Returns
            % ---------
            % cmd : string
            %   Command which can luanch TensorBoard from command line

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
            % Returns
            % ---------
            % file : string
            %   Path to the shell script, will match `.fileShellScriptTensorboard`

            f = rc.fileShellScriptTensorboard;
            fid = fopen(f, 'w');
            fprintf(fid, rc.getTensorboardCommand());
            fclose(fid);
            chmod('uga+rx', f);
        end

        function addRun(rc, r)
            % addRun(run)
            % Adds a :ref:`LFADS_Run` to the collection. Note that :ref:`LFADS_Run` instances are generally added to their Run Collection upon construction, so calling this method is likely unnecessary for the end user.
            %
            % Parameters
            % ------------------
            % run : :ref:`LFADS_Run`
            %   Run to add to the collection.

            if isempty(rc.runs)
                rc.runs = r;
            else
                % check for existing run and replace it
                names = arrayfun(@(oldR) oldR.name, rc.runs, 'UniformOutput', false);
                [tf, idx] = ismember(r.name, names);
                if tf
                    debug('Replacing existing run with matching name %s\n', rc.runs(idx).name);
                    rc.runs(idx) = r;
                else
                    rc.runs(end+1, :) = r;
                end
            end
            r.runCollection = rc;
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

    methods
        function rc2 = copyClearRuns(rc)
            % Make a copy of this run collection, but without the runs inside
            %
            % Returns
            % ---------
            % rc : :ref:`LFADS_RunCollection`
            %   Copy of RunCollection sans runs
            rc2 = copy(rc);
            rc2.clearRuns();
        end
    end

    methods (Access = protected)
       function header = getHeader(rc)
          if ~isscalar(rc)
             header = getHeader@matlab.mixin.CustomDisplay(rc);
          else
             className = matlab.mixin.CustomDisplay.getClassNameForHeader(rc);
             newHeader = sprintf('%s %s', className, rc.name);
             header = sprintf('%s\n  param suffix: %s\n  %d runs in %s\n', ...
                 newHeader, rc.params.generateSuffix(), rc.nRuns, rc.path);

             for s = 1:rc.nRuns
                 header = cat(2, header, sprintf('  [%2d] %s', s, rc.runs(s).getFirstLineHeader()));
             end
          end
       end
    end
end
