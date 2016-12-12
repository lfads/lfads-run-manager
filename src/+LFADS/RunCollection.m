classdef RunCollection < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
% Wraps a set of runs often sharing common parameter settings but utilizing
% different sets of datasets

    properties
        name = ''
        comment = ''
        rootPath = '';
        
        fileShellScriptTensorboard
        
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
    end
    
    methods
        function rc = RunCollection(rootPath, name, datasetCollection, runParams)
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
        
        function clearRuns(rc)
            rc.runs = []; 
        end 
        
        function filterRuns(rc, mask)
            rc.runs = rc.runs(mask);
        end
        
        function n = get.nRuns(rc)
            n = numel(rc.runs);
        end
        
        function n = get.nDatasets(rc)
            n = rc.datasetCollection.nDatasets;
        end
        
        function str = getTensorboardCommand(rc)
            runEntry = cellvec(rc.nRuns);
            for r = 1:rc.nRuns
                runEntry{r} = sprintf('%s:%s', rc.runs(r).name, rc.runs(r).pathLFADSOutput);
            end
            str = sprintf('tensorboard --logdir=%s', strjoin(runEntry, ','));
        end 
        
        function f = writeTensorboardShellScript(rc)
            f = rc.fileShellScriptTensorboard;
            fid = fopen(f, 'w');
            fprintf(fid, rc.getTensorboardCommand());
            fclose(fid);
            chmod('uga+rx', f);
        end

        function addRun(rc, r)
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
    end
    
    methods
        function rc2 = copyClearRuns(rc)
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