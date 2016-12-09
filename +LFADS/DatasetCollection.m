classdef DatasetCollection < handle & matlab.mixin.CustomDisplay
% A multi-day collection of raw data processed in a common way, e.g. Pierre
% from Eric with rehtresholded broadband

    properties
        name = ''
        comment = ''
        path = ''
    end
    
    properties(SetAccess=protected)
        datasets
    end
    
    properties(Dependent)
        nDatasets
    end
    
    methods
        function ds = DatasetCollection(path)
            ds.path = path;
            [~, ds.name] = fileparts(path);
        end

        function addDataset(dc, ds)
            if isempty(dc.datasets)
                dc.datasets = ds;
            else
                dc.datasets(end+1, :) = ds;
            end
        end
        
        function clearDatasets(dc)
            dc.datasets = []; 
        end 
        
        function n = get.nDatasets(dc)
            n = numel(dc.datasets);
        end
        
        function loadInfo(dc)
            prog = ProgressBar(dc.nDatasets, 'Loading info');
            for i = 1:dc.nDatasets
                prog.update(i, 'Loading info for dataset %s', dc.datasets(i).name);
                dc.datasets(i).loadInfo();
            end
            prog.finish();
        end
        
        function filterDatasets(dc, mask)
            dc.datasets = dc.datasets(mask);
        end
        
        function datasets = matchDatasetsByName(dc, names)
            [tf, which] = ismember(names, {dc.datasets.name});
            assert(all(tf), 'Missing datasets %s', LFADS.Utils.strjoin(names(~tf), ', '));
            datasets = dc.datasets(which);
        end

        function t = getDatasetInfoTable(dc)
            dc.loadInfo();
            rowNames = arrayfun(@(ds) ds.name, dc.datasets, 'UniformOutput', false);
            date = arrayfun(@(ds) ds.datestr, dc.datasets, 'UniformOutput', false);
            saveTags = arrayfun(@(ds) strjoin(ds.saveTags, ','), dc.datasets, 'UniformOutput', false);
            nChannels = arrayfun(@(ds) ds.nChannels, dc.datasets, 'UniformOutput', true);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets, 'UniformOutput', true);
            
            t = table(date, saveTags, nTrials, nChannels, nChannelsHighSNR, 'RowNames', rowNames);
        end
        
    end
    
    methods (Access = protected)
       function header = getHeader(dc)
          if ~isscalar(dc)
             header = getHeader@matlab.mixin.CustomDisplay(dc);
          else
             className = matlab.mixin.CustomDisplay.getClassNameForHeader(dc);
             newHeader = sprintf('%s %s', className, dc.name);
             header = sprintf('%s\n  %d datasets in %s\n',newHeader, dc.nDatasets, dc.path);
             
             for s = 1:dc.nDatasets
                 header = cat(2, header, sprintf('  %2d  %s\n', s, dc.datasets(s).name));
             end
          end
       end
    end
end
