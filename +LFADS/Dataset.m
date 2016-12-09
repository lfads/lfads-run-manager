classdef Dataset < handle & matlab.mixin.CustomDisplay
% A single-day of raw data processed in a common way

    properties
        name = ''
        comment = ''
        relPath = ''
        collection
    end
    
    properties
        infoLoaded = false;
        subject = ''
        saveTags
        datenum
        nChannels
        nTrials
    end
    
    properties(Dependent)
        path
        datestr
    end
    
    methods(Abstract)
        loadInfoFromData(data);
    end
    
    methods
        function ds = Dataset(collection, relPath)
            [~, ds.name] = fileparts(relPath);
            ds.relPath = relPath;
            ds.collection = collection;
        end
        
        function p = get.path(ds)
            if isempty(ds.collection)
                p = ds.relPath;
            else
                p = fullfile(ds.collection.path, ds.relPath);
            end
        end
        
        function data = loadData(ds)
            in = load(ds.path);
            data = in.data;
        end
        
        function reloadInfo(ds)
            ds.infoLoaded = false;
            ds.loadInfo();
        end
        
        function loadInfo(ds)
            if ds.infoLoaded, return; end
            
            data = ds.loadData();
            
            ds.loadInfoFromData(data);
            
            ds.infoLoaded = true;
        end
        
        function ds = get.datestr(ds)
            if isempty(ds.datenum)
                ds = '';
            else
                ds = datestr(ds.datenum, 'yyyy-mm-dd'); %#ok<CPROP>
            end
        end
    end
    
    methods (Access = protected)
       function header = getHeader(ds)
          if ~isscalar(ds)
             header = getHeader@matlab.mixin.CustomDisplay(ds);
          else
             className = matlab.mixin.CustomDisplay.getClassNameForHeader(ds);
             newHeader = sprintf('%s %s in %s', className, ds.name, ds.relPath);
             header = sprintf('%s\n',newHeader);
          end
       end
    end
    
end
