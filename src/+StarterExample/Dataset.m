classdef Dataset < LFADS.Dataset
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
        end
        
        function loadInfoFromData(ds, data)
            % modify this to extract the metadata loaded from the data file
            ds.subject = data.subject;
            ds.saveTags = data.saveTags;
            ds.datenum  = data.datenum;
            ds.nChannels = data.nChannels;
            ds.nTrials = data.nTrials;
        end
        
        function data = loadData(ds)
            % load this dataset's data file from .path
            in = load(ds.path);
            data = in.data;
        end
        
    end
end
