classdef Dataset < LFADS.Dataset
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
        end

        function data = loadData(ds)
            % load this dataset's data file from .path
            in = load(ds.path);
            data = in.data;
        end

        function loadInfo(ds)
            % Load this Dataset's metadata if not already loaded

            if ds.infoLoaded, return; end

            % modify this to extract the metadata loaded from the data file
            % data = ds.loadData();
            % ds.subject = data.subject;
            % ds.saveTags = data.saveTags;
            % ds.datenum  = data.datenum;
            % ds.nChannels = data.nChannels;
            % ds.nTrials = numel(data.trials);

            ds.infoLoaded = true;
        end

    end
end
