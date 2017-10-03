classdef Dataset < LFADS.Dataset
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
            % you might also wish to set ds.name here,
            % possibly by adding a third argument to the constructor
            % and assigning it to ds.name
        end

        function data = loadData(ds)
            % load this dataset's data file from .path
            data = load(ds.path);
        end

        function loadInfo(ds)
            % Load this Dataset's metadata if not already loaded

            if ds.infoLoaded, return; end

            % modify this to extract the metadata loaded from the data file
            data = ds.loadData();
            ds.subject = data.subject;
            ds.saveTags = 1;
            ds.datenum  = datenum(data.datetime);
            ds.nChannels = size(data.spikes, 2);
            ds.nTrials = size(data.spikes, 1);

            ds.infoLoaded = true;
        end

    end
end
