classdef RunCollection < LFADS.RunCollection
    % no need to modify anything here, but feel free to add useful methods
    % and properties as useful
    
    methods
        function rc = RunCollection(rootPath, name, datasetCollection, varargin)
            rc@LFADS.RunCollection(rootPath, name, datasetCollection, varargin{:});
        end
    end
end