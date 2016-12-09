classdef RunCollection < LFADS.RunCollection
    methods
        function rc = RunCollection(rootPath, name, datasetCollection, runParams)
            rc = rc@LFADS.RunCollection(rootPath, name, datasetCollection, runParams);
        end
        
        function addRun(rc, r)
            assert(isa(r, 'PierreEricLFADS.Run'), 'Must be PierreEricLFADS.Run instance');
            addRun@LFADS.RunCollection(rc, r);
        end
    end
end