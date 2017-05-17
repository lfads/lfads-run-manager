classdef DatasetCollection < LFADS.DatasetCollection    
    methods
        function ds = DatasetCollection(path)
            ds = ds@LFADS.DatasetCollection(path);
        end
        
        function filterHavingMinimumTrials(dc, minTrials)
            % example of a function that will filter down datasets based on
            % their metadata.
            nTrials = cat(1, dc.datasets.nTrials);
            
            % filterDatasets is provided by DatasetCollection
            dc.filterDatasets(nTrials >= minTrials);
        end
    end
end
