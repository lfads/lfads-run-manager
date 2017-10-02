classdef RunSpec < LFADS.RunSpec

    methods
        function r = RunSpec(name, datasetCollection, datasetIndicesOrNames)
            % run = RunSpec(name, datasetCollection, datasetIndicesOrNames)
            %
            % Parameters
            % ------------
            % name : string
            %   Unique name for this run within the collection
            %
            % datasetCollection : :ref:`LFADS_DatasetCollection` instance
            %   DatasetCollection in which this run will find its datasets
            %
            % datasetIndicesOrNames : vector of integers or cellstr
            %   if numeric vector, will select datasets by index in datasetCollection. If cellstr of names, will
            %   search for datasets by name within datasetCollection.

            runClassName = sprintf('%s.Run', LFADS.Utils.getPackage());
            r = LFADS.RunSpec(name, runClassName, datasetCollection, datasetIndicesOrNames);
        end
    end
end
