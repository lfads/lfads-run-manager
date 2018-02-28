classdef RunSpec < LFADS.RunSpec

    methods
        function r = RunSpec(name, datasetCollection, datasetIndicesOrNames)
            % run = RunSpec(name, datasetCollection, datasetIndicesOrNames)
            % Passes along to LFADS.RunSpec constructor but inserts
            % automatically determined Run class name, e.g.
            % 'StarterExample.Run'
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

            r = r@LFADS.RunSpec(name, '', datasetCollection, datasetIndicesOrNames);
            r.runClassName = strrep(class(r), 'RunSpec', 'Run');
        end
    end
end
