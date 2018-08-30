classdef Dataset < LFADS.Dataset
    properties
        % subpop ids are conceptually similar to region ids but are contiguous
        % whereas regionIds were specified as is in the raw data
        subpopulationIds % N x 1 numeric array for each neuron
    end
    
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
            % you might also wish to set ds.name here,
            % possibly by adding a third argument to the constructor
            % and assigning it to ds.name
        end
    end
end
