classdef RunSpec < handle & matlab.mixin.CustomDisplay
    % Represents a set of specifications for an LFADS experiment on a specific subset of datasets
    % within a :ref:`LFADS_DatasetCollection`. RunSpec instances are used to create
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.

    methods(Abstract)
        cls = getRunClassName(r);
        % Return the name of the subclass of :ref:`LFADS_Run` to be used
        % when generating runs
        % 
        % Returns:
        %   cls : string 
        %     Class name which subclasses from :ref:`LFADS_Run`
        
    end
    
    properties
        name char = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk

        comment char = '' % Textual comment for convenience

        version uint32 = 2; % Internal versioning allowing for graceful evolution of path settings
    end

    properties
        datasets % Array of :ref:`LFADS_Dataset` instances which this particular Run will utilize

        datasetCollection % Dataset collection used by this run (and all runs in the same RunCollection)
    end

    properties(Dependent)
        nDatasets % Number of datasets used by this run
    end

    methods
        function r = RunSpec(name, datasetCollection, datasetIndicesOrNames)
            % run = RunSpec(name, datasetCollection[, datasetIndicesOrNames])
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

            r.name = name;
            r.datasetCollection = datasetCollection;

            if nargin > 2
                r.selectDatasets(datasetIndicesOrNames);
            end
        end

        function n = get.nDatasets(r)
            n = numel(r.datasets);
        end

        function selectDatasets(r, datasetIndicesOrNames)
            % selectDatasets(datasetIndicesOrNames)
            %
            % Parameters
            % -----------
            % datasetIndicesOrNames : vector of integers or cellstr
            %   if numeric vector, will select datasets by index in datasetCollection. If cellstr of names, will
            %   search for datasets by name within datasetCollection.
            %
            if isnumeric(datasetIndicesOrNames) || islogical(datasetIndicesOrNames)
                r.selectDatasetsByIndex(datasetIndicesOrNames);
            else
                r.selectDatasetsByName(datasetIndicesOrNames);
            end
        end

        function selectDatasetsByIndex(r, idx)
            % Specify the datasets that this Run should use from the set of datasets in its DatasetCollection by indices
            %
            % Parameters
            % ------------
            % idx : logical mask or indices
            %   Selection applied to this run's DatasetCollection's array of datasets

            r.datasets = r.datasetCollection.datasets(idx);
        end

        function selectDatasetsByName(r, names)
            % Specify the datasets that this Run should use from the set of datasets in its DatasetCollection using name matching
            %
            % Parameters
            % ------------
            % names : string or cellstr
            %   Name or names to search for within this run's DatasetCollection's array of datasets

            r.datasets = r.datasetCollection.matchDatasetsByName(names);
        end
    end

    methods(Hidden)
        function h = getFirstLineHeader(r)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(r);
            h = sprintf('%s %s (%d datasets)\n', className, r.name, r.nDatasets);
        end
    end

    methods (Access = protected)
       function header = getHeader(r)
          if ~isscalar(r)
             header = getHeader@matlab.mixin.CustomDisplay(r);
          else
             rc = r.runCollection;
             header = sprintf('%s\n  %d datasets in %s\n', r.getFirstLineHeader(), r.nDatasets, rc.path);
             for s = 1:r.nDatasets
                 header = cat(2, header, sprintf('    [%2d] %s', s, r.datasets(s).getHeader()));
             end
          end
       end
    end
end
