
classdef RunSpec < handle & matlab.mixin.CustomDisplay
    % Represents a set of specifications for an LFADS experiment on a specific subset of datasets
    % within a :ref:`LFADS_DatasetCollection`. RunSpec instances are used to create
    % instances, and all runs in a collection share the same parameter settings, which are represented by a shared
    % :ref:`LFADS_RunParams` instance.

    properties
        runClassName = '' % Class name that will be used to create LFADS.Run instances. Or function_handle to class constructor
        
        name char = '' % Name of this run unique within its RunCollection, will be used as subfolder on disk

        comment char = '' % Textual comment for convenience

        version uint32 = 2; % Internal versioning allowing for graceful evolution of path settings
    end

    properties
        datasets % Array of :ref:`LFADS_Dataset` instances which this particular Run will utilize
        datasetIndsInCollection % indices of each dataset into datasetCollection.datasets
        
        datasetCollection % Dataset collection used by this run (and all runs in the same RunCollection)
    end

    properties(Dependent)
        nDatasets % Number of datasets used by this run
    end

    methods
        function r = RunSpec(varargin)
            % run = RunSpec(name, runClassName, datasetCollection[, datasetIndicesOrNames])
            %
            % Parameters
            % ------------
            % name : string
            %   Unique name for this run within the collection
            %
            % runClassName : string
            %   Class name that will be used to create LFADS.Run instances using
            % 
            % datasetCollection : :ref:`LFADS_DatasetCollection` instance
            %   DatasetCollection in which this run will find its datasets
            %
            % datasetIndicesOrNames : vector of integers or cellstr
            %   if numeric vector, will select datasets by index in datasetCollection. If cellstr of names, will
            %   search for datasets by name within datasetCollection.

            p = inputParser();
            p.addOptional('name', '', @ischar);
            p.addOptional('runClassName', '', @ischar);
            p.addOptional('datasetCollection', [], @(x) isa(x, 'LFADS.DatasetCollection'));
            p.addOptional('datasetIndicesOrNames', [], @(x) true);
            p.parse(varargin{:});
            
            r.name = p.Results.name;
            r.runClassName = p.Results.runClassName;
            r.datasetCollection = p.Results.datasetCollection;

            if ~isempty(p.Results.datasetIndicesOrNames)
                r.selectDatasets(p.Results.datasetIndicesOrNames);
            end
        end
        
        function tf = eq(a, b)
            % Overloaded == operator to enable equality if name,
            % datasetCollection, and datasets all match
            
            tf = LFADS.Utils.bsxfun_object(@compare, a, b);
                
            function tf = compare(a,b)
                tf = strcmp(a.name, b.name) ...
                    && isequal(a.datasets, b.datasets) && isequal(a.datasetCollection, b.datasetCollection);
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
            if isempty(datasetIndicesOrNames)
                r.selectDatasets([]);
            elseif isnumeric(datasetIndicesOrNames) || islogical(datasetIndicesOrNames)
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
            r.datasetIndsInCollection = LFADS.Utils.vectorMaskToIndices(idx);
        end

        function selectDatasetsByName(r, names)
            % Specify the datasets that this Run should use from the set of datasets in its DatasetCollection using name matching
            %
            % Parameters
            % ------------
            % names : string or cellstr
            %   Name or names to search for within this run's DatasetCollection's array of datasets

            [r.datasets, r.datasetIndsInCollection] = r.datasetCollection.matchDatasetsByName(names);
        end
        
        function text = generateSummaryText(r, indent, index)
            % Generate a string representation of this Run Spec's info
            %
            % Returns:
            %   text : string
            if nargin < 2
                indent = 0;
            end
            if nargin > 2
                indexStr = sprintf('[runSpec %d] ', index);
            else
                indexStr = '';
            end
            indent = blanks(indent);
            text = sprintf('%s%s%s\n', indent, indexStr, r.getFirstLineHeader());
            for s = 1:r.nDatasets
                text = cat(2, text, sprintf('%s  [ds %d] %s\n', indent, r.datasetIndsInCollection(s), r.datasets(s).getFirstLineHeader()));
            end
        end
    end

    methods(Hidden)
        function h = getFirstLineHeader(r)
            className = class(r);
            h = sprintf('%s "%s" (%d datasets)', className, r.name, r.nDatasets);
        end
    end

    methods (Access = protected)
       function header = getHeader(r)
          if ~isscalar(r)
             header = getHeader@matlab.mixin.CustomDisplay(r);
          else
             header = sprintf('%s\n  %d datasets from collection %s\n', r.getFirstLineHeader(), r.nDatasets, r.datasetCollection.name);
             for s = 1:r.nDatasets
                 header = cat(2, header, sprintf('    [ds %d] %s\n', r.datasetIndsInCollection(s), r.datasets(s).getFirstLineHeader()));
             end
          end
       end
    end
end
