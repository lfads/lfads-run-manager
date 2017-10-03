classdef DatasetCollection < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
    % A collection of multiple Datasets to be processed by LFADS as a cohesive group, either using stitching to
    % incorporate multiple datasets simultaneously, or individually to multiple independent LFADS fits.

    properties
        % Information about this DatasetCollection's name and location

        name char = '' % Name of the dataset collection, will be used to construct folder paths on disk
        comment char = '' % Textual comment for convenience
        path char = '' % path to this dataset collection's root location on disk, if applicable
    end

    properties(SetAccess=protected)
        datasets % Array of :ref:`LFADS_Dataset` instances
    end

    properties(Dependent)
        nDatasets % Number of datasets in this collection
    end

    methods
        function ds = DatasetCollection(path, name)
            % ds = DatasetCollection(path, [name=leaf of path])
            % Parameters
            % ------------------
            % path : string
            %   Path to root of dataset collection. Can be blank if you've overriden `loadData` and aren't loading
            %   datasets from disk directly as `.mat` files
            % name : string
            %   Dataset collection name. Defaults to the leaf folder in `path`

            if nargin > 0
                ds.path = path;
                if nargin < 2
                    [~, ds.name] = fileparts(path);
                else
                    ds.name = name;
                end
            end
        end

        function addDataset(dc, ds)
            % addDataset(dataset)
            % Adds a dataset to the collection. Note that :ref:`LFADS_Dataset` instances are added to their dataset collection upon construction, so calling this method is likely unnecessary for the end user.
            %
            % Parameters
            % ------------------
            % dataset : :ref:`LFADS_Dataset`
            %   Dataset to add to the collection.

            if isempty(dc.datasets)
                dc.datasets = ds;
            else
                % check for existing dataset and replace it
                names = arrayfun(@(old) old.name, dc.datasets, 'UniformOutput', false);
                [tf, idx] = ismember(ds.name, names);
                if tf
                    fprintf('Replacing existing dataset with matching name %s\n', ds.name);
                    dc.datasets(idx) = ds;
                else
                    dc.datasets(end+1, :) = ds;
                end
            end
            ds.collection = dc;
        end

        function clearDatasets(dc)
            % Flush all datasets from this collection
            dc.datasets = [];
        end

        function n = get.nDatasets(dc)
            n = numel(dc.datasets);
        end
        
        function reloadInfo(dc)
            % Call `reloadInfo` on each dataset in this collection

            prog = LFADS.Utils.ProgressBar(dc.nDatasets, 'Loading info');
            for i = 1:dc.nDatasets
                prog.update(i, 'Loading info for dataset %s', dc.datasets(i).name);
                dc.datasets(i).reloadInfo();
            end
            prog.finish();
        end

        function loadInfo(dc)
            % Call `loadInfo` on each dataset in this collection

            prog = LFADS.Utils.ProgressBar(dc.nDatasets, 'Loading info');
            for i = 1:dc.nDatasets
                prog.update(i, 'Loading info for dataset %s', dc.datasets(i).name);
                dc.datasets(i).loadInfo();
            end
            prog.finish();
        end

        function filterDatasets(dc, mask)
            % Retain a selected subsets of datasets, effectively doing datasets = datasets(mask)
            %
            % Parameters
            % ------------------
            % mask : logical or indices
            %   Selection applied to `.datasets`

            dc.datasets = dc.datasets(mask);
        end
        
        function filterHavingMinimumTrials(dc, minTrials)
            nTrials = cat(1, dc.datasets.nTrials);
            dc.filterDatasets(nTrials >= minTrials);
        end
        
        function filterHavingMinimumTrialsForBatchSize(dc, runParams)
            minTrials = ceil(runParams.c_batch_size * (runParams.trainToTestRatio+1));
            dc.filterHavingMinimumTrials(minTrials);
        end

        function [datasets, idx] = matchDatasetsByName(dc, names)
            % Returns the subset of datasets in this collection matching a name in names.
            %
            % Args:
            %   names : string or cellstr
            %     Name or names of datasets to find
            %
            % Returns:
            %   datasets : LFADS.Dataset array
            %   idx : list of indices into datasets array

            [tf, idx] = ismember(names, {dc.datasets.name});
            assert(all(tf), 'Missing datasets %s', LFADS.Utils.strjoin(names(~tf), ', '));
            datasets = dc.datasets(idx);
        end

        function t = getDatasetInfoTable(dc)
            % Build a `table` of datasets and associated metadata within this dataset. This will call `loadInfo` and load info for datasets whose metadata has not yet been loaded.

            dc.loadInfo();
            rowNames = arrayfun(@(ds) ds.name, dc.datasets, 'UniformOutput', false);
            date = arrayfun(@(ds) datetime(ds.datenum, 'ConvertFrom','datenum'), dc.datasets, 'UniformOutput', false);
            saveTags = arrayfun(@(ds) LFADS.Utils.strjoin(ds.saveTags, ','), dc.datasets, 'UniformOutput', false);
            nChannels = arrayfun(@(ds) ds.nChannels, dc.datasets, 'UniformOutput', true);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets, 'UniformOutput', true);

            t = table(date, saveTags, nTrials, nChannels, 'RowNames', rowNames);
        end
    end

    methods (Access = protected)
       function header = getHeader(dc)
          if ~isscalar(dc)
             header = getHeader@matlab.mixin.CustomDisplay(dc);
          else
             className = class(dc);
             newHeader = sprintf('%s "%s"', className, dc.name);
             header = sprintf('%s\n  %d datasets in %s\n',newHeader, dc.nDatasets, dc.path);

             for s = 1:dc.nDatasets
                 header = cat(2, header, sprintf('  [%2d] %s\n', s, dc.datasets(s).getFirstLineHeader()));
             end
          end
       end
    end
end
