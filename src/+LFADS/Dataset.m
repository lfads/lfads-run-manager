classdef Dataset < handle & matlab.mixin.CustomDisplay
    % A single-day collection of raw data to be processed by LFADS

    methods(Abstract)
        % These methods you must define in order to retrieve whatever metadata are needed from this dataset

       loadInfoFromData(data);
       % This method should load any metadata about the dataset. There are no specific requirements. This simply
       % enables an easy way to ensure that this potentially time-consuming computation is performed only once.
    end

    properties
        % Information about the dataset's name and location

        name = ''
        % Unique, semi-permanent identifier for the dataset. This will be used for paths and filenames and
        % variables within the LFADS model, so it's best to pick a scheme and stick to it.

        comment = ''
        % Textual comment for convenience

        relPath = ''
        % Path to the Matlab data file that will be loaded when loadData is called, relative to the path of the
        % :ref:`LFADS_DatasetCollection` to which this Dataset belongs. If you override loadData, this can be left blank.

        collection
        % :ref:`LFADS_DatasetCollection` to which this Dataset belongs
    end

    properties
        % Useful metadata about this dataset to be loaded by loadInfoFromData. These are not used explicitly
        % and can be left blank.

        infoLoaded = false;
        % Has loadInfo already been called and the info fields populated?

        subject = ''
        % Dataset subject or participant name

        saveTags
        % Data grouping identifiers

        datenum
        % Matlab datenum identifying the collection time of this dataset

        nChannels
        % Number of spike channels recorded in this dataset

        nTrials
        % Number of behavioral trials recorded in this dataset

    end

    properties(Dependent)
        path
        % Full path to data which will be loaded by load data, a concatenation of the DatasetCollection path and relPath

        datestr
        % a string version of datenum
    end

    methods
        function ds = Dataset(collection, relPath)
            % ds = Dataset(collection, relPath)
            % Parameters
            % ------------
            % collection : :ref:`LFADS_DatasetCollection`
            %   DatasetCollection to which this Dataset belongs
            %
            % relPath : string
            %   Relative path to data from `collection.path`

            [~, ds.name] = fileparts(relPath);
            ds.relPath = relPath;
            collection.addDataset(ds);
        end

        function p = get.path(ds)
            if isempty(ds.collection)
                p = ds.relPath;
            else
                p = fullfile(ds.collection.path, ds.relPath);
            end
        end

        function data = loadData(ds)
            % Load this Dataset's data file from .path
            in = load(ds.path);
            data = in.data;
        end

        function loadInfo(ds)
            % Load this Dataset's metadata if not already loaded

            if ds.infoLoaded, return; end
            data = ds.loadData();
            ds.loadInfoFromData(data);
            ds.infoLoaded = true;
        end

        function reloadInfo(ds)
            % Load or reload this Dataset's metadata even if already loaded
            ds.infoLoaded = false;
            ds.loadInfo();
        end

        function ds = get.datestr(ds)
            if isempty(ds.datenum)
                ds = '';
            else
                ds = datestr(ds.datenum, 'yyyy-mm-dd'); %#ok<CPROP>
            end
        end
    end

    methods (Access = protected)
       function header = getHeader(ds)
          if ~isscalar(ds)
             header = getHeader@matlab.mixin.CustomDisplay(ds);
          else
             className = matlab.mixin.CustomDisplay.getClassNameForHeader(ds);
             newHeader = sprintf('%s %s in %s', className, ds.name, ds.relPath);
             header = sprintf('%s\n',newHeader);
          end
       end
    end

end
