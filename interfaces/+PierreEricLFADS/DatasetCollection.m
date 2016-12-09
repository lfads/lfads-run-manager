classdef DatasetCollection < LFADS.DatasetCollection
    methods
        function ds = DatasetCollection(path)
            ds = ds@LFADS.DatasetCollection(path);
        end
        
        function autoDetectDatasets(dc)
            dc.clearDatasets;
            
            % automatically find all files within and build datasets
            files = dir(dc.path);
            for iF = 1:numel(files)
                if strncmp(files(iF).name, '.', 1), continue, end
                info = files(iF);
                [~, ~, ext] = fileparts(info.name);
                if ~strcmp(ext, '.mat'), continue; end
                ds = PierreEricLFADS.Dataset(dc, info.name);
                dc.addDataset(ds);
            end
        end
        
        function filterHasHighSNRChannels(dc)
            dc.loadInfo();
            mask = arrayfun(@(ds) ds.nChannelsHighSNR > 0, dc.datasets);
            dc.filterDatasets(mask);
        end
        
        function filterBestSaveTagEachDate(dc)
            % for days with multiple saveTags - to avoid overlap, we only take
            %    the savetag with the most trials
            allDates = arrayfun(@(ds) ds.datenum, dc.datasets);
            [uniqueDays, ~, udAssignments] = unique(allDates);
            
            maskKeep = falsevec(dc.nDatasets);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets);
            for nd = 1:numel(uniqueDays)
                thisDayDatasetInds = find(udAssignments == nd);
                
                numTrials = nTrials(thisDayDatasetInds);
                [~, whichSetToKeep] = max(numTrials);
                maskKeep(thisDayDatasetInds(whichSetToKeep)) = true;
            end
            
            dc.filterDatasets(maskKeep);
        end
        
        function addDataset(dc, ds)
            assert(isa(ds, 'PierreEricLFADS.Dataset'), 'Must be PierreEricLFADS.Dataset instance');
            addDataset@LFADS.DatasetCollection(dc, ds);
        end
        
        function t = getDatasetInfoTable(dc)
            dc.loadInfo();
            rowNames = arrayfun(@(ds) ds.name, dc.datasets, 'UniformOutput', false);
            date = arrayfun(@(ds) ds.datestr, dc.datasets, 'UniformOutput', false);
            saveTags = arrayfun(@(ds) strjoin(ds.saveTags, ','), dc.datasets, 'UniformOutput', false);
            nChannels = arrayfun(@(ds) ds.nChannels, dc.datasets, 'UniformOutput', true);
            nChannelsHighSNR = arrayfun(@(ds) ds.nChannelsHighSNR, dc.datasets, 'UniformOutput', true);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets, 'UniformOutput', true);
            
            t = table(date, saveTags, nTrials, nChannels, nChannelsHighSNR, 'RowNames', rowNames);
        end
    end
end
