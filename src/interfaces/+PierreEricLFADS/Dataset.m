classdef Dataset < LFADS.Dataset
% A single-day of raw data processed in a common way
    properties(SetAccess=protected)
        nChannelsHighSNR
    end
    
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
        end
        
        function loadInfoFromData(ds, data)
            ds.subject = strtok(data.dataset, ' ');
            
            % find save tag
            match = regexp(data.dataset, 'saveTag(Group)? (?<saveTag>[\d,])', 'names');
            ds.saveTags = sscanf(match.saveTag, '%d,')';
            % find datestr
            match = regexp(data.dataset, '(?<datestr>\d{4}-\d{2}-\d{2})', 'names');
            ds.datenum  = datenum(match.datestr);
            
            ds.nChannels = data.nUnits;
            ds.nChannelsHighSNR = numel(data.chSNR >= 2);
            ds.nTrials = data.nTrials;
        end
    end
end
