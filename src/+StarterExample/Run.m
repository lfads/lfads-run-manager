classdef Run < LFADS.Run
    methods
        function r = Run(varargin) 
           r@LFADS.Run(varargin{:});
        end
        
        function seq = convertDatasetToSequenceStruct(r, dataset)
            % Overwrite this to convert dataset into the format required by
            % sequence data, with fields `y`, `y_time`, `binWidthMs`, and
            % `conditionId`.
            
            nTrials = numel(dataset);
            
            for iT = 1:nTrials
                seq(iT).y = dataset(iT).y;
                seq(iT).y_time = dataset(iT).y_time;
                seq(iT).binWidthMs = dataset(iT).binWidthMs;
                seq(iT).conditionId = dataset(iT).conditionId;
            end
                
            seq = LFADS.Utils.makecol(seq);
        end
    end
end
