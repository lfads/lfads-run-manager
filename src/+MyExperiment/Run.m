classdef Run < LFADS.Run
    methods
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end

        function [counts, timeVecMs] = generateRatesForDataset(r, dataset, mode, varargin) %#ok<INUSD>
            % Generate binned spike count tensor for a single dataset.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices. You can
            %   include a different subset of the data (or different time
            %   windows) for the alignment process separately from the actual
            %   data exported to LFADS, or return the same for both. Alignment
            %   is only relevant for multi-dataset models. If you wish to use
            %   separate data for alignment, override the method usesDifferentSequenceDataForAlignment
            %   to return true as well. 
            %
            % Returns
            % ----------
            % counts : nTrials x nChannels x nTime tensor
            %   spike counts in time bins in trials x channels x time. These
            %   should be total counts, not normalized rates, as they will be
            %   added during rebinning.
            % timeVecMs: nTime x 1 vector 
            %   of timepoints in milliseconds associated with each time bin. You can start this
            %   wherever you like, but timeVecMs(2) - timeVecMs(1) will be
            %   treated as the spike bin width used when the data are later
            %   rebinned to match run.params.spikeBinMs
            
            data = dataset.loadData();
            
            counts = data.spikes;
            timeVecMs = data.timeMs;
        end
    end
end
