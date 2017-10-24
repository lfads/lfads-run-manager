classdef Run < LFADS.Run
    methods
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end

        function out = generateCountsForDataset(r, dataset, mode, varargin) %#ok<INUSL,INUSD>
            % Generate binned spike count tensor for a single dataset.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            %
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices. You can
            %   include a different subset of the data (or different time
            %   windows) for the alignment process separately from the actual
            %   data exported to LFADS, or return the same for both. Alignment
            %   is only relevant for multi-dataset models. If you wish to use
            %   separate data for alignment, override the method usesDifferentDataForAlignment
            %   to return true as well.
            %
            % Returns
            % ----------
            % out: a scalar struct with the following fields:
            %
            % .counts : nTrials x nChannels x nTime tensor
            %   spike counts in time bins in trials x channels x time. These
            %   should be total counts, not normalized rates, as they will be
            %   added during rebinning.
            %
            % .timeVecMs: nTime x 1 vector
            %   of timepoints in milliseconds associated with each time bin. You can start this
            %   wherever you like, but timeVecMs(2) - timeVecMs(1) will be
            %   treated as the spike bin width used when the data are later
            %   rebinned to match run.params.spikeBinMs
            %
            % .conditionId: nTrials x 1 vector
            %   of unique conditionIds. Can be cell array of strings or
            %   vector of unique integers.

            data = dataset.loadData();

            out.counts = data.spikes;
            out.timeVecMs = data.timeMs;
            out.conditionId = data.conditionId;
            out.truth = data.true_rates;
        end
    end
end
