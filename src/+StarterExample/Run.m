classdef Run < LFADS.Run
    methods
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end

        function convertDatasetToSequenceStruct(r, dataset, mode, varargin)
            % Overwrite this to convert dataset into the format required by
            % sequence data, with fields `y`, `y_time`, `binWidthMs`, and
            % (optionally) `conditionId`.
            %
            % Converts the loaded data within a dataset into a sequence struct. The sequence data returned will be a
            % struct array where each element corresponds to a trial. You can include any metadata or information fields
            % that you like for future reference or analysis. At a minimum, you must include field `.y`, `.y_time`, and `.params.dtMS`
            %
            % For each trial:
            %   - `.y` will contain spike data as an `nNeurons` x `nTimeBins` array of spike counts. The spike binning
            %       for this is determined by you, and can be left at 1 ms. Later,
            %       this data will be rebinned according to RunParams .spikeBinMs field.
            %   - `.y_time` provides a time vector corresponding to the time
            %       bins of `.y`, which should be identical on each trial
            %   - `.binWidthMs` specifies the time bin width for `.y`
            %   - optionally: `.conditionId` specifies the condition to which each trial
            %       belongs. This information isn't passed to LFADS. It is used
            %       only when building the alignment matrices for multi-session
            %       stitching, if trial-averaging is employed.
            %
            % Parameters
            % ------------
            % dataset : :ref:`LFADS_Dataset`
            %   The :ref:`LFADS_Dataset` instance from which data were loaded
            % mode (string) : typically 'export' indicating sequence struct
            %   will be exported for LFADS, or 'alignment' indicating that this
            %   struct will be used to generate alignment matrices
            %
            % Returns
            % ----------
            % seq : struct Array
            %   sequence formatted data. A struct array where each elemnt corresponds to a specific trial.

%            nTrials = numel(dataset);

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
