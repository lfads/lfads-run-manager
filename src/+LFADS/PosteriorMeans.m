classdef PosteriorMeans
    % Wrapper around loaded posterior mean sampling data performed on a
    % trained LFADS model
    
    properties
        kind char = ''; % either posterior_sample_and_average or posterior_push_mean
        
        time % nTime x 1 vector of times in ms for each of the timeseries below
        controller_outputs % `nControllerOutputs x T x nTrials` controller outputs to generator inputs
        factors % `nFactors x T x nTrials` factor trajectories
        
        post_g0_mean % c_ic_enc_dim x nTrials initial conditions at encoder (before projecting out to generator units)
        post_g0_logvar % c_ic_enc_dim x nTrials initial conditions at encoder (before projecting out to generator units)
        
        generator_ics % nGeneratorUnits x nTrials` generator initial conditions
        generator_states % nGeneratorUnits x T x nTrials 
        rates % nNeurons x T x nTrials (this in spikes/sec!)  
        costs % nTrials x 1 
        nll_bound_vaes % nTrials x 1
        nll_bound_iwaes % nTrials x 1
        validInds % list of validation trial indices
        trainInds % list of training trial indices
        params % :ref:`LFADS_RunParams` instance
        
        conditionIds % nTrials x 1 vector of condition ids
        rawCounts % nNeurons x T x nTrials array of raw spike counts (this is in spikes/bin!)
        externalInputs % nExtInputs x T x nTrials array of external inputs
        
        datasetIndexByTrial % nTrials x 1 vector indicating which dataset this corresponds to
    end
    
    properties(Dependent)
        binWidthMs
        scaleFactorCountsToRates
        
        isValid % contains valid, loaded data, false means empty
        nControllerOutputs
        nGeneratorUnits % number of units in the generator RNN
        nFactors % number of output factors
        nNeurons % number of actual neurons as provided in the training and validation datasets
        T % number of timepoints 
        nTrials
    end
    
    methods
        function pm = PosteriorMeans(validFile, trainFile, validInds, trainInds, run, varargin)
            % pm = PosteriorMeans(validFileList, trainFileList, validInds, trainInds, run, time, conditionIds, rawCounts, kind)
            % Construct instance by loading the posterior means files
            %
            % Args:
            %   run (LFADS.Run): run
            %   time (vector): time vector for all timeseries
            p = inputParser();
            p.addParameter('datasetIndex', NaN, @(x) isempty(x) || isscalar(x));
            p.addParameter('time', [], @isvector);
            p.addParameter('conditionIds', [], @isvector)
            p.addParameter('rawCounts', [], @isnumeric);
            p.addParameter('externalInputs', [], @(x) isnumeric(x) || islogical(x));
            p.addParameter('kind', 'posterior_sample_and_average', @ischar);
            p.parse(varargin{:});
            
            if nargin == 0
                return;
            end
            
            % defer to utility to do the h5 loading and trial splicing
            pms = LFADS.Utils.loadPosteriorMeans(validFile, trainFile, validInds, trainInds);

            if isfield(pms, 'controller_outputs')
                pm.controller_outputs = pms.controller_outputs;
            else
                pm.controller_outputs = [];
            end
            vars_to_transfer = { 'factors', 'post_g0_mean', 'post_g0_logvar', ...
                                'generator_ics', 'generator_states', 'costs', ...
                                'nll_bound_vaes', 'nll_bound_iwaes', 'validInds', ...
                               'trainInds' };
            for nv = 1:numel( vars_to_transfer )
                v = vars_to_transfer{ nv };
                if isfield( pms, v )
                    pm.( v ) = pms. ( v );
                else
                    warning('PosteriorMeans: couldn''t find variable %s', v);
                    pm.( v ) = [];
                end
            end

            % convert rates into spikes / sec
            pm.rates = pms.rates * 1000 / run.params.spikeBinMs;
            % store the times

            if isempty(p.Results.time)
                pm.time = (1:size(pm.rates, 2)) * run.params.spikeBinMs;
            else
                pm.time = p.Results.time(1:size(pm.rates, 2));
            end
            pm.params = run.params;
                
            pm.conditionIds = p.Results.conditionIds;
            pm.rawCounts = p.Results.rawCounts;
            pm.externalInputs = p.Results.externalInputs;
            pm.kind = p.Results.kind;
            
            pm.datasetIndexByTrial = repmat(p.Results.datasetIndex, pm.nTrials, 1);
        end
    end
    
    methods
        function tf = get.isValid(pm)
            tf = ~isempty(pm.factors);
        end
        
        function n = get.nControllerOutputs(pm)
            n = size(pm.controller_outputs, 1);
        end
        
        function n = get.nGeneratorUnits(pm)
            n = size(pm.generator_states, 1);
        end
        
        function n = get.nFactors(pm)
            n = size(pm.factors, 1);
        end
        
        function n = get.nNeurons(pm)
            n = size(pm.rates, 1);
        end
        
        function n = get.T(pm)
            n = size(pm.factors, 2);
        end
        
        function n = get.nTrials(pm)
            n = size(pm.factors, 3);
        end
        
        function dt = get.binWidthMs(pm)
            dt = pm.time(2) - pm.time(1);
        end
        
        function gain = get.scaleFactorCountsToRates(pm)
            gain = 1000 / pm.binWidthMs;
        end
    end
    
    
    methods
        function avg = getConditionAveragedFieldValues(pm, field, conditionIds)
            % avg = getConditionAveragedFieldValues(field, conditionIds)
            % Using the conditionId values, average over trials (the last dim) within each
            % condition for a particular field in this class.
            %
            % Args:
            %   field (string): name of field in LFADS.PosteriorMeans class
            %   conditionIds (vector) : condition labels with size nTrials x 1. The results of
            %      unique(conditionIds) will determine the ordering of the
            %      conditions within avg
            %
            % Returns:
            %   avg (numeric):
            %     Averaged values, will have same size as the corresponding
            %     field in PosteriorMeans, except the last dim's size will
            %     be nConditions == numel(unique(conditionIds)) rather than
            %     nTrials.
            
            data = pm.(field);
            if isvector(data)
                dim = 1;
            else
                dim = ndims(data);
            end
            
            [uc, ~, whichC] = unique(conditionIds);
            nC = numel(uc);
            sz = size(data);
            sz(dim) = nC;
            
            avg = nan(sz, 'like', data);
            
            for iC = 1:nC
                if dim == 1
                    avg(iC) = mean(data(whichC == iC), dim);
                elseif dim == 2
                    avg(:, iC) = mean(data(:, whichC == iC), dim);
                elseif dim == 3
                    avg(:, :, iC) = mean(data(:, :, whichC == iC), dim);
                end
            end                      
        end
        
        function exportToHDF5(pm, filename)
            pathTo = fileparts(filename);
            if ~exist(pathTo, 'dir')
                LFADS.Utils.mkdirRecursive(pathTo);
            end
            if exist(filename, 'file')
                delete(filename);
            end
            
            props = pm.listAllProperties();
            for iP = 1:numel(props)
                value = pm.(props{iP});
                if isempty(value)
                    continue;
                end
                if islogical(value)
                    value = uint8(value);
                end
                if ~isnumeric(value), continue, end
                value = LFADS.Utils.reverseDims(value); % reorder dims since python expects row major
                h5create(filename, ['/' props{iP}], size(value), 'DataType', class(value));
                h5write(filename, ['/' props{iP}], value);
            end
        end
        
        function props = listAllProperties(pm, includeDependent)
            if nargin < 2
                includeDependent = true;
            end
            meta = metaclass(pm);
            
            props = cell(numel(meta.PropertyList), 1);
            mask = false(numel(meta.PropertyList), 1);
            for i = 1:numel(meta.PropertyList)
                prop = meta.PropertyList(i);
                name = prop.Name;
                props{i} = name;
                
                if ~prop.Constant && ~prop.Hidden && (~prop.Dependent || includeDependent)
                    mask(i) = true;
                end
            end
            
            props = props(mask);
        end
    end
    
    methods
        function pmcat = concatenateOverTrials(pms)
            constructorFn = str2func(class(pms(1)));
            pmcat = constructorFn();
            
            props = pms(1).listAllProperties(false);
            
            % find one PM with at least 2 trials so we can figure out
            % dimensions
            idxUse = NaN;
            for iPM = 1:numel(pms)
                if pms(iPM).nTrials > 1
                    idxUse = iPM;
                    break;
                end
            end
            assert(~isnan(idxUse), 'No PosteriorMeans have at least 2 trials');
            pmTest = pms(idxUse);
                        
            for iP = 1:numel(props)
                prop = props{iP};
                
                dim = find(size(pmTest.(prop)) == pmTest.nTrials, 1, 'last');
                pmcat.(prop) = cat(dim, pms.(prop));
            end
        end
        
    end
end
