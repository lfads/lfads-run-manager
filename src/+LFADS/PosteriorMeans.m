classdef PosteriorMeans
    % Wrapper around loaded posterior mean sampling data performed on a
    % trained LFADS model
    
    properties
        time % nTime x 1 vector of times in ms for each of the timeseries below
        controller_outputs % `nControllerOutputs x T x nTrials` controller outputs to generator inputs
        factors % `nFactors x T x nTrials` factor trajectories
        
        post_g0_mean % c_ic_enc_dim x nTrials initial conditions at encoder (before projecting out to generator units)
        post_g0_logvar % c_ic_enc_dim x nTrials initial conditions at encoder (before projecting out to generator units)
        
        generator_ics % nGeneratorUnits x nTrials` generator initial conditions
        generator_states % nGeneratorUnits x T x nTrials 
        rates % nNeurons x T x nTrials        
        costs % nTrials x 1
        nll_bound_vaes % nTrials x 1
        nll_bound_iwaes % nTrials x 1
        validInds % list of validation trial indices
        trainInds % list of training trial indices
        params % :ref:`LFADS_RunParams` instance
        
        conditionIds % nTrials x 1 vector of condition ids
        rawCounts % nNeurons x T x nTrials array of raw spike counts
    end
    
    properties(Dependent)
        isValid % contains valid, loaded data, false means empty
        nControllerOutputs
        nGeneratorUnits % number of units in the generator RNN
        nFactors % number of output factors
        nNeurons % number of actual neurons as provided in the training and validation datasets
        T % number of timepoints 
        nTrials
    end
    
    methods
        function pm = PosteriorMeans(pms, params, time, conditionIds, rawCounts)
            % pm = PosteriorMeans(pms, params, seq)
            % Construct instance by copying struct fields 
            %
            % Args:
            %   pms (struct): Posterior mean struct loaded from disk
            %   params (LFADS.RunParams): Run parameters 
            %   time (vector): time vector for all timeseries
            
            if nargin > 0
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
                        warning( sprintf( 'PosteriorMeans: couldn''t find variable %s', ...
                                          v ) );
                        pm.( v ) = [];
                    end

                end
                
                % convert rates into spikes / sec
                pm.rates = pms.rates * 1000 / params.spikeBinMs;
                % store the times
                pm.time = time;
                
            end
            if nargin > 1
                pm.params = params;
            end
            if nargin > 3
                pm.conditionIds = conditionIds;
            end
            if nargin > 4
                pm.rawCounts = rawCounts;
            end
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
            if exist(filename)
                delete(filename);
            end
            
            props = pm.listAllProperties();
            for iP = 1:numel(props)
                value = pm.(props{iP});
                if ~isnumeric(value) continue, end
                h5create(filename, ['/' props{iP}], size(value), 'DataType', class(value));
                h5write(filename, ['/' props{iP}], value);
            end
        end
        
        function props = listAllProperties(pm)
            meta = metaclass(pm);
            
            props = cell(numel(meta.PropertyList), 1);
            mask = false(numel(meta.PropertyList), 1);
            for i = 1:numel(meta.PropertyList)
                prop = meta.PropertyList(i);
                name = prop.Name;
                props{i} = name;
                
                if ~prop.Constant && ~prop.Hidden
                    mask(i) = true;
                end
            end
            
            props = props(mask);
        end
    end    
end
