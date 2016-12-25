classdef PosteriorMeans
    % Wrapper around loaded posterior mean sampling data performed on a
    % trained LFADS model
    
    properties
        controller_outputs % `nControllerOutputs x T x nTrials` controller outputs to generator inputs
               factors % `nFactors x T x nTrials` factor trajectories
         generator_ics % nGeneratorUnits x nTrials` generator initial conditions
      generator_states % nGeneratorUnits x T x nTrials 
                 rates % nNeurons x T x nTrials
             validInds % list of validation trial indices
             trainInds % list of training trial indices
                params % :ref:`LFADS_RunParams` instance
    end
    
    properties
        isValid % contains valid, loaded data, false means empty
        nControllerOutputs
        nGeneratorUnits % number of units in the generator RNN
        nFactors % number of output factors
        nNeurons % number of actual neurons as provided in the training and validation datasets
        T % number of timepoints 
        nTrials
    end
    
    methods
        function pm = PosteriorMeans(pms, params)
            % Construct instance by copying struct fields 
            %
            % Args:
            %   pms (struct) : Posterior mean struct loaded from disk
            %   params (ref:`LFADS_RunParams`) : Run parameters 
            
            if nargin > 0
                pm.controller_outputs = pms.controller_outputs;
                pm.factors = pms.factors;
                pm.generator_ics = pms.generator_ics;
                pm.generator_states = pms.generator_states;
                pm.rates = pms.rates;
                pm.validInds = pms.validInds;
                pm.trainInds = pms.trainInds;
            end
            if nargin > 1
                pm.params = params;
            end
        end
    end
    
    methods
        function tf = get.isValid(pm)
            tf = ~isempty(pm.controller_outputs);
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
            n = size(pm.controller_outputs, 3);
        end
    end
end
