classdef HyperParamsSet < handle
    % Create the HyperParamsSet object
    % For setting a parameter as the copy of another one set: DistrType = 'copy',
    % DistrParam = <name of the source parameter to be copied>
    properties
        ParamsSet
        ParamsSamples
    end
    
    methods
        function HPs = HyperParamsSet()
            HPs.ParamsSet = [];                 % Name of the parameter
            HPs.ParamsSamples = [];
        end
        
        function addHyperParam(HPs, varargin)
            newHP = hpOpt.HyperParam(varargin{:});
            
            HPs.ParamsSet = cat(1, HPs.ParamsSet, newHP);
        end
        
        function sampleSet(HPs, nsample)
            % Sample each HyperParameter
            for hp = HPs.ParamsSet(:)'
                paramSamp(nsample).(hp.ParamName) = [];     % initialize the array
                if strcmp(hp.ParamDistrType, 'copy')     % If the param is a copy of another param
                    hpVals = [paramSamp(:).(hp.ParamDistrParams)];
                else                                % Else sample from Distribution
                    hpVals = hp.sample(nsample);
                end
                for i = 1:length(hpVals)
                    paramSamp(i).(hp.ParamName) = hpVals(i);
                end
            end
            HPs.ParamsSamples = paramSamp;
        end
        
        % Make Grid
        function makeGrid(HPs)
            i = 0;
            for hp = HPs.ParamsSet(:)'
                if ~strcmp(hp.ParamDistrType, 'copy')
                    i = i + 1;
                    hpVals{i} = hp.sample([]);
                end
            end
            
            cmdStr = '';
            for i = 1:length(hpVals)
                cmdStr = sprintf('%s hpVals{%i},', cmdStr, i);
            end
            
            cmdStr = ['combvec(' cmdStr(1:end-1) ')'];
            paramVecSamp = eval(cmdStr);
            
            nSamp = size(paramVecSamp, 2);
            j = 0;
            for hp = HPs.ParamsSet(:)'
                j = j + 1;
                paramSamp(nSamp).(hp.ParamName) = [];     % initialize the array
                if strcmp(hp.ParamDistrType, 'copy')     % If the param is a copy of another param
                    hpVals = [paramSamp(:).(hp.ParamDistrParams)];
                else                                % Else sample from params vector
                    hpVals = paramVecSamp(j, :);
                end
                
                for i = 1:nSamp
                    paramSamp(i).(hp.ParamName) = hpVals(i);
                end
            end
            HPs.ParamsSamples = paramSamp;
            
        end
    end
    
end