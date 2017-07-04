classdef HyperParam < handle
    properties
        ParamName               % Name of the parameter
        ParamType               % Type of the parameter, e.g. int, cont
        ParamDistrType          % Distribution type
        ParamDistrParams        % Parameters for distribution
        ParamValueVec           % Vectors
    end
    
    methods
        function HP = HyperParam(varargin)
            p = inputParser();
            p.addParameter('ParamName', [], @ischar)
            p.addParameter('ParamType', 'cont', @ischar)
            p.addParameter('DistrType', 'uniform', @ischar)
            p.addParameter('DistrParam', [])
            p.addParameter('ValueVec', [])
            p.parse(varargin{:});

            HP.ParamName = p.Results.ParamName;
            HP.ParamType = p.Results.ParamType;
            HP.ParamDistrType = p.Results.DistrType;  
            HP.ParamDistrParams = p.Results.DistrParam; 
            HP.ParamValueVec = p.Results.ValueVec;
        end
        
        function Samples = sample(HP, nsample)
                switch HP.ParamDistrType
                    case 'copy'
                        % This will be taked care of in HyperParamSet
                        % (higher level)
                        % For setting a parameter as the copy of another one set: 
                        % DistrType = 'copy',
                        % DistrParam = <name of the source parameter to be copied>
                    case 'constant'
                        V = HP.ParamValueVec;
                        if length(V) > 1
                           Samples = V; % use the input vector by the user
                        else
                            Samples = V * ones(1, nsample);
                        end
                        
                    case 'grid'
                        numGrids = HP.ParamDistrParams;
                        Vmin = HP.ParamValueVec(1);
                        Vmax = HP.ParamValueVec(2);
                        Samples = linspace(Vmin, Vmax, numGrids);
                        
                    case 'uniform'
                        Vmin = HP.ParamValueVec(1);
                        Vmax = HP.ParamValueVec(2);
                        % Generate samples from unifrom distribution in
                        % [0,1]
                        rndSamples = rand(1, nsample);
                        
                        % Shift and scale to match the desired range
                        Samples = (Vmax - Vmin) * rndSamples + Vmin;                        
                        
                    case 'loguniform'
                        ePar = HP.ParamDistrParams;
                        Vmin = HP.ParamValueVec(1);
                        Vmax = HP.ParamValueVec(2);
                        % Generates samples from log uniform probability 
                        % distribution in [0,1]
                        rndSamples = exp(ePar * rand(1, nsample)) / (exp(ePar) + exp(0));
                        
                        % Shift and scale to match the desired range
                        Samples = (Vmax - Vmin) * rndSamples + Vmin;  
                        
                    otherwise
                        error('Undefined Parameter Distribution Type! Choose one of these values: ''constant'', ''grid'', ''uniform'', ''loguniform''');
                end
                if strcmp(HP.ParamType, 'int')
                    Samples = round(Samples);
                end
        end
    end
    
end