classdef RunParams < LFADS.RunParams
    properties
        c_do_subpop_readin logical = false;
        c_do_subpop_readout logical = false;
        
        numSubfactorsBySubpopulation uint16; % vector nSubpopulations x 1 of subfactor counts
    end
    
    methods
        function par = RunParams(varargin)
            par@LFADS.RunParams(varargin{:});
        end
    end
    
    methods(Access=protected)
        function groups = generateStandardPropertyGroups(obj)
            groupTitle = 'Two-Stage Readin and Readout';
            propList = {'c_do_subpop_readin', 'c_do_subpop_readout', 'numSubfactorsBySubpopulation'};
            insert = matlab.mixin.util.PropertyGroup(propList, groupTitle);

            groups = generateStandardPropertyGroups@LFADS.RunParams(obj);
            groups = [insert, groups];
        end
    end
end
