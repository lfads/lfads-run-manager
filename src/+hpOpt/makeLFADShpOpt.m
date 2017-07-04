% Gets the LFADS RunCollection, defaultParams (RunParams obj) and HyperParamsSet Objects and 
% adds the parameters to LFADS RunCollection

function rc = makeLFADShpOpt(rc, defaultParams, HPSetobj)

   ParamName = {HPSetobj.ParamsSet(:).ParamName};
   for hp = HPSetobj.ParamsSamples
        for Pname = ParamName
            defaultParams.(Pname{1}) = hp.(Pname{1});
        end
        rc.addParams(defaultParams)
   end