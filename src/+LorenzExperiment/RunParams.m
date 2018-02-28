classdef RunParams < LFADS.RunParams
   properties
       % Add additional custom parameters here. The default you assign to
       % them will be used when computing the hash value. Any params whose value
       % differs from the default will be included in the hash value, to allow new
       % parameters to be added without invalidating old hashes. So choose
       % the default once and don't change it. If you decide to use another
       % value later by default, override it in the constructor instead.
   end
   
   methods
       function par = RunParams()
          % adjust whatever defaults you like. Modificiations made here
          % function identically to modifications made by the user, meaning
          % that the value will be incorporated in the hash only if the
          % value differs from the property definition in the class.
          
       end 
   end
end