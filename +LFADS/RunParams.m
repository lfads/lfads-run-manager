classdef RunParams < handle
   properties
       spikeBinMs = 2;
       batchSize = 200; % trials
       regularizerIncreaseSteps = 900;
       learningRateDecayFactor = 0.98;
       keepProb = 0.95; % dropout of units in the network
   end
   
   methods(Abstract)
        suffix = generateSuffix(p);
   end
end