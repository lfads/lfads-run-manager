classdef RunParams < handle
    % Collection of parameters which are common to all :ref:`LFADS_Run` instances in a :ref:`LFADS_RunCollection`. You
    % must create a subclass of RunParams in which you specify the serialized representation of the parameters that
    % will be used in paths on disk.

    methods(Abstract)
         suffix = generateSuffix(p);
         % Generates a string representation of the parameters that is used for reading off the parameters from a
         % folder path on disk. You can include as many or as few fields in this serialized representation as you see
         % fit.
         %
         % Returns
         % ---------
         % suffix : string
         %   serialized string suffix suitable for inclusion in file paths
    end

   properties
       spikeBinMs = 2; % Spike bin width in ms
       batchSize = 200; % Number of trials in one mini-batch
       regularizerIncreaseSteps = 900; % Number of steps over which the regularizer costs increase
       learningRateDecayFactor = 0.98; % Decay rate of the learning rate
       keepProb = 0.95; % Dropout of units in the network

       % by default, no need to create an alignment matrix
       useAlignmentMatrix = false;

       % added by CP
       gen_dim = 100; % generator network size
       ci_enc_dim = 128; % network size for controller input encoder
       ic_enc_dim = 128; % network size for IC encoder
       con_dim = 128; %controller dimensionality
       factors_dim = 30;
       co_dim = 4;
       ic_dim = 64;
       in_factors_dim = 50; % i don't beieve this is used in non-multiday 
   end


end
