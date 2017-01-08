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
       % data or high-level params
       spikeBinMs = 2; % Spike bin width in ms
       trainToTestRatio = 4; % how many train v. test trials?
                             % defaults to 4:1 ratio

       useAlignmentMatrix = false; % by default, no need to create an alignment matrix

       % DJO params
       useDJOparams = true; % default to use these params set
                            % initially by DJO. otherwise use the
                            % command-line versions specified below

       batchSize = 200; % Number of trials in one mini-batch
       regularizerIncreaseSteps = 900; % Number of steps over which the regularizer costs increase
       learningRateDecayFactor = 0.98; % Decay rate of the learning rate
       keepProb = 0.95; % Dropout of units in the network

       % added by CP
       %  these are (a subset of the) command-line params to run_lfads.py
       %  start each with c_ so it can be differentiated from other params
       c_cell_clip_value double = 5; % used to avoid stepping too far during training
       c_factors_dim uint16 = 50;
       c_in_factors_dim uint16 = 50;
       c_ic_enc_dim uint16 = 128; % network size for IC encoder
       c_ci_enc_dim uint16 = 128; % network size for controller input encoder
       c_gen_dim uint16 = 100; % generator network size
       c_keep_prob double = 0.95; % randomly drop units during each training pass
       c_learning_rate_decay_factor double = 0.98; % how quickly to decrease the learning rate
       c_device char = '/gpu:0'; % which visible GPU/CPU to use
       c_co_dim uint16 = 4;
       c_do_causal_controller logical = false; % restrict input
                                               % encoder from seeing the future?
       c_l2_gen_scale double = 500; % how much to weight the
                                    % generator l2 cost
       c_l2_con_scale double = 500; % how much to weight the
                                    % controller l2 cost
       c_batch_size uint16 = 256; % number of trials to use during
                                  % each training pass
       c_kl_increase_steps uint16 = 900; % Number of steps over which the kl costs increase
       c_l2_increase_steps uint16 = 900; % Number of steps over which the l2 costs increase
       c_controller_input_lag uint16 = 1; 
       c_ic_dim uint16 = 64; % dimensionality of the initial conditions
       c_con_dim uint16 = 128; %controller dimensionality
   end


   methods
       function str = generateCommandLineOptionsString(p)
       % function str = generateCommandLineOptionsString(p)
       %
       %   generates a string of all the command line options to
       %   pass directly into run_lfads.py

           % get all the parameters
           f = fields(p);
           % only keep the command line params
           %  these begin with 'c_'
           keepFields = false(numel(f), 1);
           for nf=1:numel(f)
               if numel(f{nf}) > 2 && strcmp(f{nf}(1:2), 'c_')
                   keepFields(nf) = true;
               end
           end
           f = f(keepFields);

           % build the output string
           str = '';
           for nf = 1:numel(f)
               clear fieldstr
               thisField = f{nf};
               thisVal = p.(thisField);
               
               % argument formatted differently for each class
               switch class(thisVal)
                 case {'uint16', 'uint32', 'int16', 'int32'}
                   fieldstr = sprintf('%i', thisVal);
                 case {'logical'}
                   if thisVal, fieldstr = 'true', 
                   else fieldstr = 'false'; end
                 case {'double','single'}
                   fieldstr = sprintf('%.3f', thisVal);
                 case {'char'}
                   fieldstr = thisVal;
                 otherwise
                   error(['don''t know this type: ' class(thisVal)]);
               end
               str = sprintf('%s --%s=%s', str, thisField(3:end), ...
                             fieldstr);
           end
           
       end
   end

end
