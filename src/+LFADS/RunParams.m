classdef RunParams < matlab.mixin.CustomDisplay
    % Collection of parameters which are common to all :ref:`LFADS_Run` instances in a :ref:`LFADS_RunCollection`. You
    % must create a subclass of RunParams in which you specify the serialized representation of the parameters that
    % will be used in paths on disk.
    %
    % When you add additional parameters to this list, don't prefix yours
    % with c_, as these have special meaning to the LFADS code and will
    % mess things up.

    % Not yet implemented:
    %  do_train_io_only - needs more integration into run manager as a post-hoc run
    %  do_train_encoder_only - - needs more integration into run manager as a post-hoc run
    %  do_reset_learning_rate - not sure it makes sense here since it would
    %    be for post hoc runs only

   methods(Access = protected)
       function groups = generateStandardPropertyGroups(obj) %#ok<MANU>
             % this is separate to allow subclasses to override this method
             groupTitle = 'Computed hashes';
             propList = {'paramHash', 'paramHashString', 'dataHash', 'dataHashString'};
             groups(1) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Run Manager logistics and data processing';
             propList = {'name', 'version', 'spikeBinMs'};
             groups(2) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'TensorFlow Logistics';
             propList = {'c_allow_gpu_growth', 'c_max_ckpt_to_keep', 'c_max_ckpt_to_keep_lve', 'c_device'};
             groups(3) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Optimization';
             propList = {'c_learning_rate_init', 'c_learning_rate_decay_factor', 'c_learning_rate_n_to_compare', 'c_learning_rate_stop', 'c_max_grad_norm', 'trainToTestRatio', 'c_batch_size', 'c_cell_clip_value'};
             groups(4) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Overfitting';
             propList = {'c_temporal_spike_jitter_width', 'c_keep_prob', 'c_l2_gen_scale', 'c_l2_con_scale', 'c_co_mean_corr_scale'};
             groups(5) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Underfitting';
             propList = {'c_kl_ic_weight', 'c_kl_co_weight', 'c_kl_start_step', 'c_kl_increase_steps', 'c_l2_start_step', 'c_l2_increase_steps', 'scaleIncreaseStepsWithDatasets', 'scaleStartStepWithDatasets'};
             groups(6) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'External inputs';
             propList = {'c_ext_input_dim', 'c_inject_ext_input_to_gen'};
             groups(7) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Controller and inferred inputs';
             propList = {'c_co_dim', 'c_prior_ar_atau', 'c_do_train_prior_ar_atau', 'c_prior_ar_nvar', 'c_do_train_prior_ar_nvar', 'c_do_causal_controller', 'c_do_feed_factors_to_controller', ...
               'c_feedback_factors_or_rates', 'c_controller_input_lag', 'c_ci_enc_dim', 'c_con_dim', ...
               'c_co_prior_var_scale'};
             groups(8) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Encoder and initial conditions for generator';
             propList = {'c_num_steps_for_gen_ic', 'c_ic_dim', 'c_ic_enc_dim', 'c_ic_prior_var_min', 'c_ic_prior_var_scale', 'c_ic_prior_var_max', 'c_ic_post_var_min'};
             groups(9) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Generator network, factors, rates';
             propList = {'c_cell_weight_scale', 'c_gen_dim', 'c_gen_cell_input_weight_scale', ...
               'c_gen_cell_rec_weight_scale', 'c_factors_dim', 'c_output_dist'};
             groups(10) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Stitching multi-session models';
             propList = {'c_do_train_readin', 'useAlignmentMatrix', 'useSingleDatasetAlignmentMatrix', 'alignmentApproach', 'alignmentExtraArgs'};
             groups(11) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

             groupTitle = 'Posterior sampling';
             propList = {'posterior_mean_kind', 'num_samples_posterior'};
             groups(12) = matlab.mixin.util.PropertyGroup(propList,groupTitle);
             
             groupTitle = 'TensorFlow Debugging';
             propList = {'c_tf_debug_cli', 'c_tf_debug_tensorboard', 'c_tf_debug_tensorboard_hostport', ...
                'c_tf_debug_dump_root', 'c_debug_verbose', 'c_debug_reduce_timesteps_to', ...
                'c_debug_print_each_step'};
             groups(13) = matlab.mixin.util.PropertyGroup(propList,groupTitle);
       end
       
       function groups = mergePropertyGroups(obj, groupsOld, groupsNew) %#ok<INUSL>
            % for each group in groupsNew
            %   if found in groupsOld, add properties to end of corresponding property lists
            %   if not found, add to beginning of property groups
   
            oldTitles = {groupsOld.Title};
            
            alreadyExist = false(numel(groupsNew), 1);
            for iN = 1:numel(groupsNew)
                [alreadyExist(iN), idx] = ismember(groupsNew(iN).Title, oldTitles);
                if alreadyExist(iN)
                    groupsOld(idx).PropertyList = cat(2, groupsOld(idx).PropertyList, groupsNew(iN).PropertyList);
                end 
            end
            
            groups = [groupsNew(~alreadyExist), groupsOld];
       end
       
       function groups = getPropertyGroups(obj)
          if isscalar(obj)
             groups = obj.generateStandardPropertyGroups();
             
             propsPrinted = arrayfun(@(grp) grp.PropertyList', groups, 'UniformOutput', false);
             propsPrinted = cat(1, propsPrinted{:});
             
             customProps = obj.listNonTransientProperties('ignoreProperties', propsPrinted, ...
                'onlyRootClassProperties', false, ...
                'ignoreHidden', true);

             groupTitle = sprintf('%s Specific Properties', class(obj));
             customGroup = matlab.mixin.util.PropertyGroup(customProps,groupTitle);
             groups = [customGroup, groups];
          else
             % Nonscalar case: call superclass method
             groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
          end
       end
    end

    properties
        % Command-line params to run_lfads.py start with c_ so they can be
        % differentiated from other params internal to the run manager.
        % Avoid naming other parameters the same as these without the c_
        % since the prefix will be stripped from the serialization

        % These properties DO NOT affect the param_HASH or the data_HASH
        name char = ''; % convenient name for displaying this param
        version uint32 = 20171107; % Used for graceful evolution of path settings

        % Data processing
        spikeBinMs double = 2; % Spike bin width in ms

        % TensorFlow logistics
        c_allow_gpu_growth logical = true; %whether to allow the GPU to dynamically allocate memory. default (false) is to allocate all the memory initially
        c_max_ckpt_to_keep uint16 = 5; % Max # of checkpoints to keep (rolling)
        c_max_ckpt_to_keep_lve uint16 = 5; % Max # of checkpoints to keep for lowest validation error models (rolling)
        c_device char = '/gpu:0'; % which visible GPU/CPU to use

        % Optimization
        c_learning_rate_init double = 0.01; % Initial learning rate
        c_learning_rate_decay_factor double = 0.98; % how quickly to decrease the learning rate
        c_learning_rate_n_to_compare uint16 = 6; % "Number of previous costs current cost has to be worse than, in order to lower learning rate
        c_learning_rate_stop double = 0.00001; % when the learning rate reaches this threshold, stop training
        c_max_grad_norm double = 200; % Max norm of gradient before clipping
        trainToTestRatio uint16 = 4; % how many train v. test trials, defaults to 4:1 ratio
        c_batch_size uint16 = 256; % number of trials to use during each training pass
        c_cell_clip_value double = 5; % Max value recurrent cell can take before being clipped

        % Overfitting
        c_temporal_spike_jitter_width uint16 = 0; % jittering spike times during training, in number of spike bins
        c_keep_prob double = 0.95; % randomly drop units during each training pass
        c_l2_gen_scale double = 500; % L2 regularization cost for the generator only.
        c_l2_con_scale double = 500; % L2 regularization cost for the controller only.
        c_co_mean_corr_scale double = 0 % Cost of correlation (thru time)in the means of controller output

        % Underfitting:
        c_kl_ic_weight double = 1; % Strength of KL weight on initial conditions KL penalty
        c_kl_co_weight double = 1; % Strength of KL weight on controller output KL penalty
        c_kl_start_step uint32 = 0; % Start increasing KL weight after this many steps
        c_kl_increase_steps uint32 = 900; % Number of steps over which the kl costs increase
        c_l2_start_step uint32 = 0; % Start increasing L2 weight after this many steps
        c_l2_increase_steps uint32 = 900; % Number of steps over which the l2 costs increase
        scaleIncreaseStepsWithDatasets logical = true; % If true, c_kl_increase_steps and c_l2_increase_steps will be multiplied by the number of datasets in a stitching Run
        scaleStartStepWithDatasets logical = true; % If true, c_kl_start_step and c_l2_start_step will be multiplied by the number of datasets in a stitching Run

        % External inputs`
        c_ext_input_dim uint16 = 0; % Number of external inputs
        c_inject_ext_input_to_gen logical = false; % should observed inputs be input to model via encoders (false) or injected directly into generator (true)?

        % Controller / Inferred inputs
        c_co_dim uint16 = 4; % number of inferred inputs (controller outputs)
        c_prior_ar_atau double = 10; % Initial autocorrelation of AR(1) priors (in time bins)
        c_do_train_prior_ar_atau logical = true; % Is the value for atau an initial value (true) or the constant value (false)
        c_prior_ar_nvar double = 0.1; % Initial noise variance for AR(1) priors
        c_do_train_prior_ar_nvar logical = true; % Is the value for the noise var an initial value (true) or the constant value (false)
        c_do_causal_controller logical = false; % restrict input encoder from seeing the future?DO_FEED_FACTORS_TO_CONTROLLER
        c_do_feed_factors_to_controller logical = true; %  Feed the factors back to the controller?
        c_feedback_factors_or_rates char = 'factors'; % Feedback the factors or the rates to the controller? Acceptable values: 'factors' or 'rates'
        c_controller_input_lag uint16 = 1;
        c_ci_enc_dim uint16 = 128; % network size for controller input encoder
        c_con_dim uint16 = 128; %controller dimensionality
        c_co_prior_var_scale = 0.1; % Variance of control input prior distribution

        % Encoder / Generator initial condition
        c_num_steps_for_gen_ic uint32 = intmax('uint32'); % Number of steps to train the generator initial condition.
        c_ic_dim uint16 = 64; % dimensionality of the initial conditions
        c_ic_enc_dim uint16 = 128; % network size for IC encoder
        c_ic_prior_var_min double = 0.1; % Minimum variance of IC prior distribution
        c_ic_prior_var_scale double = 0.1; % Variance of IC prior distribution
        c_ic_prior_var_max double = 0.1 % Maximum variance of IC prior distribution
        c_ic_post_var_min double = 0.0001; % Minimum variance of IC posterior distribution

        % Generator network
        c_cell_weight_scale double = 1.0; % Input scaling for input weights in generator
        c_gen_dim uint16 = 100; % generator network size
        c_gen_cell_input_weight_scale double = 1.0; % Input scaling for input weights in generator, which will be divided by sqrt(#inputs)
        c_gen_cell_rec_weight_scale double = 1.0; % Input scaling for recurrent weights in generator.

        % Factors from generator
        c_factors_dim uint16 = 50;

        % Output distribution for rates
        c_output_dist char = 'poisson'; % or gaussian, type of output distribution

        % Stitching and alignment matrices
        c_do_train_readin logical = true; % for stitching models, make the readin matrices trainable (true) or fix them to equal the alignment matrices (false)
        useAlignmentMatrix logical = false; % Whether to use an alignment matrix when stitching datasets together.
        useSingleDatasetAlignmentMatrix logical = false;  % Whether to use an alignment matrix using a single dataset, for dimensionality reduction upstream of the encoder
        alignmentApproach char = 'regressGlobalPCs'; % one of 'regressGlobalPCs' or 'ridgeRegressGlobalPCs'. see LFADS.Run/prepareAlignmentMatrices
        alignmentExtraArgs cell = {}; % extra args to pass to MutlisessionAlignmentTool

        % Posterior mean sampling
        % These properties DO NOT affect the param_HASH or the data_HASH
        posterior_mean_kind char = 'posterior_sample_and_average'; % or 'posterior_push_mean'
        num_samples_posterior = 512; % number of samples when using posterior_sample_and_average
        
        % Debugging (not included in hash)
        c_tf_debug_cli logical = false;
        c_tf_debug_tensorboard logical = false;
        c_tf_debug_tensorboard_hostport char = 'localhost:6064'
        c_tf_debug_dump_root = ''
        c_debug_verbose = false
        c_debug_reduce_timesteps_to = []
        c_debug_print_each_step = false
    end

    % Retired properties that should be kept around for hash value purposes
    % but no longer output to LFADS
    properties(Hidden)
        c_in_factors_dim uint16 = 50;
        setInFactorsMatchDataForSingleDataset logical = false; % if true, c_in_factors_dim will be set to the dimensionality of the data when only a single dataset is used
    end

    properties(Dependent)
        paramHash
        paramHashString
        dataHash
        dataHashString
    end

    methods
        function paramArray = generateSweep(p, varargin)
            % Generates an array of RunParams objects that sweep the
            % specified parameter values.
            %
            % Usage: params = p.generateSweep('prop1', valueList1, 'prop2', valueList2, ...);

            assert(mod(numel(varargin), 2) == 0, 'Inputs must be ''propertyName'', value pairs');

            nProp = numel(varargin) / 2;
            props = varargin(1:2:end);
            vals = varargin(2:2:end);
            outSz = cellfun(@numel, vals);
            if nProp == 1
                outSz(2) = 1;
            end

            % check for valid properties
            assert(iscellstr(props),  'Inputs must be ''propertyName'', value pairs');
            allProps = properties(p);
            for iProp = 1:nProp
                assert(ismember(props{iProp}, allProps), '%s is not a valid property of %s', props{iProp}, class(p));
                assert(isvector(vals{iProp}), 'Value list for %s is not a vector', props{iProp});
            end

            paramArray = repmat(p, outSz);

            valInds = cell(nProp, 1);

            for iParam = 1:numel(paramArray)
                [valInds{:}] = ind2sub(outSz, iParam);
                for iProp = 1:nProp
                    if iscell(vals{iProp})
                        val = vals{iProp}{valInds{iProp}};
                    else
                        val = vals{iProp}(valInds{iProp});
                    end
                    paramArray(iParam).(props{iProp}) = val;
                end
            end
        end

        function tf = eq(a, b)
            % Overloaded == operator to enable equality if name,
            % datasetCollection, and datasets all match

            tf = isequal(a.generateHash, b.generateHash);
        end

        function [props, propMeta] = listNonTransientProperties(p, varargin)
            % Return a list of properties in this class that are not
            % Dependent, Constant, Transient, or Hidden
            %
            % Args:
            %   ignoreProperties : cellstr
            %     list of properties to ignore
            %   onlyRootClassProperties : bool (False)
            %     if true, only include properties declared in
            %     LFADS.RunParams. if false, include properties declared in
            %     subclasses.
            %
            % Returns:
            %   props : cellstr
            %     list of properties
            %   propMeta : meta.property array
            %      meta.property list from metaclass

            parser = inputParser();
            parser.addParameter('ignoreProperties', {}, @iscellstr);
            parser.addParameter('onlyRootClassProperties', false, @islogical);
            parser.addParameter('ignoreHidden', false, @islogical);

            parser.parse(varargin{:});

            if parser.Results.onlyRootClassProperties
                meta = ?LFADS.RunParams;
            else
                meta = metaclass(p);
            end

            % always exclude these fields
            alwaysIgnore = {'version', 'name'};

            propertyList = meta.PropertyList;
            props = cell(numel(propertyList), 1);
            mask = false(numel(propertyList), 1);
            for i = 1:numel(propertyList)
                prop = propertyList(i);
                name = prop.Name; %#ok<*PROPLC>
                props{i} = name;
                if ismember(name, parser.Results.ignoreProperties)
                    continue;
                end
                if ismember(name, alwaysIgnore)
                    continue;
                end
                % skip properties that are Dependent, Constant, Transient,
                % or Hidden. Serialize the value if it differes from the
                if ~prop.Dependent && ~prop.Constant && ~prop.Transient && (~parser.Results.ignoreHidden || ~prop.Hidden)
                    mask(i) = true;
                end
            end

            props = props(mask);
            propMeta = propertyList(mask);
        end

        function out = getPropertyValueSubset(p, varargin)
            % Generates a struct of non-transient property values that
            % differ from their default values.
            %
            % Args:
            %   ignoreProperties (cellstr) : list of properties to ignore
            %   filterDiffersFromDefault (bool=true) :
            %     if true, only properties whose value differs from the
            %     default value
            %   onlyRootClassProperties (bool = false): if true, only include properties declared in
            %     LFADS.RunParams. if false, include properties declared in
            %     subclasses.
            %   defaultsFromClassDefinition (bool = true):
            %     if true, default value comes from the class definition,
            %     next to each property's definition. If false, the default
            %     values come from constructing a new class of the same
            %     type with no arguments. Default true.
            % Returns:
            %   differ (struct) : property names and values that differ

            parser = inputParser();
            parser.addParameter('onlyDifferentFromDefault', false, @islogical);
            parser.addParameter('ignoreProperties', {}, @iscellstr);
            parser.addParameter('onlyRootClassProperties', false, @islogical);
            parser.addParameter('defaultsFromClassDefinition', true, @islogical);
            parser.addParameter('ignoreHiddenUnlessDifferentFromDefault', false, @islogical);
            parser.addParameter('omitFields', {}, @iscellstr);

            parser.parse(varargin{:});

            [props, propMeta] = p.listNonTransientProperties('ignoreProperties', parser.Results.ignoreProperties, ...
                'onlyRootClassProperties', parser.Results.onlyRootClassProperties, ...
                'ignoreHidden', false); % get hidden props too

            if ~parser.Results.defaultsFromClassDefinition
                defaultInstance = eval(class(p));
            end

            out = struct();
            for i = 1:numel(props)
                prop = props{i};
                value = p.(prop);

                % skips prop equal to default (either all props or just
                % hidden on
                if (propMeta(i).Hidden && parser.Results.ignoreHiddenUnlessDifferentFromDefault) || parser.Results.onlyDifferentFromDefault
                    if parser.Results.defaultsFromClassDefinition
                        if propMeta(i).HasDefault
                            def = propMeta(i).DefaultValue;
                        else
                            def = [];
                        end
                    else
                        def = defaultInstance.(prop);
                    end
                    if isequal(def, value)
                        continue;
                    end
                end
                out.(prop) = value;
            end
        end

        function list = getListPropertiesNotAffectingInputDataHash(p)
            % list non-transient properties that do not affect the
            % data_HASH. currently includes all c_* properties except
            % c_factors_dim as well as the properties listed by getListPropertiesNotAffectingHash

            props = p.listNonTransientProperties('onlyRootClassProperties', true);

            if p.version < 20171107
                mask = cellfun(@(x) strncmp('c_', x, 2), props);
            else
                % c_factors_dim now affects the data hash because it
                % changes the saved alignment matrix
                doesAffectInputData = {'c_factors_dim'};
                mask = cellfun(@(x) strncmp('c_', x, 2) && ~ismember(x, doesAffectInputData), props);
            end
            list = props(mask);
            list = union(list, p.getListPropertiesNotAffectingHash());
        end

        function list = getListPropertiesNotAffectingHash(p) %#ok<MANU>
            % this provides a list of all properties in the class that
            % should not affect the resulting param_HASH, regardless of
            % their values.
            list = {'version', 'num_samples_posterior', 'posterior_mean_kind', ...
            'c_tf_debug_cli', 'c_tf_debug_tensorboard', 'c_tf_debug_tensorboard_hostport', ...
            'c_tf_debug_dump_root', 'c_debug_verbose', 'c_debug_reduce_timesteps_to', ...
            'c_debug_print_each_step'};
        end
        
        function list = getListPropertiesScaleWithNumDatasets(p)
            if p.scaleIncreaseStepsWithDatasets
                list1 = {'c_kl_increase_steps', 'c_l2_increase_steps'};
            else
                list1 = {};
            end
            if p.scaleStartStepWithDatasets
                list2 = {'c_kl_start_step', 'c_l2_start_step'};
            else
                list2 = {};
            end
            list = union(list1, list2);
        end

        function hash = generateHash(p)
            % Generate a short hash of this RunParams non-transient
            % properties that can be used in a directory name.
            %
            % Args:
            %   length (int): number of characters to truncate the hash value to
            %
            % Returns:
            %   hash (string): hash string

            %parser = inputParser();
            %parser.addParameter('length', 6, @isscalar);
            %parser.KeepUnmatched = true;
            %parser.parse(varargin{:});

            data = p.generateStructForHash();
            length = 6;
            
            % we include AutoDownsizeIntegers=true to allow for graceful
            % expansion of integer types (some used to be 16 bit that are
            % now 32 bit)
            hash = LFADS.Utils.DataHash(data, struct('Format', 'base64', 'AutoDownsizeIntegers', true));
            hash = strrep(strrep(hash, '/', '_'), '+', '-'); % https://tools.ietf.org/html/rfc3548#page-6
            if numel(hash) > length
                hash = hash(1:length);
            end
        end
        
        function data = generateStructForHash(p)
            ignore = p.getListPropertiesNotAffectingHash();
            data = p.getPropertyValueSubset('ignoreProperties', ignore, 'onlyDifferentFromDefault', true);
        end

        function hash = generateInputDataHash(p)
            % Generate a short hash of this RunParams non-transient
            % properties that includes only those properties that affect
            % the data preprocessing that are fed into LFADS, excluding
            % parameters that affect the internal operation of LFADS only.
            % This is used to generate the input directory
            %
            % Args:
            %   length (int): number of characters to truncate the hash value to
            %
            % Returns:
            %   hash (string): hash string

            %parser = inputParser();
            %parser.addParameter('length', 6, @isscalar);
            %parser.KeepUnmatched = true;
            %parser.parse(varargin{:});

            length = 6;
            data = p.generateStructForInputDataHash();
            hash = LFADS.Utils.DataHash(data, struct('Format', 'base64', 'AutoDownsizeIntegers', true));
            hash = strrep(strrep(hash, '/', '_'), '+', '-'); % https://tools.ietf.org/html/rfc3548#page-6
            if numel(hash) > length
                hash = hash(1:length);
            end
        end
        
        function data = generateStructForInputDataHash(p)
            propsIgnore = p.getListPropertiesNotAffectingInputDataHash();
            data = p.getPropertyValueSubset('ignoreProperties', propsIgnore, 'onlyDifferentFromDefault', true);
        end

        function str = generateHashName(p)
            % Generate a short hash name like 'param_HASH' where HASH is
            % generated by generateHash()

            str = sprintf('param_%s', p.generateHash());
        end

        function str = generateInputDataHashName(p)
            % Generate a short hash name like 'data_HASH' where HASH is
            % generated by generateInputDataHash()

            str = sprintf('data_%s', p.generateInputDataHash());
        end

        function str = generateString(p, varargin)
            % Generates a string representation of all parameters, with custom
            % strings inserted between property names and values.
            %
            % Args:
            %   beforeProp : char
            %     default ''
            %   betweenPropValue : char
            %     default '='
            %   afterValue : char
            %     default ''
            %   betweenProps : char
            %     default ' '
            %   ignoreProperties : cellstr
            %     list of properties to ignore
            %   onlyRootClassProperties : bool (False)
            %     if true, only include properties declared in
            %     LFADS.RunParams. if false, include properties declared in
            %     subclasses.
            %
            % Also accepts arguments for getPropertyValueSubset
            %
            % Returns:
            %   str : string
            %     serialized string suffix suitable for inclusion in file paths
            %

            parser = inputParser();
            parser.addParameter('beforeProp', '', @ischar);
            parser.addParameter('betweenPropValue', '=', @ischar);
            parser.addParameter('afterValue', '', @ischar);
            parser.addParameter('betweenProps', ' ', @ischar);
            parser.addParameter('onlyDifferentFromDefault', false, @islogical);
            parser.addParameter('ignoreHiddenUnlessDifferentFromDefault', true, @islogical);
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});

            data = p.getPropertyValueSubset('onlyDifferentFromDefault', parser.Results.onlyDifferentFromDefault, ...
                'ignoreHiddenUnlessDifferentFromDefault', parser.Results.ignoreHiddenUnlessDifferentFromDefault, parser.Unmatched);

            props = fieldnames(data);

            if isempty(props) && parser.Results.onlyDifferentFromDefault
                str = 'default';
            else
                str = '';
                for i = 1:numel(props)
                    prop = props{i};
                    value = data.(prop);

                    this = sprintf('%s%s%s%s%s%s', parser.Results.beforeProp, prop, parser.Results.betweenPropValue, ...
                        p.serializePropertyValue(prop, value), parser.Results.afterValue, parser.Results.betweenProps);
                    str = cat(2, str, this);
                end
                if ~isempty(str)
                    str = str(1:end-numel(parser.Results.betweenProps));
                end
            end
        end

        function str = generateShortDifferencesString(p)
            str = p.generateString('onlyDifferentFromDefault', true, ...
                'ignoreHiddenUnlessDifferentFromDefault', true, ...
                'defaultsFromClassDefinition', true, ... % can change this to false if many options are changed in the constructor. They'll still factor into the hash though.
                'beforeProp', '', 'betweenPropValue', '=', 'afterValue', '', 'betweenProps', ' ');
        end

        function str = generateSummaryText(p, indent, paramIndex)
            if nargin < 2
                indent = 0;
            end
            className = class(p);
            if nargin > 2
                indexStr = sprintf('[%d %s %s]', paramIndex, p.paramHashString, p.dataHashString);
            else
                indexStr = sprintf('[%s %s]', p.paramHashString, p.dataHashString);
            end
            header = sprintf('%s%s %s "%s"\n%s%s\n\n', blanks(indent), indexStr, className, p.name, blanks(indent+2), p.generateShortDifferencesString());
            text = p.generateString('onlyDifferentFromDefault', false, ...
                'ignoreHiddenUnlessDifferentFromDefault', true, ...
                'beforeProp', blanks(indent+2), 'betweenPropValue', ': ', 'afterValue', '', 'betweenProps', sprintf('\n')); %#ok<SPRINTFN>
            str = cat(2, header, text, sprintf('\n')); %#ok<SPRINTFN>
        end

        function valstr = serializePropertyValue(p, prop, value) %#ok<INUSL>
            % Generates a string representation like name_value
            % Strips "c\_" from beginnning of name
            %
            % Args:
            %  prop : string
            %    which property
            %  value : Matlab built in type

            threshDots = 6;
            switch class(value)
                case {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64'}
                    if isscalar(value)
                        valstr  = sprintf('%i', value);
                    elseif numel(value) < threshDots
                        valstr = mat2str(value);
                    else
                        valstr = '[...]';
                    end
                    
                case {'logical'}
                    if isscalar(value)
                        if value
                            valstr = 'true';
                        else
                            valstr = 'false';
                        end
                    elseif numel(value) < threshDots
                        valstr = mat2str(value);
                    else
                        valstr = '[...]';
                    end
                    
                case {'double','single'}
                    if isscalar(value)
                        valstr = sprintf('%g', value);
                    elseif numel(value) < threshDots
                        valstr = mat2str(value);
                    else
                        valstr = '[...]';
                    end
                    
                case {'char'}
                    valstr = value;
                    
                case {'cell'}
                    assert(isempty(value) || isvector(value), 'Cell value for property %s must be vector');
                    parts = cell(numel(value), 1);
                    for i = 1:numel(value)
                        parts{i} = p.serializePropertyValue('', value{i});
                    end
                    valstr = sprintf('{%s}', strjoin(parts, ', '));

                otherwise
                    error(['don''t know this type: ' class(value)]);
            end
        end
    end

    methods
        function str = generateCommandLineOptionsString(p, run, varargin)
            % str = generateCommandLineOptionsString(p)
            % Generates a string of all the command line options to
            % pass directly into run_lfads.py. Also takes care of scaling
            % c_kl_increase_steps and c_l2_increase_steps by the number of
            % datasets in a Run when scaleIncreaseStepsWithDatasets == true, 
            % and the same for scaleStartStepWithDatasets
            %
            % Args:
            %   run (LFADS.Run) : run instance for which these options will
            %   be generated
            %
            % Returns:
            %   str : char

            parser = inputParser();
            parser.addParameter('omitFields', {}, @iscellstr);
            parser.parse(varargin{:});

            values = p.getPropertyValueSubset('ignoreProperties', parser.Results.omitFields, ...
                'ignoreHidden', true);

            % get all the parameters except the fields which are skipped
            f = fieldnames(values);

            % only keep the command line params
            %  these begin with 'c_'
            keepFields = false(numel(f), 1);
            for nf=1:numel(f)
                if numel(f{nf}) > 2 && strcmp(f{nf}(1:2), 'c_')
                    keepFields(nf) = true;
                end
            end
            f = f(keepFields);
            
            propsScaleWithDatasets = p.getListPropertiesScaleWithNumDatasets();

            % build the output string
            str = '';
            for nf = 1:numel(f)
                clear fieldstr
                thisField = f{nf};
                thisVal = p.(thisField);
                
                % skip empty values
                if isempty(thisVal)
                    continue;
                end

                % modify specific param values here based on other
                % properties
                if ismember(thisField, propsScaleWithDatasets)
                    % scale thisVal by nDatasets
                    thisVal = thisVal * run.nDatasets;
                end

                % argument formatted differently for each class
                switch class(thisVal)
                    case {'uint16', 'uint32', 'int16', 'int32'}
                        fieldstr = sprintf('%i', thisVal);
                    case {'logical'}
                        if thisVal
                            fieldstr = 'true';
                        else
                            fieldstr = 'false';
                        end
                    case {'double','single'}
                        fieldstr = sprintf('%f', thisVal);
                    case {'char'}
                        fieldstr = thisVal;
                    otherwise
                        error(['don''t know this type: ' class(thisVal)]);
                end
                str = sprintf('%s --%s=%s', str, thisField(3:end), ...
                    fieldstr);
            end

        end

        function str = get.paramHash(p)
            str = p.generateHash();
        end

        function str = get.dataHash(p)
            str = p.generateInputDataHash();
        end

        function str = get.paramHashString(p)
            str = ['param_' p.generateHash()];
        end

        function str = get.dataHashString(p)
            str = ['data_' p.generateInputDataHash()];
        end
    end

    methods(Hidden)
        function p = copy(p)
            % RunParams was once a handle class, this enables copy to go
            % through as intended
        end

        function h = getFirstLineHeader(p)
            className = class(p);
            if isempty(p.name)
                nameStr = '';
            else
                nameStr = sprintf(' "%s"', p.name);
            end
            h = sprintf('%s param_%s data_%s%s\n%s', className, p.generateHash, p.generateInputDataHash, nameStr, p.generateShortDifferencesString());
        end
    end

    methods (Access = protected)
       function header = getHeader(p)
          if ~isscalar(p)
             header = getHeader@matlab.mixin.CustomDisplay(p);
          else
             header = sprintf('%s\n', p.getFirstLineHeader());
          end
       end
    end
    
    methods(Static)
        function diff = diffTwoParams(p1, p2)
            d1 = p1.generateStructForHash();
            d2 = p2.generateStructForHash();
            
            
        end
    end
end
