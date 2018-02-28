classdef RunParams < matlab.mixin.CustomDisplay 
    % Collection of parameters which are common to all :ref:`LFADS_Run` instances in a :ref:`LFADS_RunCollection`. You
    % must create a subclass of RunParams in which you specify the serialized representation of the parameters that
    % will be used in paths on disk. 
    % 
    % When you add additional parameters to this list, don't prefix yours
    % with c_, as these have special meaning to the LFADS code and will
    % mess things up.
    
    properties
        % High-level params unrelated to LFADS internals
        spikeBinMs double = 2; % Spike bin width in ms
        trainToTestRatio uint16 = 4; % how many train v. test trials, defaults to 4:1 ratio    
        useAlignmentMatrix logical = false; % Whether to use an alignment matrix when stitching datasets together.
        useSingleDatasetAlignmentMatrix logical = false;  % Whether to use an alignment matrix using a single dataset, for dimensionality reduction upstream of the encoder 
        scaleIncreaseStepsWithDatasets logical = true; % If true, c_kl_increase_steps and c_l2_increase_steps will be multiplied by the number of datasets in a Run
    end
    
    properties
        % Subset of the command-line params to run_lfads.py
        % Each starts with c_ so it can be differentiated from other params
        % Avoid naming other parameters the same as these without the c_
        % since the prefix will be stripped from the serialization
        
        c_cell_clip_value double = 5; % used to avoid stepping too far during training
        c_factors_dim uint16 = 50;
        c_ic_enc_dim uint16 = 128; % network size for IC encoder
        c_ci_enc_dim uint16 = 128; % network size for controller input encoder
        c_gen_dim uint16 = 100; % generator network size
        c_keep_prob double = 0.95; % randomly drop units during each training pass
        c_learning_rate_decay_factor double = 0.98; % how quickly to decrease the learning rate
        c_device char = '/gpu:0'; % which visible GPU/CPU to use
        
        c_co_dim uint16 = 4;
        c_do_causal_controller logical = false; % restrict input encoder from seeing the future?DO_FEED_FACTORS_TO_CONTROLLER
        c_do_feed_factors_to_controller logical = true; % restrict input encoder from seeing the future?
        c_feedback_factors_or_rates char = 'factors'; % Feedback the factors or the rates to the controller? Acceptable values: 'factors' or 'rates'
        c_controller_input_lag uint16 = 1;
        
        c_do_train_readin logical = true; % for stitching models, make the readin matrices trainable (true) or fix them to equal the alignment matrices (false)
        
        c_l2_gen_scale double = 500; % how much to weight the generator l2 cost
        c_l2_con_scale double = 500; % how much to weight the controller l2 cost
        c_batch_size uint16 = 256; % number of trials to use during each training pass
        c_kl_increase_steps uint16 = 900; % Number of steps over which the kl costs increase
        c_l2_increase_steps uint16 = 900; % Number of steps over which the l2 costs increase
        c_ic_dim uint16 = 64; % dimensionality of the initial conditions
        c_con_dim uint16 = 128; %controller dimensionality
        
        c_learning_rate_stop = 0.00001; % when the learning rate reaches this threshold, stop training
        c_temporal_spike_jitter_width uint16 = 0; % jittering spike times during training, in units of bin size
        
        c_allow_gpu_growth logical = true; %whether to allow the GPU to dynamically allocate memory. default (false) is to allocate all the memory initially
        c_kl_ic_weight double = 1; % how much to weight the generator l2 cost
        c_kl_co_weight double = 1; % how much to weight the controller l2 cost
        
        c_inject_ext_input_to_gen logical = false; % should observed inputs be input to model via encoders (false) or injected directly into generator (true)?
    end
    
    properties
        % These properties DO NOT affect the param_HASH, as they are
        % returned by getListPropertiesNotAffectingParamHash below.
        % They also do not affect data_HASH as they do not begin with c_
        
        num_samples_posterior = 512; % number of samples 
        version uint32 = 20171107; % Used for graceful evolution of path settings
        name char = ''; % convenient name for displaying this param
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
            
            props = cell(numel(meta.PropertyList), 1);
            mask = false(numel(meta.PropertyList), 1);
            for i = 1:numel(meta.PropertyList)
                prop = meta.PropertyList(i);
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
            propMeta = meta.PropertyList(mask);
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
                        def = propMeta(i).DefaultValue;
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
            list = {'version', 'num_samples_posterior'};
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
            
            length = 6;
            ignore = p.getListPropertiesNotAffectingHash();
            data = p.getPropertyValueSubset('ignoreProperties', ignore, 'onlyDifferentFromDefault', true);
            hash = LFADS.Utils.DataHash(data, struct('Format', 'base64'));
            hash = strrep(strrep(hash, '/', '_'), '+', '-'); % https://tools.ietf.org/html/rfc3548#page-6
            if numel(hash) > length
                hash = hash(1:length);
            end
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
            propsIgnore = p.getListPropertiesNotAffectingInputDataHash();
            data = p.getPropertyValueSubset('ignoreProperties', propsIgnore, 'onlyDifferentFromDefault', true);
            hash = LFADS.Utils.DataHash(data, struct('Format', 'base64'));
            hash = strrep(strrep(hash, '/', '_'), '+', '-'); % https://tools.ietf.org/html/rfc3548#page-6
            if numel(hash) > length
                hash = hash(1:length); 
            end
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
            
            switch class(value)
                case {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64'}
                    valstr  = sprintf('%i', value);
                case {'logical'}
                    if value
                        valstr = 'true';
                    else
                        valstr = 'false';
                    end
                case {'double','single'}
                    valstr = sprintf('%g', value);
                case {'char'}
                    valstr = value;
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
            % datasets in a Run when scaleIncreaseStepsWithDatasets == true
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
            
            % build the output string
            str = '';
            for nf = 1:numel(f)
                clear fieldstr
                thisField = f{nf};
                thisVal = p.(thisField);
                
                % modify specific param values here based on other
                % properties
                if ismember(thisField, {'c_kl_increase_steps', 'c_l2_increase_steps'}) && ...
                        p.scaleIncreaseStepsWithDatasets
                    % scale thisVal by nDatasets
                    thisVal = thisVal * run.nDatasets;
                end
                
%                 if ismember(thisField, {'c_in_factors_dim'}) && ...
%                         p.setInFactorsMatchDataForSingleDataset && run.nDatasets == 1
%                     thisVal = zeros(1, 'like', thisVal); % to auto scale
%                 end
                
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
    
end
