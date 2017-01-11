classdef RunParams < handle & matlab.mixin.CustomDisplay
    % Collection of parameters which are common to all :ref:`LFADS_Run` instances in a :ref:`LFADS_RunCollection`. You
    % must create a subclass of RunParams in which you specify the serialized representation of the parameters that
    % will be used in paths on disk.
    
    properties
        % High-level params unrelated to LFADS internals
        
        spikeBinMs uint16 = 2; % Spike bin width in ms
        trainToTestRatio uint16 = 4; % how many train v. test trials, defaults to 4:1 ratio    
        useAlignmentMatrix logical = false; % Whether to use an alignment matrix when stitching datasets together. Default = false.
    end
    
    properties
        % Subset of the command-line params to run_lfads.py
        % Each starts with c_ so it can be differentiated from other params
        % Avoid naming other parameters the same as these without the c_
        % since the prefix will be stripped from the serialization
        
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
        c_do_causal_controller logical = false; % restrict input encoder from seeing the future?
        c_l2_gen_scale double = 500; % how much to weight the generator l2 cost
        c_l2_con_scale double = 500; % how much to weight the controller l2 cost
        c_batch_size uint16 = 256; % number of trials to use during each training pass
        c_kl_increase_steps uint16 = 900; % Number of steps over which the kl costs increase
        c_l2_increase_steps uint16 = 900; % Number of steps over which the l2 costs increase
        c_controller_input_lag uint16 = 1;
        c_ic_dim uint16 = 64; % dimensionality of the initial conditions
        c_con_dim uint16 = 128; %controller dimensionality
    end
    
    methods
        % Methods that you may wish to override in custom subclasses
        
        function str = generateString(p, varargin)
            % Generates a string representation of the parameters that is used for reading off the parameters from a
            % folder path on disk. You can include as many or as few fields in this serialized representation as you see
            % fit. The default implementation will compare all property values
            % (including those defined in subclasses) to the inline default property
            % values defined in the class that defines them. Any property
            % values that differ will be included in the serialization.
            %
            % By default this will be a list of properties and values
            % separated by double underscores, e.g.
            % 
            % .. code::
            %
            %     prop1Name_value1__prop2Name_value2__prop3Name_value3
            %
            % Args:
            %   ignoreProperties : cellstr
            %     list of properties to ignore when serializing, useful
            %     when called from subclasses
            %   onlyRootClassProperties : bool (False)
            %     if true, only include properties declared in
            %     LFADS.RunParams. if false, include properties declared in
            %     subclasses.
            %
            % Returns:
            %   str : string
            %     serialized string suffix suitable for inclusion in file paths
            %
            
            parser = inputParser();
            parser.addParameter('ignoreProperties', {}, @iscellstr);
            parser.addParameter('onlyRootClassProperties', false, @islogical);
            parser.parse(varargin{:});
            
            if parser.Results.onlyRootClassProperties
                meta = ?LFADS.RunParams;
            else
                meta = metaclass(p);
            end
            
            str = '';
            for i = 1:numel(meta.PropertyList)
                prop = meta.PropertyList(i);
                name = prop.Name;
                if ismember(name, parser.Results.ignoreProperties)
                    continue;
                end
                % skip properties that are Dependent, Constant, Transient,
                % or Hidden. Serialize the value if it differes from the 
                if ~prop.Dependent && ~prop.Constant && ~prop.Transient && ~prop.Hidden
                    if ~isequal(prop.DefaultValue, p.(name))
                        str = cat(2, str, p.serializeProperty(name, p.(name)), '__');
                    end
                end
            end
            if isempty(str)
                str = 'default';
            end
            % strip __ off the end
            if strcmp(str(end-1:end), '__'), str = str(1:end-2); end
        end
        
        function str = serializeProperty(p, name, value)
            % Generates a string representation like name_value
            % Strips "c\_" from beginnning of name
            %
            % Args:
            %  name : string
            %  value : Matlab built in type
            
            switch class(value)
                case {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64'}
                    valstr  = sprintf('%i', value);
                case {'logical'}
                    if thisVal
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
            if strncmp(name, 'c_', 2)
                name = name(3:end);
            end
            str = sprintf('%s_%s', name, valstr);
        end
    end
    
    methods
        function str = generateCommandLineOptionsString(p)
            % str = generateCommandLineOptionsString(p)
            % Generates a string of all the command line options to
            % pass directly into run_lfads.py
            %
            % Returns:
            %   str : char
            
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
                        if thisVal
                            fieldstr = 'true';
                        else
                            fieldstr = 'false';
                        end
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
    
    
    methods(Hidden)
        function h = getFirstLineHeader(p)
            className = class(p);
            h = sprintf('%s "%s"', className, p.generateString());
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
