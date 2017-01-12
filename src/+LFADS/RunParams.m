classdef RunParams < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
    % Collection of parameters which are common to all :ref:`LFADS_Run` instances in a :ref:`LFADS_RunCollection`. You
    % must create a subclass of RunParams in which you specify the serialized representation of the parameters that
    % will be used in paths on disk.
    
    properties
        % High-level params unrelated to LFADS internals
        
        spikeBinMs double = 2; % Spike bin width in ms
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
            
            parser = inputParser();
            parser.addParameter('ignoreProperties', {}, @iscellstr);
            parser.addParameter('onlyRootClassProperties', false, @islogical);
            parser.parse(varargin{:});
            
            if parser.Results.onlyRootClassProperties
                meta = ?LFADS.RunParams;
            else
                meta = metaclass(p);
            end
            
            props = cell(numel(meta.PropertyList), 1);
            mask = false(numel(meta.PropertyList), 1);
            for i = 1:numel(meta.PropertyList)
                prop = meta.PropertyList(i);
                name = prop.Name;
                props{i} = name;
                if ismember(name, parser.Results.ignoreProperties)
                    continue;
                end
                % skip properties that are Dependent, Constant, Transient,
                % or Hidden. Serialize the value if it differes from the 
                if ~prop.Dependent && ~prop.Constant && ~prop.Transient && ~prop.Hidden
                    mask(i) = true;
                end
            end
            
            props = props(mask);
            propMeta = meta.PropertyList(mask);
        end

        function hash = generateHash(p, varargin)
            % Generate a short hash of this RunParams non-transient
            % properties that can be used in a directory name.
            %
            % Args:
            %   length : int
            %     number of characters to truncate the hash value to
            %   ignoreProperties : cellstr
            %     list of properties to ignore
            %   onlyRootClassProperties : bool (False)
            %     if true, only include properties declared in
            %     LFADS.RunParams. if false, include properties declared in
            %     subclasses.
            parser = inputParser();
            parser.addParameter('length', 6, @isscalar);
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            
            props = p.listNonTransientProperties(parser.Unmatched);
            data = struct();
            for i = 1:numel(props)
                prop = props{i};
                data.(prop) = p.(prop);
            end
            
            hash = LFADS.Utils.DataHash(data, struct('Format', 'base64'));
            hash = lower(hash(1:parser.Results.length));
        end
        
        function str = generateHashName(p, varargin)
            % Generate a short hash name like 'param_HASH' where HASH is
            % generated by generateHash()
            
            str = sprintf('param_%s', p.generateHash(varargin{:}));
        end
            
        
        function str = generateString(p, varargin)
            % Generates a string representation of all parameters, with custom
            % strings inserted between property names and values. 
            %
            % Args:
            %   filterDiffersFromDefault : bool
            %     if true, only properties whose value differs from the
            %     default value 
            %   defaultsFromClassDefinition : bool
            %     if true, default value comes from the class definition,
            %     next to each property's definition. If false, the default
            %     values come from constructing a new class of the same
            %     type with no arguments. Default true.
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
            % Returns:
            %   str : string
            %     serialized string suffix suitable for inclusion in file paths
            %
            
            parser = inputParser();
            parser.addParameter('onlyDifferentFromDefault', false, @islogical);
            parser.addParameter('defaultsFromClassDefinition', true, @islogical);
            parser.addParameter('beforeProp', '', @ischar);
            parser.addParameter('betweenPropValue', '=', @ischar);
            parser.addParameter('afterValue', '', @ischar);
            parser.addParameter('betweenProps', ' ', @ischar);
            parser.parse(varargin{:});
            
            [props, propMeta] = p.listNonTransientProperties(parser.Unmatched);
            str = '';
            
            defaultInstance = eval(class(p));
                
            for i = 1:numel(props)
                prop = props{i};
                value = p.(prop);
                
                if parser.Results.onlyDifferentFromDefault
                    if parser.Results.defaultsFromClassDefinition
                        def = propMeta(i).DefaultValue;
                    else
                        def = defaultInstance.(prop);
                    end
                    if isequal(def, value)
                        continue;
                    end
                end
                
                this = sprintf('%s%s%s%s%s%s', parser.Results.beforeProp, prop, parser.Results.betweenPropValue, ...
                    p.serializePropertyValue(prop, value), parser.Results.afterValue, parser.Results.betweenProps);
                str = cat(2, str, this);
            end
            if ~isempty(str)
                str = str(1:end-numel(parser.Results.betweenProps));
            end
            if parser.Results.onlyDifferentFromDefault && isempty(str)
                str = 'default';
            end
        end
        
        function str = generateShortDifferencesString(p)
            str = p.generateString('onlyDifferentFromDefault', true, ...
                'defaultsFromClassDefinition', false, ...
                'beforeProp', '', 'betweenPropValue', '=', 'afterValue', '', 'betweenProps', ' ');
        end
        
        function str = generateSummaryText(p, indent, paramIndex)
            if nargin < 2
                indent = 0;
            end
            className = class(p);
            if nargin > 2
                indexStr = sprintf('[%d param_%s]', paramIndex, p.generateHash);
            else
                indexStr = sprintf('[param_%s]', p.generateHash);
            end
            header = sprintf('%s%s %s\n%sDiff: %s\n\n', blanks(indent), indexStr,className, blanks(indent+2), p.generateShortDifferencesString());
            text = p.generateString('onlyDifferentFromDefault', false, ...
                'beforeProp', blanks(indent+2), 'betweenPropValue', ': ', 'afterValue', '', 'betweenProps', sprintf('\n'));
            str = cat(2, header, text, sprintf('\n'));
        end
        
        function valstr = serializePropertyValue(p, prop, value)
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
            h = sprintf('%s %s', className, p.generateShortDifferencesString());
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
