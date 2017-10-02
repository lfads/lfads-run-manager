function package = getPackage(varargin)

    p = inputParser;
    p.addParameter('stackOffset', 0, @isscalar);
    p.addParameter('up', 0, @isscalar); % how many parents to go up to
    p.parse(varargin{:});

    files = dbstack('-completenames');
    if numel(files) == 2
        % must be executing from cell mode
        path = fileparts(matlab.desktop.editor.getActiveFilename());
    else
        path = fileparts(files(2 + p.Results.stackOffset).file);
    end
    package = '';

    upRemaining = p.Results.up;

    while ~isempty(path)
        [path, name] = fileparts(path);
        if upRemaining > 0
            upRemaining = upRemaining - 1;
            continue;
        end
        if strcmp(name(1), '+')
            if isempty(package)
                package = name(2:end);
            else
                package = sprintf('%s.%s', name(2:end), package);
            end
        else
            break;
        end
    end

end