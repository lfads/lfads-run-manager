function [status, result] = chmod(permissions, fileList, varargin)
    % runs chmod permissions fileList
    p = inputParser;
    p.addRequired('permissions', @(x) isempty(x) || ischar(x));
    p.addRequired('fileList', @(x) ischar(x) || iscellstr(x));
    p.addParamValue('recursive', false, @islogical);
    p.addParamValue('printError', true, @islogical);
    p.parse(permissions, fileList, varargin{:});
    printError = p.Results.printError;
    recursive = p.Results.recursive;

    if isempty(permissions)
        status = 0;
        result = '';
        return;
    end
    
    if ~iscell(fileList)
        fileList = {fileList};
    end

    fileList = cellfun(@LFADS.Utils.GetFullPath, fileList, 'UniformOutput', false);
    fileListEscaped = cellfun(@(path) strrep(path, ' ', '\ '), fileList, 'UniformOutput', false);
    fileListString = strjoin(fileListEscaped, ' ');

    if recursive
        flags = '-R ';
    else
        flags = '';
    end
    cmd = sprintf('chmod %s%s %s', flags, permissions, fileListString);

    [status, result] = system(cmd);
    
    if status && printError
        tcprintf('light red', 'Error running chmod:\n');
        tcprintf('light red', '%s\n', result);
    end

end
