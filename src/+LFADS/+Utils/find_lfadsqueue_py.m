function folder = find_lfadsqueue_py()
    if nargin < 1
        warn = false;
    end
    
    stack = dbstack('-completenames');
    thisFilePath = fileparts(stack(1).file);
    
    folder = fileparts(fileparts(thisFilePath));

end
