function mkdirRecursive(dirPath, cdTo)
% like mkdir -p : creates intermediate directories as required

s = warning('off', 'MATLAB:MKDIR:DirectoryExists');

if exist(dirPath, 'dir')
    if nargin >= 2 && cdTo
        cd(dirPath);
    end
    return;
else
    parent = fileparts(dirPath);
    if ~isempty(parent)
        mkdirRecursive(parent);
    end

    mkdir(dirPath);
end

if nargin < 2
    cdTo = false;
end
if cdTo
    cd(dirPath);
end
warning(s);

