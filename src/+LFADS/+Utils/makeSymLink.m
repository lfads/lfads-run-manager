function success = makeSymLink(src, link, expandSource)
% makeSymLink(src, linkDest, expandPaths)

    if nargin < 3
        expandSource = true;
    end

    if expandSource % leave false for relative symlink
        src = LFADS.Utils.GetFullPath(src);
    end
        
    link = LFADS.Utils.GetFullPath(link);
    LFADS.Utils.mkdirRecursive(fileparts(link));

    cmd = sprintf('ln -sfn "%s" "%s"', src, link);
    [status, output] = unix(cmd);
    
    if status
        fprintf('Error creating symlink: \n');
        fprintf('%s\n', output);
    end

    success = ~status;
end
