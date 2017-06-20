function success = makeSymLink(src, link)
% makeSymLink(src, linkDest)

    src = LFADS.Utils.GetFullPath(src);
    link = LFADS.Utils.GetFullPath(link);
    LFADS.Utils.mkdirRecursive(fileparts(link));
    if exist(link, 'file')
        delete(link);
    end
    cmd = sprintf('ln -s "%s" "%s"', src, link);
    [status, output] = unix(cmd);
    
    if status
        fprintf('Error creating symlink: \n');
        fprintf('%s\n', output);
    end

    success = ~status;
end
