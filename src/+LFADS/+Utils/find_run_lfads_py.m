function folder = find_run_lfads_py(warn)
    if nargin < 1
        warn = false;
    end
    [~, s] = system('which run_lfads.py');

    if contains(s, 'not found')
        if warn
            warning('run_lfads.py was not found on the system PATH');
        end
        folder = '';
    else
        folder = fileparts(strip(s));
    end

end
