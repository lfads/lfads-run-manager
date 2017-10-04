function out = read_parameters(lfdir)

f = fullfile(lfdir,'hyperparameters-0.txt');
s = warning('off', 'MATLAB:namelengthmaxexceeded');
out = loadjson(LFADS.Utils.GetFullPath(f));
warning(s);
