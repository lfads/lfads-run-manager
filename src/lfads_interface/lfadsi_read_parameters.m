function out = lfadsi_read_parameters(lfdir)

f = fullfile(lfdir,'hyperparameters-0.txt');
s = warning('off', 'MATLAB:namelengthmaxexceeded');
out = loadjson(f);
warning(s);
