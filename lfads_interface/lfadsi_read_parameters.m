function out = lfadsi_read_parameters(lfdir)

f = fullfile(lfdir,'hyperparameters-0.txt');
out = loadjson(f);
