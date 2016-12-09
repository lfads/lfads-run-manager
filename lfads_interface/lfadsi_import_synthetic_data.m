function out = lfads_import_synthetic_data(infile)
%% imports results from lfads
%%
%% inputs: 
%%
%%   infile: filename with stored data


    in_info = h5info(infile);
    for nn = 1:numel(in_info.Datasets)
        d = in_info.Datasets(nn).Name;
        tmp = h5read(infile, sprintf('/%s', d));
        out.(d) = tmp;
    end
