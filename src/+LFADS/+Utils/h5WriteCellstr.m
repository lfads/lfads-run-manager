function h5WriteCellstr(filename, name, cstr)
% based entirely on https://www.mathworks.com/matlabcentral/fileexchange/24091-hdf5-read-write-cellstr-example

% Generate a file
if ~exist(filename, 'file')
    fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
else
    fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
end

% Set variable length string type
VLstr_type = H5T.copy('H5T_C_S1');
H5T.set_size(VLstr_type,'H5T_VARIABLE');

% Create a dataspace for cellstr
H5S_UNLIMITED = H5ML.get_constant_value('H5S_UNLIMITED');
dspace = H5S.create_simple(1,numel(cstr),H5S_UNLIMITED);

% Create a dataset plist for chunking
plist = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(plist, 2); % 2 strings per chunk

% Create dataset
dset = H5D.create(fid, name, VLstr_type,dspace,plist);

% Write data
H5D.write(dset,VLstr_type,'H5S_ALL','H5S_ALL','H5P_DEFAULT', cstr);

% Close file & resources
H5P.close(plist);
H5T.close(VLstr_type);
H5S.close(dspace);
H5D.close(dset);
H5F.close(fid);

end