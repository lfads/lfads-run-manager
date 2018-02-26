%% This script walks through running LFADS on a single Lorenz dataset

%% Generate synthetic Lorenz datasets

% build the dataset collection
datasetPath = '~/lorenz_example/datasets';

% generate demo datasets
if ~exist(fullfile(datasetPath, 'dataset001.mat'), 'file')
    LFADS.Utils.generateDemoDatasets(datasetPath, 'nDatasets', 3);
end

%% Locate and specify the datasets
dc = MyExperiment.DatasetCollection(datasetPath);
dc.name = 'lorenz_example';

% add individual datasets
MyExperiment.Dataset(dc, 'dataset001.mat');
MyExperiment.Dataset(dc, 'dataset002.mat');
MyExperiment.Dataset(dc, 'dataset003.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

% print information loaded from each dataset
dc.getDatasetInfoTable()

%% Build RunCollection - first a single run

% Run a single model for each of the datasets
runRoot = '~/lorenz_example/runs';
rc_ds1 = MyExperiment.RunCollection(runRoot, 'exampleSingleRun', dc);

% replace this with the date this script was authored as YYYYMMDD 
% This ensures that updates to lfads-run-manager won't invalidate older 
% runs already on disk and provides for backwards compatibility
rc_ds1.version = 20171107;

%% Set some hyperparameters

par = MyExperiment.RunParams;
par.spikeBinMs = 2; % rebin the data at 2 ms
par.c_co_dim = 0; % no controller --> no inputs to generator
par.c_batch_size = 150; % must be < 1/5 of the min trial count
par.c_factors_dim = 8; % and manually set it for multisession stitched models
par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN
par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc_ds1.addParams(par);

%% Create the RunSpecs

% Define a RunSpec, which indicates which datasets get included, as well as
% what name to call this run spec
runSpecName = dc.datasets(1).getSingleRunName(); % generates a simple run name from this datasets name
whichDataset = 1;
runSpec = MyExperiment.RunSpec(runSpecName, dc, 1);

% add this RunSpec to the RunCollection
rc_ds1.addRunSpec(runSpec);

% adding a return here allows you to call this script to recreate all of
% the objects here for subsequent analysis after the actual LFADS models
% have been trained. The code below will setup the LFADS training runs on 
% disk the first time around, and should be run once manually.
return;

%% Prepare LFADS input and shell scripts

% generate all of the data files LFADS needs to run everything
rc_ds1.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
% you should set display to a valid x display
% Other options are available 
rc_ds1.writeShellScriptRunQueue('display', 0, 'virtualenv', 'tensorflow');



