%% This script walks through running LFADS to stitch multiple Lorenz datasets

%% Generate synthetic Lorenz datasets

% build the dataset collection
datasetPath = '~/lorenz_example/datasets';

% generate demo datasets
if ~exist(fullfile(datasetPath, 'dataset001.mat'), 'file')
    LFADS.Utils.generateDemoDatasets(datasetPath, 'nDatasets', 3);
end

%% Locate and specify the datasets
dc = LorenzExperiment.DatasetCollection(datasetPath);
dc.name = 'lorenz_example';

% add individual datasets
LorenzExperiment.Dataset(dc, 'dataset001.mat');
LorenzExperiment.Dataset(dc, 'dataset002.mat');
LorenzExperiment.Dataset(dc, 'dataset003.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

% print information loaded from each dataset
dc.getDatasetInfoTable()

%% Set some hyperparameters

par = LorenzExperiment.RunParams;
par.spikeBinMs = 2; % rebin the data at 2 ms
par.c_co_dim = 0; % no controller --> no inputs to generator
par.c_batch_size = 150; % must be < 1/5 of the min trial count
par.c_factors_dim = 8; % and manually set it for multisession stitched models
par.useAlignmentMatrix = true; % use alignment matrices initial guess for multisession stitching

par.c_gen_dim = 64; % number of units in generator RNN
par.c_ic_enc_dim = 64; % number of units in encoder RNN

par.c_learning_rate_stop = 1e-3; % we can stop really early for the demo

%% Running a multi-dataset stitching run

% Now we'll do something slightly more interesting. We'll do a total of 4
% LFADS runs. 3 will be single-session LFADS runs on each of the datasets
% individually. The last will be a multi-session stitched dataset that
% leverages all 3 of the datasets in a common shared model.

runRoot = '~/lorenz_example/runs';
rc = LorenzExperiment.RunCollection(runRoot, 'exampleStitching', dc);

% replace this with the date this script was authored as YYYYMMDD
% This ensures that updates to lfads-run-manager won't invalidate older
% runs already on disk and provides for backwards compatibility
rc.version = 20171107;

% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);

% Add a RunSpec to train individual models for each dataset
for iR = 1:dc.nDatasets
    runSpec = LorenzExperiment.RunSpec(dc.datasets(iR).getSingleRunName(), dc, iR);
    rc.addRunSpec(runSpec);
end

% Add a RunSpec using all datasets which LFADS will then "stitch" into a
% shared dynamical model
rc.addRunSpec(LorenzExperiment.RunSpec('all', dc, 1:dc.nDatasets));

% adding a return here allows you to call this script to recreate all of
% the objects here for subsequent analysis after the actual LFADS models
% have been trained. The code below will setup the LFADS training runs on
% disk the first time around, and should be run once manually.
return;

%% Prepare LFADS input and shell scripts

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
% you should set display to a valid x display
% Other options are available
rc.writeShellScriptRunQueue('display', 0, 'virtualenv', 'tensorflow');

%% Looking at the alignment matrices used

runStitched = rc.findRuns('all', 1); % 'all' looks up the RunSpec by name, 1 refers to the first (and here, the only) RunParams

alignTool = runStitched.multisessionAlignmentTool;
if isempty(alignTool)
    runStitched.doMultisessionAlignment();
    alignTool = runStitched.multisessionAlignmentTool;
end

alignTool.plotAlignmentReconstruction();

% You want the colored traces to resemble the black "global" trace. The
% black traces are the PC scores using data from all the datasets. The
% colored traces are the best linear reconstruction of the black traces
% from each individual dataset alone. The projection which achieves this
% best reconstruction is used as the initial seed for the readin matrices
% for LFADS, which can be trainable or fixed depending on
% par.do_train_readin.

%% Run LFADS

% You should now run at the command line
% source activate tensorflow   # if you're using a virtual machine
% python ~/lorenz_example/runs/exampleRun_dataset1/run_lfadsqueue.py

% And then wait until training and posterior sampling are finished
