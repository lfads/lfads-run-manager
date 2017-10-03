%% Generate synthetic Lorenz datasets

% build the dataset collection
datasetPath = '~/lorenz_example/datasets';
nDatasets = 4;

% generate demo datasets
if ~exist(fullfile(datasetPath, sprintf('dataset%03d.mat', nDatasets)), 'file')
    LFADS.Utils.generateDemoDatasets(datasetPath, nDatasets);
end

%% Locate and specify the datasets
dc = MyExperiment.DatasetCollection(datasetPath);

% add individual datasets
MyExperiment.Dataset(dc, 'dataset001.mat');
MyExperiment.Dataset(dc, 'dataset002.mat');
MyExperiment.Dataset(dc, 'dataset003.mat');
MyExperiment.Dataset(dc, 'dataset004.mat');

% load metadata from the datasets to populate the dataset collection
dc.loadInfo;

%% Set parameters for the entire run collection

par = MyExperiment.RunParams;
par.useAlignmentMatrix = true;
par.scaleIncreaseStepsWithDatasets = true;
par.c_in_factors_dim = 8;
par.c_factors_dim = 8;

%% Build RunCollection
% Run a single model for each dataset, and one stitched run with all datasets

runRoot = '~/lorenz_example/runs';
rc = MyExperiment.RunCollection(runRoot, 'exampleRun', dc);
 
% add a single set of parameters to this run collection. Additional
% parameters can be added. LFADS.RunParams is a value class, unlike the other objects
% which are handle classes, so you can modify par freely.
rc.addParams(par);

% add each individual run
for iR = 1:dc.nDatasets
    runSpec = MyExperiment.RunSpec(dc.datasets(iR).getSingleRunName(), dc, dc.datasets(iR).name);
    rc.addRunSpec(runSpec);
end

% add the final stitching run with all datasets
rc.addRunSpec(MyExperiment.RunSpec('all', dc, 1:dc.nDatasets));

% adding a return here allows you to call this script to recreate all of
% the objects here for subsequent analysis after the actual LFADS models
% have been trained. The code below will setup the LFADS runs in the first
% place.

return; 

%% Prepare LFADS input

% generate all of the data files LFADS needs to run everything
rc.prepareForLFADS();

% write a text file summarizing the run specs and parameters used
rc.writeSummaryText();

% write a python script that will train all of the LFADS runs using a
% load-balancer against the available CPUs and GPUs
rc.writeShellScriptRunQueue('display', 42, 'maxTasksSimultaneously', 20, 'gpuList', [0 1 2 3 6 7 8]);
