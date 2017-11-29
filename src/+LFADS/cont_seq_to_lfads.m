function cont_seq_to_lfads(seqs, outpath, outfiles, varargin)
% function seq_to_lfads(seqs, varargin)
%   pass in a cell array of 'seq' structs - one element per trial
%   if there are different 'runID' elements, then the struct is split
%   
%
%   seq must contain:
%     y (ms binned spikes, N neurons x T ms)
%   seq may optionally contain:
%     x_true (true values of latents for each trial)
%     y_true (true values of underlying FRs for each trial)
%
%
%   outpath - will be appended to the default output directory
%       e.g. - 'lorenz'
%
%   varargins: trainInds, testInds, binSizeMS, inputBinSizeMS, alignment_matrix_cxf

LFADS.Utils.mkdirRecursive(outpath);

p = inputParser();
p.addParameter('trainInds', {}, @isvector);
p.addParameter('testInds', {}, @isvector);
p.addParameter('binSizeMs', 5, @isscalar);
p.addParameter('inputBinSizeMs', 1, @isscalar);
p.addParameter('conversion_factor', 0.5, @isscalar);
p.addParameter('alignment_matrix_cxf', {}, @isvector);
p.addParameter('alignment_bias_c', {}, @iscell);
p.parse(varargin{:});
trainInds = p.Results.trainInds;
testInds = p.Results.testInds;
binSizeMS = p.Results.binSizeMs;
inputBinSizeMS = p.Results.inputBinSizeMs;
conversion_factor = p.Results.conversion_factor;
alignment_matrix_cxf = p.Results.alignment_matrix_cxf;
alignment_bias_c = p.Results.alignment_bias_c;

% rem = assignopts({'trainInds', 'testInds', 'binSizeMS', 'inputBinSizeMS', ...
%                  'conversion_factor', 'alignment_matrix_cxf', 'whichChannels'}, varargin);
if isempty(trainInds) || isempty(testInds)
    error('specify trainInds and testInds');
end

% turn per-dataset params into cell arrays
if ~iscell(seqs)
    seqs = {seqs};
end
if ~iscell(trainInds)
    trainInds = {trainInds};
end
if ~iscell(testInds)
    testInds = {testInds};
end

if ~exist('whichChannels','var')
    for nset = 1:numel(seqs)
        nNeuronsTmp = unique(arrayfun(@(x) size(x.y,1), seqs{nset}));
        assert(numel(nNeuronsTmp)==1, 'trials do not have consistent #s of neurons');
        whichChannels{nset} = 1:nNeuronsTmp;
    end
end

if ~iscell(whichChannels)
    whichChannels = {whichChannels};
end


% rename the variables, will need the "trainInds" and "testInds"
%    names later...
allTrainInds = trainInds;
allTestInds = testInds;

% % get the default pointers to data / code
% tmp = comparisonDefaults();
% dataPath = tmp.dataPath;
% clear tmp;

% where to save data
% saveDir = fullfile(dataPath, outpath);
spikesDir = outpath;

% spikes that get sent into lfads
% spikesDir = fullfile(saveDir,'lfads','data');
% place where lfads should store its networks
% tfNetworksDir = fullfile(saveDir,'lfads','networks');
% a cmd to call lfads with
% cmdSave = '/home/chethan/tmp/lfadscmd';

% resampling the data


% create any dirs that don't exist
dirsToCreate = {spikesDir};
for nn = 1:numel(dirsToCreate)
    if ~isdir(dirsToCreate{nn})
        mkdir(dirsToCreate{nn});
    end
end

prog = LFADS.Utils.ProgressBar(numel(seqs), 'Writing LFADS Input files');
for ndset = 1:numel(seqs)
    prog.update(ndset);
    
    seq = seqs{ndset};
    trainInds = allTrainInds{ndset};
    testInds = allTestInds{ndset};
    whichChannelsThisSet = whichChannels{ndset};

    outfile = fullfile(outpath, outfiles{ndset});
    
%     if numel(seqs) == 1
%         outfile = fullfile(spikesDir,'lfads_spikes.h5');
%     else
%     end
    
%     delete old spikes
    if exist(outfile,'file')
%         warning(sprintf('Deleting file %s...', outfile));
        delete(outfile);
    end


    %% we are going to collapse all units into a 3-D matrix. first, verify that sizes are constant over all trials
    %nTimeMS = unique(arrayfun(@(x) size(x.y,2), seq)) * ;
    %assert(numel(nTimeMS)==1, ['trials do not have constant ' ...
    %                    'length']);
    % neural dimension
    nNeurons = numel(whichChannelsThisSet);

    % how many output time sample should there be
    nTimeSamples = size(seq(1).y, 2);
    
    % how many input timebins should we keep
    %inputTimeBinsToKeep = nTimeSamples * double( binSizeMS ) / double( inputBinSizeMS );

    nTrainTrials = numel(trainInds);
    nTestTrials = numel(testInds);

    % initialize output variables
    ytrain = zeros(nTrainTrials, nTimeSamples, nNeurons);
    ytest = zeros(nTestTrials, nTimeSamples, nNeurons);
%     if isfield(seq,'x_true')
%         xtrain_true = zeros(nTrainTrials, nTimeSamples, size(seq(1).x_true,1));
%         xtest_true = zeros(nTestTrials, nTimeSamples, size(seq(1).x_true,1));
%     end
    if isfield(seq,'y_true')
        ytrain_true = zeros(nTrainTrials, nTimeSamples, nNeurons);
        ytest_true = zeros(nTestTrials, nTimeSamples, nNeurons);
    end
    
    % compile the train data
    for it = 1:nTrainTrials
        nn = trainInds(it);
        
        ytrain(it,:,:) = seq(nn).y';
        % store down the (optional) true latents and FRs
        if isfield(seq,'x_true')
            xtrain_true(it,:,:) = seq(nn).x_true';
        end
        
        if isfield(seq,'y_true')
            ytrain_true(it,:,:) = seq(nn).y_true';
        end
    end

    % compile the test data
    for it = 1:nTestTrials
        nn = testInds(it);
        ytest(it,:,:) = seq(nn).y';
        % store down the (optional) true latents and FRs
        if isfield(seq,'x_true')
            xtest_true(it,:,:) = seq(nn).x_true';
        end
        if isfield(seq,'y_true')
            ytest_true(it,:,:) = seq(nn).y_true';
        end        
    end

%     fprintf('Each trial is %g bins\n', nTimeSamples);

    varout = {};

    if isfield(seq,'y_true')
        varout{end+1} = 'train_truth';
        varout{end+1} = ytrain_true;
        varout{end+1} = 'valid_truth';
        varout{end+1} = ytest_true;
    end
    
    if numel(varout)
        % if sending in "truth" variables, need to define a conversion factor
        varout{end+1} = 'conversion_factor';
        varout{end+1} = conversion_factor;
    end

   
    % send out an alignment matrix (and bias) if defined
    if exist('alignment_matrix_cxf','var') && ~isempty(alignment_matrix_cxf)
        varout{end+1} = 'alignment_matrix_cxf';
        varout{end+1} = alignment_matrix_cxf{ndset};
    end
    if exist('alignment_bias_c','var') && ~isempty(alignment_bias_c)
        varout{end+1} = 'alignment_bias_c';
        varout{end+1} = LFADS.Utils.makecol(alignment_bias_c{ndset});
    end

    % also store down the train and test inds, for posterity
    varout{end+1} = 'train_inds';
    varout{end+1} = trainInds;
    varout{end+1} = 'valid_inds';
    varout{end+1} = testInds;

    %% export the continuous data
    lfadsi_export_contdata(outfile, ytrain, ytest, varout{:})
end
prog.finish();


% tmp = comparisonDefaults();
% lfads_python_path = tmp.lfads_run_path_no_controller;

