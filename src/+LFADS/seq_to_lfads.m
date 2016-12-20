function seq_to_lfads(seqs, outpath, outfiles, varargin)
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
p.parse(varargin{:});
trainInds = p.Results.trainInds;
testInds = p.Results.testInds;
binSizeMS = p.Results.binSizeMs;
inputBinSizeMS = p.Results.inputBinSizeMs;
conversion_factor = p.Results.conversion_factor;
alignment_matrix_cxf = p.Results.alignment_matrix_cxf;

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

if ~exist('whichChannels','var'), 
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

% can we rebin data or is there a mismatch?
rebin = binSizeMS / inputBinSizeMS;
if mod(rebin, 1) ~= 0
    error('output binsize is not divisible by input binsize');
end

% create any dirs that don't exist
dirsToCreate = {spikesDir};
for nn = 1:numel(dirsToCreate)
    if ~isdir(dirsToCreate{nn})
        mkdir(dirsToCreate{nn});
    end
end

for ndset = 1:numel(seqs)
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
    nTimeMS = unique(arrayfun(@(x) size(x.y,2), seq)) * inputBinSizeMS;
    assert(numel(nTimeMS)==1, ['trials do not have constant ' ...
                        'length']);
    % neural dimension
    nNeurons = numel(whichChannelsThisSet);

    % how many output time bins should there be
    nTimeBins = floor(nTimeMS / binSizeMS);
    % how many input timebins should we keep
    inputTimeBinsToKeep = nTimeBins * binSizeMS / inputBinSizeMS;

    nTrainTrials = numel(trainInds);
    nTestTrials = numel(testInds);

    % initialize output variables
    ytrain = zeros(nTrainTrials, nTimeBins, nNeurons);
    ytest = zeros(nTestTrials, nTimeBins, nNeurons);
    if isfield(seq,'x_true')
        xtrain_true = zeros(nTrainTrials, nTimeBins, size(seq(1).x_true,1));
        xtest_true = zeros(nTestTrials, nTimeBins, size(seq(1).x_true,1));
    end
    if isfield(seq,'y_true')
        ytrain_true = zeros(nTrainTrials, nTimeBins, nNeurons);
        ytest_true = zeros(nTestTrials, nTimeBins, nNeurons);
    end
    
    % compile the train data
    for it = 1:nTrainTrials
        nn = trainInds(it);
        spks = seq(nn).y(whichChannelsThisSet,1:inputTimeBinsToKeep);
        if binSizeMS ~= inputBinSizeMS
            tmp = reshape(full(spks'),[],nTimeBins,size(spks,1));
            tmp2=squeeze(sum(tmp));
        else
            tmp2 = full(spks)';
        end
        ytrain(it,:,:) = int64(tmp2);
        % store down the (optional) true latents and FRs
        if isfield(seq,'x_true')
            xkeep = seq(nn).x_true(:, 1:inputTimeBinsToKeep);
            tmp = reshape(xkeep',[],nTimeBins,size(xkeep,1));
            tmp2=squeeze(sum(tmp))/binSizeMS;
            xtrain_true(it,:,:) = tmp2;
        end
        if isfield(seq,'y_true')
            ykeep = seq(nn).y_true(whichChannelsThisSet, 1:inputTimeBinsToKeep);
            tmp = reshape(ykeep',[],nTimeBins,size(ykeep,1));
            tmp2=squeeze(sum(tmp))/binSizeMS;
            ytrain_true(it,:,:) = tmp2;
        end
    end

    % compile the test data
    for it = 1:nTestTrials
        nn = testInds(it);
        spks = seq(nn).y(whichChannelsThisSet, 1:inputTimeBinsToKeep);
        if binSizeMS ~= inputBinSizeMS
            tmp = reshape(full(spks'),[],nTimeBins,size(spks,1));
            tmp2=squeeze(sum(tmp));
        else
            tmp2 = full(spks)';
        end
        ytest(it,:,:) = int64(tmp2);
        % store down the (optional) true latents and FRs
        if isfield(seq,'x_true')
            xkeep = seq(nn).x_true(:,1:inputTimeBinsToKeep);
            tmp = reshape(xkeep',[],nTimeBins,size(xkeep,1));
            tmp2=squeeze(sum(tmp))/binSizeMS;
            xtest_true(it,:,:) = tmp2;
        end
        if isfield(seq,'y_true')
            ykeep = seq(nn).y_true(whichChannelsThisSet, 1:inputTimeBinsToKeep);
            tmp = reshape(ykeep',[],nTimeBins,size(ykeep,1));
            tmp2=squeeze(sum(tmp))/binSizeMS;
            ytest_true(it,:,:) = tmp2;
        end        
    end

%     fprintf('Each trial is %g bins\n', nTimeBins);

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
        varout{end+1} = conversion_factor / rebin;
    end

    % send out an alignment matrix if defined
    if exist('alignment_matrix_cxf','var')
        varout{end+1} = 'alignment_matrix_cxf';
        varout{end+1} = alignment_matrix_cxf{ndset};
    end

    % also store down the train and test inds, for posterity
    varout{end+1} = 'train_inds';
    varout{end+1} = trainInds;
    varout{end+1} = 'valid_inds';
    varout{end+1} = testInds;

    %% export the spikes
    fprintf('Saving LFADS Input in %s\n', outfile);
    lfadsi_export_spikes(outfile, ytrain, ytest, varout{:})
end


% tmp = comparisonDefaults();
% lfads_python_path = tmp.lfads_run_path_no_controller;

