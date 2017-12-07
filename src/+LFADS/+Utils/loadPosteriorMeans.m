function pms = loadPosteriorMeans(validFile, trainFile, validInds, ...
                                  trainInds)
%function pms = loadPosteriorMeans(validFile, trainFile, validInds, ...
%                                  trainInds)

total_trials = numel(validInds) + numel(trainInds);

outputFields = {'controller_outputs','factors', 'rates', 'generator_states', 'generator_ics','costs','nll_bound_vaes','nll_bound_iwaes'};
storedVariables = {'controller_outputs','factors','output_dist_params', 'gen_states', 'gen_ics','costs','nll_bound_vaes','nll_bound_iwaes'};

tfInfo = h5info(trainFile);
tfNames = {tfInfo.Datasets.Name};
vaInfo = h5info(validFile);
vaNames = {vaInfo.Datasets.Name};

if numel(intersect(tfNames, vaNames)) ~= numel(tfNames)
    error('training and validation have different data?');
end

[vars_to_get, name_index, ~] = intersect(storedVariables, tfNames);
for nf = 1:numel(vars_to_get)
    variable = vars_to_get{nf};
    outputName = outputFields{name_index(nf)};

    valid.(outputName) = h5read(validFile, sprintf('/%s',variable));
    train.(outputName) = h5read(trainFile, sprintf('/%s',variable));
    
    % put them in the correct order
    if ismember(variable, {'costs','nll_bound_vaes','nll_bound_iwaes'})
        % trials on dim 1
        has_valid = size(valid.(outputName),1);
        has_train = size(train.(outputName),1);
    
        pms.(outputName) = nan(1, total_trials);
        pms.(outputName)(1,validInds(1:has_valid)) = valid.(outputName);
        pms.(outputName)(1,trainInds(1:has_train)) = train.(outputName);
    elseif ismember(variable, {'gen_ics'})
        % 2d, trials on dim 2
        has_valid = size(valid.(outputName),2);
        has_train = size(train.(outputName),2);
    
        pms.(outputName) = nan(size(train.(outputName),1), ...
                           total_trials);
        pms.(outputName)(:,validInds(1:has_valid)) = valid.(outputName);
        pms.(outputName)(:,trainInds(1:has_train)) = train.(outputName);
    else
        % 3d, trials on dim 3
        has_valid = size(valid.(outputName),3);
        has_train = size(train.(outputName),3);
    
        pms.(outputName) = nan(size(train.(outputName),1), ...
                           size(train.(outputName),2),...
                           total_trials);
        pms.(outputName)(:,:,validInds(1:has_valid)) = valid.(outputName);
        pms.(outputName)(:,:,trainInds(1:has_train)) = train.(outputName);
    end
end

pms.validInds = validInds';
pms.trainInds = trainInds';
