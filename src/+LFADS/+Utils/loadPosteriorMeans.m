function pms = loadPosteriorMeans(validFile, trainFile, validInds, ...
                                  trainInds)
%function pms = loadPosteriorMeans(validFile, trainFile, validInds, ...
%                                  trainInds)

total_trials = numel(validInds) + numel(trainInds);

outputFields = {'controller_outputs','rates','factors', 'rates', 'generator_states', 'generator_ics'};
storedVariables = {'controller_outputs','output_dist_params','factors', ...
                   'rates', 'gen_states', 'gen_ics'};

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
    if ismember(variable, {'gen_ics'})
        % 2d, trials on dim 2
        has_valid = size(valid.(outputName),2);
        has_train = size(train.(outputName),2);
    
        pms.(outputName) = nan(size(train.(outputName),1), ...
                           total_trials);
        pms.(outputName)(:,validInds(1:has_valid)) = valid.(outputName);
        pms.(outputName)(:,trainInds(1:has_train)) = train.(outputName);
    else
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
