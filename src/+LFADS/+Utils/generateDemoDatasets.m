function generateDemoDatasets(datasetPath, nDatasets)
    rng(1);

    %% Create temp directory for datasets
    if ~exist(datasetPath, 'dir')
        mkdir(datasetPath);
    end

    %% Generate datasets
    minChannels = 20;
    maxChannels = 30;

    nConditions = 50;
    minTrialsC = 20;
    maxTrialsC = 30;
    T = 1000;

    nTrialsCByDataset = randi([minTrialsC, maxTrialsC], nDatasets);
    nChannelsByDataset = randi([minChannels, maxChannels], nDatasets);

    meanFr = 5;
    D = 3;
    
    initialConditions = randn(3, nDatasets) * 2; % random IC per condition, but constant across datasets for stitching

    for iDS = 1:nDatasets
        nTrC = nTrialsCByDataset(iDS);
        nCh = nChannelsByDataset(iDS);

        W = (rand(nCh, D)+1) .* sign(randn(nCh, D)); % N x 3
        b = log(meanFr / 1000) * ones(nCh, 1); % N x 1

        trialIdx = 1;
        spikes = nan(nTrC * nConditions, nCh, T);
        conditionId = nan(nTrC * nConditions, 1);

        for iC = 1:nConditions
            % generate rates via lorenz
            x0 = initialConditions(:, iDS);
            X = lorenz(T, x0); % 3 x T

            rates = W*X + b; % N x T

            for iTr = 1:nTrC
                spikes(trialIdx, :, :) = poissonSpike(rates);
                conditionId(trialIdx) = iC;
                trialIdx = trialIdx + 1;
            end
        end

        datasets(iDS).spikes = spikes;
        datasets(iDS).timeMs = (0:T-1)';
        datasets(iDS).conditionId = conditionId;
        datasets(iDS).datetime = datetime('today') - nDatasets + iDS;
        datasets(iDS).subject = 'lorenz_example';
    end

    %% Save datasets to disk

    for iDS = 1:nDatasets
        fname = fullfile(datasetPath, sprintf('dataset%03d.mat', iDS));
        dataset = datasets(iDS);
        fprintf('Saving %s\n', fname);
        save(fname, '-struct', 'dataset');
    end

end

%% Utility functions

function spikes = poissonSpike(rates)
    exprates = exp(rates);
    exprates(exprates < -20) = -20;
    exprates(exprates >  20) =  20;
    
    spikes = poissrnd(exprates) > 0;
end

function X = lorenz(T, x0)
% Euler integration of lorenz attractor. X will be 3xT, x0 is 3x1

    r = 28;
    s = 10;
    b = 8/3;
    dt = 0.006;
    derivfn = @(x) [s * (x(2) - x(1)); ...
                    x(1) * (r - x(3)); ...
                    x(1) * x(2) - b * x(3) ];

    X = nan(3, T);
    X(:, 1) = x0;
    for t = 2:T
        X(:, t) = X(:, t-1) + dt * derivfn(X(:, t-1));
    end
    
    X = X - mean(X, 2);
    X = X ./ max(X, [], 2);
end