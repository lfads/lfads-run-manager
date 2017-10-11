function datasets = generateDemoDatasets(datasetPath, nDatasets)
    % generates a set of demo datasets for LFADS based on generating
    % spikes via a chaotic Lorenz attractor
    
    rng(1);

    %% Create temp directory for datasets
    if ~exist(datasetPath, 'dir')
        mkdir(datasetPath);
    end

    %% Generate datasets
    minChannels = 25;
    maxChannels = 35;

    nConditions = 65;
    minTrialsC = 20;
    maxTrialsC = 30;
    T = 1000;

    nTrialsCByDataset = randi([minTrialsC, maxTrialsC], nDatasets);
    nChannelsByDataset = randi([minChannels, maxChannels], nDatasets);

    meanFr = 5;
    D = 3;
    
    initialConditions = randn(3, nConditions) * 2; % random IC per condition, but constant across datasets for stitching
    lorenz_trajectories = nan(3, T, nConditions);
    for iC = 1:nConditions
        lorenz_trajectories(:, :, iC) = lorenz(T, initialConditions(:, iC));
    end
    
    for iDS = 1:nDatasets
        nTrC = nTrialsCByDataset(iDS);
        nCh = nChannelsByDataset(iDS);

        W = sort((rand(nCh, D)+1) .* sign(randn(nCh, D)), 1); % N x 3
        b = log(meanFr / 1000) * ones(nCh, 1); % N x 1

        trialIdx = 1;
        spikes = nan(nTrC * nConditions, nCh, T);
        conditionId = nan(nTrC * nConditions, 1);

        for iC = 1:nConditions
            % generate rates via lorenz
            X = lorenz_trajectories(:, :, iC); % 3 X T
            log_rates = W*X + b; % Nx3 x 3xT = N x T
            rates = exp(clip(log_rates, -20, 20));
            
            for iTr = 1:nTrC
                spikes(trialIdx, :, :) = poissrnd(rates) > 0;
                conditionId(trialIdx) = iC;
                trialIdx = trialIdx + 1;
            end
        end

        datasets(iDS).spikes = spikes;
        datasets(iDS).timeMs = (0:T-1)';
        datasets(iDS).conditionId = conditionId;
        datasets(iDS).datetime = datetime('today') - nDatasets + iDS;
        datasets(iDS).subject = 'lorenz_example';
        
        datasets(iDS).W = W;
        datasets(iDS).b = b;
        datasets(iDS).lorenz_trajectories = lorenz_trajectories;
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

function x = clip(x, lo, hi)
    x(x < lo) = lo;
    x(x > hi) = hi;
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