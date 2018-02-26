function datasets = generateDemoDatasets(datasetPath, varargin)
    % generates a set of demo datasets for LFADS based on generating
    % spikes via a chaotic Lorenz attractor
    
    p = inputParser();
    p.addParameter('nDatasets', 3, @isscalar);
    p.addParameter('minChannels', 25, @isscalar);
    p.addParameter('maxChannels', 35, @isscalar);
    p.addParameter('nConditions', 65, @isscalar);
    p.addParameter('minTrialsC', 20, @isscalar); % per condition
    p.addParameter('maxTrialsC', 30, @isscalar);
    p.addParameter('nTime', 1000, @isscalar);
    p.addParameter('meanFr', 5, @isscalar);
    p.addParameter('T_burnIn', 500, @isscalar);
    p.parse(varargin{:});
    
    %% Create temp directory for datasets
    if ~exist(datasetPath, 'dir')
        mkdir(datasetPath);
    end

    %% Generate datasets
    nDatasets = p.Results.nDatasets;
    minChannels = p.Results.minChannels;
    maxChannels = p.Results.maxChannels;

    nConditions = p.Results.nConditions;
    minTrialsC = p.Results.minTrialsC;
    maxTrialsC = p.Results.maxTrialsC;
    T = p.Results.nTime;
    T_burnIn = p.Results.T_burnIn;
    
    meanFr = p.Results.meanFr;
    D = 3; % dimensionality of lorenz system
    
    s = RandStream('mt19937ar','Seed', 0);

    nTrialsCByDataset = randi(s, [minTrialsC, maxTrialsC], nDatasets);
    nChannelsByDataset = randi(s, [minChannels, maxChannels], nDatasets);
    initialConditions = bsxfun(@plus, randn(s, D, nConditions), [0 0 25]'); % random IC per condition, but constant across datasets for stitching
    
    % run the lorenz system for a short time to bring the initial
    % conditions into the butterfly
    if T_burnIn > 1
        for iC = 1:nConditions
            traj = lorenz(T_burnIn, initialConditions(:, iC));
            initialConditions(:, iC) = traj(:, end);
        end
    end

    % generate the lorenz trajectories
    lorenz_trajectories = nan(D, T, nConditions);
    for iC = 1:nConditions
        lorenz_trajectories(:, :, iC) = lorenz(T, initialConditions(:, iC));
    end
    
    % transform to [0 1] range
    lorenz_trajectories = lorenz_trajectories - mean(lorenz_trajectories, 2);
    lorenz_trajectories = lorenz_trajectories ./ max(lorenz_trajectories, [], 2);
    
    figure();
    plotTrajectories(lorenz_trajectories);
    drawnow;
    
    % project the trajectories into each dataset's neurons
    for iDS = 1:nDatasets
        nTrC = nTrialsCByDataset(iDS);
        nCh = nChannelsByDataset(iDS);

        W = sort((rand(s, nCh, D)+1) .* sign(randn(s, nCh, D)), 1); % N x 3
        b = log(meanFr / 1000) * ones(nCh, 1); % N x 1

        trialIdx = 1;
        spikes = nan(nTrC * nConditions, nCh, T);
        conditionId = nan(nTrC * nConditions, 1);
        true_rates = nan(nTrC * nConditions, nCh, T);

        for iC = 1:nConditions
            % generate rates via lorenz
            X = lorenz_trajectories(:, :, iC); % 3 X T
            log_rates = bsxfun(@plus, W*X, b); % Nx3 x 3xT = N x T
            rates = exp(clip(log_rates, -20, 20));
            
            for iTr = 1:nTrC
                true_rates(trialIdx, :, :) = rates;
                spikes(trialIdx, :, :) = poissrnd(rates) > 0;
                conditionId(trialIdx) = iC;
                trialIdx = trialIdx + 1;
            end
        end

        datasets(iDS).spikes = spikes; %#ok<*AGROW>
        datasets(iDS).timeMs = (0:T-1)'; 
        datasets(iDS).conditionId = conditionId;
        datasets(iDS).datetime = datetime('today') - nDatasets + iDS;
        datasets(iDS).subject = 'lorenz_example';
        
        datasets(iDS).W = W;
        datasets(iDS).b = b;
        datasets(iDS).lorenz_trajectories = lorenz_trajectories;
        datasets(iDS).true_rates = true_rates;
    end

    %% Save datasets to disk

    for iDS = 1:nDatasets
        fname = fullfile(datasetPath, sprintf('dataset%03d.mat', iDS));
        dataset = datasets(iDS); %#ok<NASGU>
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
end

function plotTrajectories(lorenz_trajectories)
    % lorenz_trajectories is 3 x T x C
    
    cla;
    ic = squeeze(lorenz_trajectories(:, 1, :));
    C = size(lorenz_trajectories, 3);
    [~, sortIdx] = sortrows(ic');
    
    lorenz_trajectories = lorenz_trajectories(:, :, sortIdx);
    ic = ic(:, sortIdx);
    cmap = parula(C);

    for iC = 1:C
        plot3(lorenz_trajectories(1, :, iC), lorenz_trajectories(2, :, iC), lorenz_trajectories(3, :, iC), ...
            'Color', [cmap(iC, :) 0.7], 'LineWidth', 0.5);
        hold on;
    end
    
    plotICs(ic);
    
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(28.9, 7.6);
    axis equal;
    axis vis3d;
    axis off;
end
    
function plotICs(ic)
    scatter3(ic(1, :), ic(2, :), ic(3, :), 'LineWidth', 0.5, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5);
end

