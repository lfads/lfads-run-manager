function [W, b, lambda, mse] = ridge_cv(Y, X, varargin)
% fits the regression problem Y = X * W + b where coefficients in W are L^2 regularized.

% where Y is T x K outputs
%       X is T x N inputs
%       W is N x K weight matrix to map from X --> Y
%       b is 1 x K bias vector

% uses L2 penalty (ridge regression) and then k-fold cross-validation to
% optimize the ridge hyperparameter
    
    p = inputParser();
    p.addParameter('lambdas',  logspace(-5, 0, 15), @isvector);
    p.addParameter('KFold', 10, @isscalar);
    p.addParameter('normalizeEach', false, @isscalar); % zscore each column individually
    p.addParameter('plotFitAllLambdas', false, @islogical);
    p.addParameter('plotWeightsEachLambda', false, @islogical);
    
    p.parse(varargin{:});
    
    T = size(X, 1);
    assert(size(Y, 1) == T, 'X and Y must have same number or rows (observations)');
    
    % center and scale overall or per-channel variance of X
    muX = mean(X, 1, 'omitnan');
    muY = mean(Y, 1, 'omitnan');
    if p.Results.normalizeEach
        divisors = std(X, [], 1, 'omitnan');
    else
        divisors = std(X(:), 'omitnan');
    end
    X = (X - muX) ./ divisors;
    Y = Y - muY;
    
    % remove all nan-time points
    nanmask = any(isnan(Y), 2) | any(isnan(X), 2);
    X(nanmask, :) = [];
    Y(nanmask, :) = [];
    T = size(X, 1);
    
    %% pick best lambda
    lambdas = p.Results.lambdas;
    nL = numel(lambdas);
    
    if nL > 1
        nFolds = p.Results.KFold;
        cv = cvpartition(T, 'KFold', nFolds);

        loss = nan(nL, nFolds);
        for iL = 1:nL
            for iF = 1:nFolds
%                 maskTrain = cv.training(iF);
                maskTrain = mod((1:T)', nFolds) == iF-1;
                loss(iL, iF) = loss_single(Y, X, maskTrain, lambdas(iL));
            end
        end

        mse = sum(loss, 2);
        clf;
        plot(log10(lambdas), mse);
        xlabel('Log(10) lambda');
        ylabel('mse');
        
        [~, idxBestLambda] = min(mse);
        lambda = lambdas(idxBestLambda);
    else
        lambda = lambdas(1);
    end
    
    %% Plot fit to Y(:, 1) with all lambdas if requested
    if p.Results.plotFitAllLambdas
        y = Y(:, 1);
        
        figure(1);
        clf;
        plot(y, 'k-');
        hold on;
        
        cmap = parula(nL);
        for iL = 1:nL
            w = fit_single(y, X, lambdas(iL));
            yhat = X * w;
            
            if iL == idxBestLambda
                plot(yhat, 'Color', 'r');
            else
                plot(yhat, 'Color', cmap(iL, :));
            end
        end
        hold off;
    end
    
    if p.Results.plotWeightsEachLambda
        wnorm = nanvec(nL);
        y = Y(:, 1);
        
        figure(2);
        clf;
        
        cmap = parula(nL);
        for iL = 1:nL
            w = fit_single(y, X, lambdas(iL));
            
            if iL == idxBestLambda
                stem((1:numel(w)) + iL/nL*0.4, w, 'Color', 'r');
            else
                stem((1:numel(w)) + iL/nL*0.4, w, 'Color', cmap(iL, :));
            end
            hold on;
            wnorm(iL) = norm(w);
        end
        hold off;
        
        figure(3);
        clf;
        ax = plotyy(log10(lambdas), wnorm, log10(lambdas), mse);
        xlabel('Log(10) lambda');
        ylabel(ax(1), 'weights norm');
        ylabel(ax(2), 'cv mse');
    end
    
    %% Fit whole dataset with single lambda
    W = fit_single(Y, X, lambda);
    W = W ./ divisors';
    
    % correct overall shrinkage of coefficients
%     correction = mean(Y ./ (X*W));
%     W = W * correction;
   
    b = muY - muX*W;
end

function W = fit_single(Y, X, lambda)
    I = eye(size(X, 2));
    
    % solve the problem Y = X*B
    % unregularized solution is B^ = inv(X'*X) * X' * Y
    % tikhonov regularized solution is B^ = inv(X'*X + lambda*I) * X' * Y
    % which maps to:
    W = (X'*X + lambda * I) \ (X' * Y); 
end

function mse = loss_single(Y, X, maskTrain, lambda)
    W = fit_single(Y(maskTrain, :), X(maskTrain, :), lambda);
    
    maskTest = ~maskTrain;
    Yhat = X(maskTest, :) * W;
    Ytest = Y(maskTest, :);
    mse = mean((Yhat(:) - Ytest(:)).^2);
end