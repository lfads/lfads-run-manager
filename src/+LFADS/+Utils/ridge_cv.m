function [W, b] = ridge_cv(Y, X, varargin)
% fits the regression problem Y = X * W + b where coefficients in W are L^2 regularized.

% where Y is T x K outputs
%       X is T x N inputs
%       W is N x K weight matrix to map from X --> Y
%       b is 1 x K bias vector

% uses L2 penalty (ridge regression) and then k-fold cross-validation to
% optimize the ridge hyperparameter

    
    p = inputParser();
    p.addParameter('lambdas',  logspace(-5,0,15), @isvector);
    p.addParameter('KFold', 5, @isscalar);
    p.addParameter('normalizeEach', false, @isscalar); % zscore each column individually
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
                maskTrain = cv.training(iF);
                loss(iL, iF) = loss_single(Y(maskTrain, :), X(maskTrain, :), lambdas(iL));
            end
        end

        mse = sum(loss, 2);
        [~, idxBestLambda] = min(mse);
        lambda = lambdas(idxBestLambda);
    else
        lambda = lambdas(1);
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

function mse = loss_single(Y, X, lambda)
    W = fit_single(Y, X, lambda);
%     W = ridge(Y, X, lambda, 0); W = W(2:end);
    Yhat = X * W;
    mse = sum((Yhat(:) - Y(:)).^2);
end