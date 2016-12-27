function [mat, rtPredict] = threshPC1(data, rt)
% data is D x T x nTrials

% prep for PCA (T*nTrials) x nFactors
pcaMat = TensorUtils.reshapeByConcatenatingDims(data, {[2 3], 1});
pcaMat = bsxfun(@rdivide, pcaMat, range(pcaMat, 1));

[~, score, ~] = pca(pcaMat);
pcaFactors = reshape(score', size(data));

LFCIS = squeeze(pcaFactors(1, :, :))';
if mean(LFCIS(:, 1)) > mean(LFCIS(:, end))
    LFCIS = -LFCIS;
end
thresh = 0.3;
threshHigh = 0.7;
LFCIS = TensorUtils.rescaleIntervalToInterval(LFCIS, [mean(LFCIS(:, 1)), max(mean(LFCIS, 1))]);

mat = LFCIS';

[~, idxCross] = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(LFCIS', 1:size(LFCIS, 1), thresh, threshHigh);

X = [idxCross' onesvec(numel(idxCross))];
b  = regress(rt, X);

rtPredict = X * b;
