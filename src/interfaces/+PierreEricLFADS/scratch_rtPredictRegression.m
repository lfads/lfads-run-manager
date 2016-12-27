load('/Users/djoshea/data/lfadsCenterOut/rtPredict20161224.mat');
%%

for iR = 1:rc.nRuns
    rc.runs(iR).sequenceData = rc.runs(iR).modifySequenceDataPostLoading(rc.runs(iR).sequenceData);
end

%%

r = rc.runs(end);
r.addVelocityToSequenceData();
r.recomputeRTInSequenceData();
seq = r.sequenceData{7};
pm = r.posteriorMeans(7);

%%

r = rc.runs(7);
r.addVelocityToSequenceData();
r.recomputeRTInSequenceData();
seq = r.sequenceData{1};
pm = r.posteriorMeans(1);



%%

speedData = arrayfun(@(x) makecol(x.handSpeed), seq, 'UniformOutput', false);
speedTime = arrayfun(@(x) makecol(x.hand_time), seq, 'UniformOutput', false);

[speedMat, speedTime] = TrialDataUtilities.Data.embedTimeseriesInMatrix(speedData, speedTime);
speedMat = squeeze(speedMat);

rtIdx = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(speedMat, 1:numel(speedTime), 50, 100);

rt = speedTime(rtIdx);

plot(speedMat);
hold on;
plot(rtIdx, TensorUtils.selectSpecificIndicesAlongDimensionEachPosition(speedMat, 1, rtIdx), 'r.');

%%
% nFactors x T x nTrials
factors = pm.generator_states;

% prep for PCA (T*nTrials) x nFactors
pcaMat = TensorUtils.reshapeByConcatenatingDims(factors, {[2 3], 1});
pcaMat = bsxfun(@rdivide, pcaMat, range(pcaMat, 1));

[coeff, score, latent] = pca(pcaMat);
pcaFactors = reshape(score', size(factors));

%%
% nTrials x T
lfads_time = seq(1).y_time(1:10:end);
LFCIS = squeeze(pcaFactors(1, :, :))';
if mean(LFCIS(:, 1)) > mean(LFCIS(:, end))
    LFCIS = -LFCIS;
end
thresh = 0.5;
threshHigh = 0.7;
LFCIS = TensorUtils.rescaleIntervalToInterval(LFCIS, [mean(LFCIS(:, 1)), max(mean(LFCIS, 1))]);

[lfads_crossTime, lfads_idxCross] = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(LFCIS', lfads_time, thresh, threshHigh);
lfads_crossTime = lfads_crossTime';
lfads_idxCross = lfads_idxCross';

clf;
plot(lfads_time, LFCIS', 'k-');
hold on;
plot(lfads_crossTime, TensorUtils.selectSpecificIndicesAlongDimensionEachPosition(LFCIS, 2, lfads_idxCross), 'rx')

%%

mask = ~isnan(rt) & ~isnan(lfads_crossTime);
lfads_rho = corr(lfads_crossTime(mask), rt(mask));

clf
scatter(lfads_crossTime(mask), rt(mask))
lfads_rho

%% Try a velocity regression

% 100 x T x nTrials --> T x nTrials
startIdx = 25;
speedMat(isnan(speedMat)) = 0;
downsampledSpeedMat = resample(speedMat, 1, 10);
speedTrain = downsampledSpeedMat(startIdx:end, :);
lfadsTrain = pm.generator_states(:, startIdx:end, :);
[mdl, predY] = LFADS_PierreExport.fitrlinear_timeseries(lfadsTrain, speedTrain);

clf;
plot(speedTrain, 'k-');
hold on;
plot(predY, 'r-');

corr(speedTrain(:), predY(:))

%%
rtPredRegression = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(predY, 1:size(predY, 1), 50, 150);
rtPredRegression = rtPredRegression';

mask = ~isnan(rt) & ~isnan(rtPredRegression);
corr(rt(mask), rtPredRegression(mask))

% 0.78 prediction of rt, 0.89 prediction of velocity

%%

ra = rc.runs(end);
% prog = ProgressBar(rc.nRuns-1,'Predicting RT from LFADS');
for iR = 1:rc.nRuns-1
%     prog.update(iR);
    debug('Dataset %d\n', iR);
    [speedRho(iR), rtRho(iR)] = rc.runs(iR).predictSpeed(1);
    
    [speedRhoAll(iR), rtRhoAll(iR)] = ra.predictSpeed(iR);
end
% prog.finish();

%%



