function [res_lfads, res_neural] = decodeKinematicsFromLFADSFactors(seq, pm, varargin)

p = inputParser();
p.addParameter('startIndHandKinematics', 400, @isscalar); % ms assuming 1 ms binning
p.addParameter('endIndHandKinematics', 800, @isscalar); % ms assuming 1 ms binning
p.addParameter('neuralKinematicsLag', 120, @isscalar); % ms
p.parse(varargin{:});
startInd = p.Results.startIndHandKinematics;
endInd = p.Results.endIndHandKinematics;

binSizeMs = pm.params.spikeBinMs;
startInd = floor(startInd/binSizeMs) * binSizeMs;
endInd = floor(endInd/binSizeMs) * binSizeMs;

% move hand kinematics to sequence structs
seq_raw=struct;
for ntr = 1:numel(seq)
    seq_raw(ntr).y = seq(ntr).y(:, startInd:endInd);
    seq_raw(ntr).X = seq(ntr).handKinematics(1:2,startInd:endInd);
    seq_raw(ntr).T = size(seq_raw(ntr).y, 2);
end

seq_lfads=seq_raw;
% create a seq where seq.y is replaced by the rates
for ntr = 1:numel(seq)
    seq_lfads(ntr).y = squeeze(pm.rates(:,(startInd/binSizeMs):(endInd/binSizeMs),ntr));
end

%% fit a kalman filter to neural data
opts = struct();
opts.neuralBinSizeMS = 1;
opts.lag = p.Results.neuralKinematicsLag;
T_raw = makeTfromR(seq_raw, binSizeMs, opts);
% fit a kalman filter
model = fitKalman(T_raw);
model = calcSteadyStateKalmanGain(model);
% actually run the decoder
T_raw_dec = kalmanIterative(model, T_raw);
% calculate R^2 on the results
xtrue = [T_raw.X];
xdec = [T_raw_dec.xk];
mdl_x = fitlm(xtrue(3,:), xdec(3,:));
mdl_y = fitlm(xtrue(4,:), xdec(4,:));

res_neural.r2X = mdl_x.Rsquared.Ordinary;
res_neural.r2Y = mdl_y.Rsquared.Ordinary;
res_neural.model = model;
res_neural.T = T_raw;
res_neural.T_dec = T_raw_dec;
res_neural.mdl_x = mdl_x;
res_neural.mdl_y = mdl_y;

%% fit a kalman filter to lfads data
opts = struct();
opts.neuralBinSizeMS = binSizeMs;
opts.lag = p.Results.neuralKinematicsLag;
T_lfads = makeTfromR(seq_lfads, binSizeMs, opts);
% fit a kalman filter
model = fitKalman(T_lfads);
model = calcSteadyStateKalmanGain(model);
% actually run the decoder
T_lfads_dec = kalmanIterative(model, T_lfads);
% calculate R^2 on the results
xtrue = [T_lfads.X];
xdec = [T_lfads_dec.xk];
mdl_x = fitlm(xtrue(3,:), xdec(3,:));
mdl_y = fitlm(xtrue(4,:), xdec(4,:));

res_lfads.r2X = mdl_x.Rsquared.Ordinary;
res_lfads.r2Y = mdl_y.Rsquared.Ordinary;
res_lfads.model = model;
res_lfads.T = T_lfads;
res_lfads.T_dec = T_lfads_dec;
res_lfads.mdl_x = mdl_x;
res_lfads.mdl_y = mdl_y;