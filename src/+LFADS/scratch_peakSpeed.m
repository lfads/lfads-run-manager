%% Try regressing against peak speed

rng(0);
r = rc.runs(1);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

seq = seqData{1};
pm = pmData{1};

% nTrials x 1
peakSpeed = cat(1, seq.peakSpeed3);

T = size(pm.generator_states, 2);
rho = nanvec(T);
for iT = 1:T
    % will be nTrials x nGenUnits
    ics = squeeze(pm.generator_states(:, iT, :))';
    % % use ICS instead
    % icsTwoDay = pm.generator_ics';

    mdl = fitrlinear(ics, peakSpeed, 'KFold', 10);

    predPeakSpeed = kfoldPredict(mdl);

    rho(iT) = corr(predPeakSpeed, peakSpeed);
%     debug('T = %d, rho = %g\n', T, rho);
end

%% From ICs
ics = pm.generator_ics';
mdl = fitrlinear(ics, peakSpeed, 'KFold', 10);
predPeakSpeed = kfoldPredict(mdl);
rhoIC = corr(predPeakSpeed, peakSpeed);