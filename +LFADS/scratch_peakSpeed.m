%% Try regressing against peak speed

rng(0);
r = rc.runs(1);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

seq = seqData{1};
pm = pmData{1};

% nTrials x 1
peakSpeed = cat(1, seq.peakSpeed3);

% nTrials x nGenUnits
% ics = squeeze(pm.generator_states(:, 3, :))';
% use ICS instead
icsTwoDay = pm.generator_ics';

mdlTwoDay = fitrlinear(icsTwoDay, peakSpeed, 'KFold', 10);

debug('TwoDay model: %g\n', kfoldLoss(mdlTwoDay));

cvpred_peakSpeed3_twoDay = kfoldPredict(mdlTwoDay);
scatter(cvpred_peakSpeed3_twoDay, peakSpeed);

debug('TwoDay rho = %g\n', corr(cvpred_peakSpeed3_twoDay, peakSpeed));

%%
rng(0);
r = rc.runs(2);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

% seq = seqData{1};
pm = pmData{1};

% % nTrials x 1
% peakSpeed3 = cat(1, seq.peakSpeed3);

% nTrials x nGenUnits
% ics = squeeze(pm.generator_states(:, 3, :))';
ics0921 = pm.generator_ics';

mdl0921 = fitrlinear(ics0921, peakSpeed, 'KFold', 10);

debug('One0921 model: %g\n', kfoldLoss(mdl0921));

cvpred_peakSpeed3_oneDay0921 = kfoldPredict(mdl0921);
hold on;
scatter(cvpred_peakSpeed3_oneDay0921, peakSpeed3);

debug('One0921 rho = %g\n', corr(cvpred_peakSpeed3_oneDay0921, peakSpeed3));

% 
% %% predict separately per condition
% 
% cnames = {'DownRight', 'Right', 'UpRight', 'Up', 'UpLeft', 'Left', 'DownLeft'};
% [~, cond] = ismember({seq.targetDirectionName}, cnames);
% 
% nC = max(cond);
% [loss0921, lossTwoDay] = nanvec(nC);
% for iC = 1:nC
%     mask = cond == iC;
%     mdl1 = fitrlinear(ics0921(mask, :), peakSpeed(mask), 'KFold', 10);
%     loss0921(iC) = kfoldLoss(mdl1);
%     mdl2 = fitrlinear(icsTwoDay(mask, :), peakSpeed(mask), 'KFold', 10);
%     lossTwoDay(iC) = kfoldLoss(mdl2);
% end
% 
% debug('Two day %g vs one day %g\n', sum(lossTwoDay), sum(loss0921));
% 


