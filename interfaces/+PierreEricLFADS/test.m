dc = PierreEricLFADS.DatasetCollection('/data2/lfads/PierreEric/export_v02_broadbandRethreshNonSorted');
dc.autoDetectDatasets();

dc.loadInfo;
dc.filterHasHighSNRChannels();
dc.filterBestSaveTagEachDate();

%%

par = PierreEricLFADS.RunParams;
par.spikeBinMs = 10;
par.batchSize = 40;
par.nTrialsKeep = 500;
par.regularizerIncreaseSteps = 100;
par.learningRateDecayFactor = 0.95;

rc = PierreEricLFADS.RunCollection('/data2/lfads/PierreEric/runs', 'oneVsTwo_0922_0921', dc, par);

%%
r = PierreEricLFADS.Run('two0921_0922', rc);
r.selectDatasetsByName({...
    'subject_Pierre.date_2016-09-21.saveTagGroup_2_export', ...
    'subject_Pierre.date_2016-09-22.saveTagGroup_1_export'});
r.prepareForLFADS();

%%

r = PierreEricLFADS.Run('one0921', rc);
r.selectDatasetsByName('subject_Pierre.date_2016-09-21.saveTagGroup_2_export');



%%

r = LFADS.Run('one0922');
rc.addRun(r);
r.selectDatasets(7);
r.prepareForLFADS();

%%

r = LFADS.Run('four0921_0922_0923_0926');
rc.addRun(r);
r.selectDatasets([6 7 8 9]);
r.prepareForLFADS();

%%

r = LFADS.Run('all');
rc.addRun(r);
r.selectDatasets(1:12);
r.prepareForLFADS();


%%
rc.writeTensorboardShellScript

%%

r = rc.runs(1);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

seq = seqData{1};
pm = pmData{1};

% which condition is each trial
[cnames, ~, cond] = unique({seq.targetDirectionName});
nC = numel(cnames);

%% plot the factors colored by condition

clf
cmap = TrialDataUtilities.Color.hslmap(nC);
% factors is nFactors x nTimeBins x nTrials
nFactors = size(pm.factors, 1);
for iF = 1:nFactors
    subtightplot(nC, 1, iF);
    for iC = 1:nC
       h = plot(squeeze(pm.factors(iF, :, cond == iC)), 'Color', cmap(iC, :));
       TrialDataUtilities.Plotting.setLineOpacity(h, 0.3);
       hold on;
       axis tight, box off;
    end
end


%% plot the controller outputs colored by condition

clf
cmap = TrialDataUtilities.Color.hslmap(nC);
% factors is nFactors x nTimeBins x nTrials
nInputs = size(pm.controller_outputs, 1);
for iF = 1:nInputs
    subtightplot(nInputs, 1, iF);
    for iC = 1:nC
       h = plot(squeeze(pm.controller_outputs(iF, :, cond == iC)), 'Color', cmap(iC, :));
       TrialDataUtilities.Plotting.setLineOpacity(h, 0.3);
       hold on;
       axis tight, box off;
    end
end

%% try a decode using two day

r = rc.runs(1);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

seq = seqData{1};
pm = pmData{1};
[res_lfads_twoDay, res_neural_twoDay] = LFADS.decodeKinematicsFromLFADSFactors(seq, pm);

figUnique(r.name);
LFADS.compareKinematicDecodes(seq, {res_lfads_twoDay, res_neural_twoDay}, 'colormap', {'b', 'r'});

%% decode using one day

r = rc.runs(2);

seqData = r.loadSequenceFiles();
pmData = r.loadPosteriorMeanSamples();

seq = seqData{1};
pm = pmData{1};
[res_lfads_0921, res_neural_0921] = LFADS.decodeKinematicsFromLFADSFactors(seq, pm);

figUnique(r.name);
LFADS.compareKinematicDecodes(seq, {res_lfads_0921, res_neural_0921}, 'colormap', {'b', 'r'});

%% grab generator state a few timesteps in as the intiial condition

% nGenUnits x nTrials
ics = squeeze(pm.generator_states(:, 3, :));

% which condition is each trial
cnames = {'DownRight', 'Right', 'UpRight', 'Up', 'UpLeft', 'Left', 'DownLeft'};
[~, cond] = ismember({seq.targetDirectionName}, cnames);

nC = numel(cnames);

clf;
% ics_tsne = tsne(ics', cond)';
ics_tsne = tsne(ics', [])'; % 2 x nTrials

for iC = 1:nC
   h = plot(ics_tsne(1, cond == iC), ics_tsne(2, cond == iC), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmap(iC, :));
   TrialDataUtilities.Plotting.showFirstInLegend(h, cnames{iC});
   hold on;
end
legend(gca, 'show');


