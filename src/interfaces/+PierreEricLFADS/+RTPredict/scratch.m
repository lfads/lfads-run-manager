%%

load('/Users/djoshea/data/lfadsCenterOut/rt20161227.mat')
%%

% take stitched data
% pm = rc.runs(4).posteriorMeans;
% pma = rc.runs(end).posteriorMeans(4);
% seq = rc.runs(end).sequenceData{4};

pm = rc.runs(7).posteriorMeans;
pma = rc.runs(end).posteriorMeans(7);
seq = rc.runs(end).sequenceData{7};
r = rc.runs(7);
ra = rc.runs(end);

%%

r.predictRTFromLFADS(1, 'normalizedThresh', 0.5);

%%
ra.predictRTFromLFADS(7, 'normalizedThresh', 0.5);

%%
% threshold pc 1
rtReal = cat(1, seq.rtThresh);
[mata, rtA] = PierreEricLFADS.RTPredict.threshPC1(pma.generator_states, rtReal);
[mat, rt] = PierreEricLFADS.RTPredict.threshPC1(pm.generator_states, rtReal);

nancorr = @(x, y) corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)));

%%

i = 104;

neural_time = seq(i).y_time(1:10:end);

clf;
plot(neural_time, mat(:, i), 'r-');
hold on;
plot(neural_time, mata(:, i), 'b-');

plot(seq(i).hand_time,  seq(i).handSpeed / 1000, 'k-');

vertLine(rtA(i), 'Color', 'r');
vertLine(rt(i), 'Color', 'b');
vertLine(rtReal(i), 'Color', 'k');
