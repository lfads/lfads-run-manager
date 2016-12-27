% take stitched data
pm = rc.runs(4).posteriorMeans;
pma = rc.runs(end).posteriorMeans(4);
seq = rc.runs(end).sequenceData{4};

% threshold pc 1
rtReal = cat(1, seq.rtThresh);
[mata, rtA] = PierreEricLFADS.RTPredict.threshPC1(pma.generator_states, rtReal);
[mat, rt] = PierreEricLFADS.RTPredict.threshPC1(pm.generator_states, rtReal);



%%
i = 55;

neural_time = seq(i).y_time(1:10:end);

clf;
plot(neural_time, mat(:, i), 'r-');
hold on;
plot(neural_time, mata(:, i), 'b-');

plot(seq(i).hand_time,  seq(i).handSpeed / 1000, 'k-');

vertLine(rtA(i), 'Color', 'r');
vertLine(rt(i), 'Color', 'b');
vertLine(rtReal(i), 'Color', 'k');
