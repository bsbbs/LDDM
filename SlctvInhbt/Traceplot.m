% mean trace of inhibitory neurons as a function of stimrate, sorted at
% stimOff
h = figure; hold on;
lowstim = stimrate == 6 | stimrate == 7;
efftrial = ~isnan(outcomes) & outcomes > -3;
inhbt = inhibitRois_pix == 1;
vals = stimOffAl.traces(:,inhbt,lowstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3),2);
plot(stimOffAl.time, meanval);

highstim = stimrate == 26 | stimrate == 27;
vals = stimOffAl.traces(:,inhbt,highstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3),2);
plot(stimOffAl.time, meanval);

% mean trace of inhibitory neurons as a function of stimrate, sorted at
% goTone
h = figure; hold on;
lowstim = stimrate == 6 | stimrate == 7;
efftrial = ~isnan(outcomes) & outcomes > -3;
inhbt = inhibitRois_pix == 1;
vals = goToneAl.traces(:,inhbt,lowstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3,'omitnan'),2);
plot(goToneAl.time, meanval);

highstim = stimrate == 26 | stimrate == 27;
vals = goToneAl.traces(:,inhbt,highstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3,'omitnan'),2);
plot(goToneAl.time, meanval);