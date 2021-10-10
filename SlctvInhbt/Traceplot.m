% mean trace of inhibitory neurons as a function of stimrate, sorted at
% stimOff
h = figure; hold on;
lowstim = stimrate == 6 | stimrate == 7;
efftrial = ~isnan(outcomes) & outcomes > -1;
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
efftrial = ~isnan(outcomes) & outcomes > -1;
inhbt = inhibitRois_pix == 1;
vals = goToneAl.traces(:,inhbt,lowstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3,'omitnan'),2);
plot(goToneAl.time, meanval);

highstim = stimrate == 26 | stimrate == 27;
vals = goToneAl.traces(:,inhbt,highstim & efftrial'); % time, neuron, trials
meanval = mean(mean(vals,3,'omitnan'),2);
plot(goToneAl.time, meanval);

%% time series of regression weights of activity on stimrate
% sorted at stimOff
efftrial = ~isnan(outcomes) & outcomes > -1;
stimrateeff = stimrate(efftrial);
inhbt = inhibitRois_pix == 1;
vals = stimOffAl.traces(:,inhbt, efftrial'); % time, cells, trials
srates = unique(stimrateeff);
meantrc = [];
for ri = 1:numel(srates) % loop on stimrate
    srate = srates(ri);
    for ci = 1:size(vals,2) % loop on cells
        meantrc(:,ci,ri) = mean(vals(:,ci,stimrateeff == srate),3);
    end
end
b = []; % fitted slope of each cell at time t
for ci = 1:size(vals,2)
    for ti = 1:size(meantrc,1) % loop on time
        f = polyfit(srates,squeeze(meantrc(ti,ci,:)),1);
        %     h = figure; hold on;
        %     plot(srates,squeeze(meantrc(ti,ci,:)),'.');
        %     plot(srates,srates*f(1) + f(2),'-');
        b(ti,ci) = f(1);
    end
end
h = figure; hold on;
for ci = 1:size(vals,2)
    plot(b(:,ci),'-');
end
h = figure;
plot(mean(b,1),'.');

% sorted at goTone
vals = goToneAl.traces(:,inhbt, efftrial'); % time, cells, trials
meantrc = [];
for ri = 1:numel(srates) % loop on stimrate
    srate = srates(ri);
    for ci = 1:size(vals,2) % loop on cells
        meantrc(:,ci,ri) = mean(vals(:,ci,stimrateeff == srate),3);
    end
end
bg = []; % fitted slope of each cell at time t
for ci = 1:size(vals,2)
    for ti = 1:size(meantrc,1) % loop on time
        f = polyfit(srates,squeeze(meantrc(ti,ci,:)),1);
        %     h = figure; hold on;
        %     plot(srates,squeeze(meantrc(ti,ci,:)),'.');
        %     plot(srates,srates*f(1) + f(2),'-');
        bg(ti,ci) = f(1);
    end
end
h = figure; hold on;
plot(mean(b,1),mean(bg,1),'.');
xlabel('\beta before stimOff');
ylabel('\beta after goTone');
