% Fig 8. adapt to different timeline
%% RT task, LDDM - parameter settings
a = a0*eye(2);
b = b0*eye(2);
w = ones(2);
predur = .8;
presentt = 0;
% %% Fig. 7a, Reaction time task - stepwised action signal
% dur = 2.1;
% stimdur = dur-presentt;
% triggert = presentt;
% sgm = 0;
% thresh = 25; %45;
% stoprule = 1;
% V0 = [30,30];
% h = figure;
% filename = 'Fig7II';
% subplot(4,1,1); hold on;
% for vi = 1:5
%     initialvals = zeros(3,2);
%     [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
%     Vinput = VmatDiag(vi,:);
%     initialvals = [R0(end,:);G0(end,:);I0(end,:)];
%     [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%     R = [R0; R];
%     lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
%     lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
% end
% plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
% text(1200,thresh*1.1,'threshold');
% 
% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+predur+presentt+dur)/dt],[0,0,2,2]-5,'k-','LineWidth',lwd/1.5); % target
% plot([1, (prepresentt+prestimdur+presentt)/dt, (prepresentt+prestimdur+presentt)/dt, (prepresentt+predur+presentt+dur)/dt],[0,0,2,2]-10,'k-','LineWidth',lwd/1.5); % motion
% plot([1, (prepresentt+prestimdur+triggert)/dt, (prepresentt+prestimdur+triggert)/dt, (prepresentt+predur+presentt+dur)/dt],[0,0,2,2]-15,'k-','LineWidth',lwd/1.5); % go signal
% ylim([-15,max(ylim)*1.1]);
% yticks([0]);
% yticklabels({'0'});
% ylabel('           Activity (a.u.)');
% xticks([]);
% drawaxis(gca, 'x', 0, 'movelabel', 1);
% xlim([1, (prepresentt+predur+presentt+dur)/dt]);
% % xticks([prepresentt/dt, (prepresentt+prestimdur)/dt]);
% % xticklabels({'Targets','Motion'});
% % xlabel('Time (a.u.)');
% xlabel(' ');
% legend([lgd1, lgd2],{'R_1', 'R_2'},...
%     'Location','NorthEast','FontSize',fontsize-5, 'FontName','Times New Roman', ...
%     'FontAngle','italic','NumColumns',1,'Box','off');
% savefig(h, filename, outdir, fontsize, aspect15);

%% Fig. 7aII, Hybrid model - Reaction time task, urgency signal.
dur = 2.1;
stimdur = dur - presentt;
triggert = presentt;
% triggert = dur;
% triggert1 = -3; %presentt+stimdur-.2;
% triggert2 = triggert;
% b1 = b;
% b2 = b;
sgm = 0;
thresh = 25; %45;
stoprule = 1;
V0 = [30,30];
h = figure;
filename = 'Fig7III';% hold on;
subplot(4,1,1); hold on;
% h = figure; hold on;
% filename = 'Fig7d';
for vi = 1:5
    initialvals = zeros(3,2);
    [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
    Vinput = VmatDiag(vi,:);
    initialvals = [R0(end,:);G0(end,:);I0(end,:)];
    
    % [R, ~, ~, ~, ~] = LcDsInhbtUrgencyGate(Vinput, w, a, b1, b2, sgm, Tau, dur, dt, presentt, triggert1, triggert2, thresh, initialvals, stimdur, stoprule);
    [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    R = [R0; R];
    lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
end
plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-'); % thresh
% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+prestimdur)/dt],...
%     [0,0,2,2]-7,'k-','LineWidth',lwd/1.5);
% plot([(predur)/dt, (predur+presentt+dur)/dt],...
%     [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs (target & stimuli)

plot([1, prepresentt/dt prepresentt/dt],...
    [0,0,2]-7,'k-','LineWidth',lwd/1.5); 
plot([prepresentt/dt, (prepresentt+prestimdur)/dt],...
    [2,2]-7,'k--','LineWidth',lwd/1.5); % input, target, pre-motion
plot([(predur)/dt, (predur+presentt+dur)/dt],...
    [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion


% plot([1, (prepresentt)/dt, (prepresentt)/dt, (predur)/dt, (predur)/dt, (predur+presentt+dur)/dt],[2,2,2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
% plot([1, (predur+presentt)/dt, (predur+triggert2)/dt, (predur+triggert2)/dt, (predur+presentt+dur)/dt],[0,0,3*b1(1)/b2(1),3,3]-16,'k-','LineWidth',lwd/1.5); % disinhibition
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-11.5,'k-','LineWidth',lwd/1.5); % disinhibition
ylim([-16.5,27]);
yticks([0]);
yticklabels({'0'});
ylabel('             Activity (a.u.)');
xticks([]);
drawaxis(gca, 'x', 0, 'movelabel', 1);
xlim([1, (prepresentt+predur+presentt+dur)/dt]);
xlabel(' ');
savefigs(h, filename, outdir, fontsize, aspect15);

%% Fig. 8b, LDDM - Fixed duration, no delay.
dur = 2.1;
stimdur = .8;
triggert = presentt+stimdur;
sgm = 0;
thresh = 25; %45;
stoprule = 1;
V0 = [30,30];
% h = figure; hold on;
% filename = 'Fig7IIb';
subplot(4,1,2); hold on;
for vi = 1:5
    initialvals = zeros(3,2);
    [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
    Vinput = VmatDiag(vi,:);
    initialvals = [R0(end,:);G0(end,:);I0(end,:)];
    [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    R = [R0; R];
    lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
end
plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
plot([1, prepresentt/dt prepresentt/dt],...
    [0,0,2]-7,'k-','LineWidth',lwd/1.5); 
plot([prepresentt/dt, (prepresentt+prestimdur)/dt],...
    [2,2]-7,'k--','LineWidth',lwd/1.5); % input, target, pre-motion
plot([(predur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+dur)/dt],...
    [2,2,0,0]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion
% 
% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+prestimdur)/dt],...
%     [0,0,2,2]-7,'k-','LineWidth',lwd/1.5);
% plot([(predur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+dur)/dt],...
%     [2,2,0,0]-7,'k-','LineWidth',lwd/1.5); % inputs (target & stimuli)
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-16,'k-','LineWidth',lwd/1.5); % disinhibition
ylim([-16.5,27]);
yticks([0]);
yticklabels({'0'});
ylabel('           Activity (a.u.)');
xticks([]);
drawaxis(gca, 'x', 0, 'movelabel', 1);
xlim([1, (prepresentt+predur+presentt+dur)/dt]);
xlabel(' ');
savefig(h, filename, outdir, fontsize, aspect15);

% Fig. 8c, LDDM - Fixed duration, with delay.
dur = 2.1;
stimdur = .8;
triggert = presentt+1.5;
sgm = 0;
thresh = 25; %45;
stoprule = 1;
V0 = [30,30];
subplot(4,1,3); hold on;
% h = figure; hold on;
% filename = 'Fig7IIc';
for vi = 1:5
    initialvals = zeros(3,2);
    [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
    Vinput = VmatDiag(vi,:);
    initialvals = [R0(end,:);G0(end,:);I0(end,:)];
    [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    R = [R0; R];
    lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
end
plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
plot([1, prepresentt/dt prepresentt/dt],...
    [0,0,2]-7,'k-','LineWidth',lwd/1.5); 
plot([prepresentt/dt, (prepresentt+prestimdur)/dt],...
    [2,2]-7,'k--','LineWidth',lwd/1.5); % input, target, pre-motion
%plot([(predur)/dt, (predur+presentt+dur)/dt],...
%     [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion
plot([(predur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+dur)/dt],...
    [2,2,0,0]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion
% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+prestimdur)/dt],...
%     [0,0,2,2]-7,'k-','LineWidth',lwd/1.5);
% plot([(predur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+dur)/dt],...
%     [2,2,0,0]-7,'k-','LineWidth',lwd/1.5); % inputs (target & stimuli)
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-16,'k-','LineWidth',lwd/1.5); % disinhibition
ylim([-16.5,27]);
yticks([0]);
yticklabels({'0'});
ylabel('           Activity (a.u.)');
xticks([]);
drawaxis(gca, 'x', 0, 'movelabel', 1);
xlim([1, (prepresentt+predur+presentt+dur)/dt]);
xlabel(' ');
savefig(h, filename, outdir, fontsize, aspect15);

% % Fig. 8d, Hybrid model - Fixed duration, with delay, preparation signal.
% dur = 2.1;
% stimdur = .8;
% triggert = presentt+1.5;
% triggert1 = 0;%presentt+stimdur-.2;
% triggert2 = triggert;
% b1 = eye(2)*.75;
% b2 = b;
% sgm = 0;
% thresh = 25; %45;
% stoprule = 1;
% V0 = [30,30];
% subplot(4,1,4); hold on;
% % h = figure; hold on;
% % filename = 'Fig7IId';
% for vi = 1:5
%     initialvals = zeros(3,2);
%     [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
%     Vinput = VmatDiag(vi,:);
%     initialvals = [R0(end,:);G0(end,:);I0(end,:)];
%     
%     [R, ~, ~, ~, ~] = LcDsInhbtUrgencyGate(Vinput, w, a, b1, b2, sgm, Tau, dur, dt, presentt, triggert1, triggert2, thresh, initialvals, stimdur, stoprule);
%     % [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%     R = [R0; R];
%     lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
%     lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
% end
% plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+prestimdur)/dt],...
%     [0,0,2,2]-7,'k-','LineWidth',lwd/1.5);
% plot([(predur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+stimdur)/dt, (predur+presentt+dur)/dt],...
%     [2,2,0,0]-7,'k-','LineWidth',lwd/1.5); % inputs (target & stimuli)
% plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
% plot([1, (predur+presentt)/dt, (predur+presentt)/dt, (predur+triggert2)/dt, (predur+triggert2)/dt, (predur+presentt+dur)/dt],[0,0,3*b1(1)*abs(triggert1)/(triggert2-triggert1),3*b1(1)/b2(1),3,3]-16,'k-','LineWidth',lwd/1.5); % disinhibition
% ylim([-16.5,27]);
% yticks([0]);
% yticklabels({'0'});
% ylabel('           Activity (a.u.)');
% xticks([]);
% drawaxis(gca, 'x', 0, 'movelabel', 1);
% xlim([1, (prepresentt+predur+presentt+dur)/dt]);
% xlabel(' ');
% savefig(h, filename, outdir, fontsize, aspect15);