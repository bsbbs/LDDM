%% Binary choices
%% Common parameters
N = 3;
w = w0*ones(N);
a = a0*eye(N);
b = b0*eye(N);
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
mygray = [.62, .62, .62; 0, 0, 0];
mygreen = [.62, 1, .62; 0, 1, 0];
myred = [1, .62, .62; 1, 0, 0];
myblue = [.62, .62, 1; 0, 0, 1];
myalpha = [0.7, 0];
cp = 1;
stoprule = 1;

%% Rate-based Plasticity 
filename = sprintf("PlasticSimulate_N%i.mat",N);
simrslt = fullfile(Simdir,filename);
if ~exist(simrslt,'file')
    Vprior = 0*ones(1,N);
    Vinput = 100*ones(1,N);
    BR = 30;
    BG = 30;
    dur = 1500; %s
    predur = 0;
    presentt = 0;
    stimdur = Inf;
    triggert = Inf;
    R0 = BR*ones(1,N);
    D0 = 0*ones(1,N);
    G0 = BG*ones(1,N);
    initialvals = [R0; G0; D0];
    winput0 = ones(1,N);
    wrg0 = rand(N,N)*0;
    wgr0 = rand(1,N)*0;
    sgmR = 0; %12;
    sgmG = 0; %12;
    sgmInput = 0; %66;
    [R, G, D, Vcourse, winputc, wrgc, wgrc,ac] = LDDM_RndInputPlasticSfN25(Vprior, Vinput, BR, BG, winput0, wrg0, wgr0, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    winputc(end,:) % 3.1250
    wrgc(end,:,:) % 0.5106, .7352; .6357
    wgrc(end,:) % 0.1548, .2229; .1927
    ac(end,:,:) % 0
    save(simrslt, "R","G","winputc","wrgc","wgrc","ac");
else
    load(simrslt);
end
%%
h = figure; 
filename = sprintf('MeanFR_N%i',N);
rate = 1000;
x1 = 1:rate:(6*60/dt);
x2 = (24*60/dt):rate:(25*60/dt);
x = [x1, x2];
hold on;
ldg = [];
ldg(1) = plot(x1*dt,R(x1,1), 'g-', 'LineWidth',lwd);
% ldg(2) = plot(x1*dt,R(x1,2), 'g--', 'LineWidth',lwd);
ldg(2) = plot(x1*dt,G(x1,1), 'r-', 'LineWidth',lwd);
% ldg(4) = plot(x1*dt,G(x1,2), 'r--', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,R(x2,1), 'g-', 'LineWidth',lwd);
% plot(50+(x2-min(x2)+max(x1))*dt,R(x2,2), 'g--', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,G(x2,1), 'r-', 'LineWidth',lwd);
% plot(50+(x2-min(x2)+max(x1))*dt,G(x2,2), 'r--', 'LineWidth',lwd);
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ylim([28,65]);
legend(ldg, {'R_1','G_1'}, "Location","east");
xlabel('Time (mins)');
ylabel('Firing rates (Hz)');
mysavefig(h, filename, plotdir, fontsize, [1.45,1.15], .75); % 70/128.77

%%
h = figure; 
filename = sprintf('WDynamic_N%i',N);
hold on;
ldg = [];
yyaxis left;
ldg(1) = plot(x1*dt, winputc(x1,1), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
% ldg(2) = plot(x1*dt, winputc(x1,2), '--', 'Color', [.38,.38,.38], 'LineWidth',lwd);
ldg(2) = plot(x1*dt, wrgc(x1,1,1), 'g-', 'LineWidth',lwd);
% ldg(4) = plot(x1*dt, wrgc(x1,1,2), 'g--', 'LineWidth',lwd);
%ldg(5) = plot(x1*dt, wrgc(x1,2,1), 'g-', 'LineWidth',lwd);
% ldg(6) = plot(x1*dt, wrgc(x1,2,2), 'g--', 'LineWidth',lwd);
% ldg(3) = plot(x1*dt, ac(x1,1,1), 'b-', 'LineWidth',lwd);
% ldg(10) = plot(x1*dt, ac(x1,2,2), 'b--', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,winputc(x2,1), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wrgc(x2,1,1), 'g-', 'LineWidth',lwd);
% plot(50+(x2-min(x2)+max(x1))*dt,ac(x2,1,1), 'b-', 'LineWidth',lwd);
ylabel('Excitatory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ax = gca;
ax.YColor = (1-.851)*[1,1,1];
yyaxis right;
ldg(3) = plot(x1*dt, wgrc(x1,1), 'r-', 'LineWidth',lwd);
% ldg(8) = plot(x1*dt, wgrc(x1,2), 'r--', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wgrc(x2,1), 'r-', 'LineWidth',lwd);
ylabel('Inhibitory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ax.YColor = (1-.851)*[1,1,1];
legend(ldg, {'w_{Input}','w_{R to G}', 'w_{G to R}'}, "Location","east");
xlabel('Time (mins)');
mysavefig(h, filename, plotdir, fontsize, [1.7,1.15], .75); % 70/128.77


%% Input values testing representation 
% load(fullfile(out_dir, 'LIPMeanData_TrinaryConditions.mat'));
% vi = LIPMeanData_TrinaryConditions.Vi;
% vo = LIPMeanData_TrinaryConditions.Vo;
vi = [0:50:200];
vo = [0:50:300, 500];
h = figure;
filename = 'InptMtx_Rprsnt';
hold on;
vi = [0:50:200];
vo = 250;
n = length(vi);
greens = {'#D3FFD3', '#A6FFA6', '#4CD1A0', '#26BA8A', '#009964'};
for i = 1:n
    plot(vo, vi(i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', greens{i}, 'MarkerEdgeColor','none');
end
vi = 100;
vo = [0:50:300, 500];
n = length(vo);
reds = {'#F4D9DD','#FFCCCC', '#E9B4BA','#DE8E97','#D36974','#C84451','#BD1E2E','#8B0000'};
for i = 1:n
    plot(vo(i), vi, 'o', 'MarkerSize', 6, 'MarkerFaceColor', reds{i}, 'MarkerEdgeColor','none');
end
vi = [0:50:200];
vo = [0:50:300, 500];
[vi, vo] = meshgrid(vi, vo);
plot(vo, vi, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','k');
ylim([-5, 200]);
xlim([-5, 500]);
xlabel('V_2 (a.u.)');
ylabel('V_1 (a.u.)');
axis equal;
mysavefig(h, filename, plotdir, fontsize, [1.2917, 1.1667]*1.1, .75);

%% Example trace for no feedback inhibition network
N = 2;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
task = sprintf('Trace_N2o');
h = figure;
filename = sprintf('%s',task);
sgmR = 0;
sgmG = 0;
sgmInput = 0;

winput = mean(winputo(end,:))*ones(1,N);
wrg = mean(wrgo(end,:))*ones(N,N);
wgr = mean(wgro(end,:))*ones(1,N);
a = eye(N)*0;
b = eye(N)*b0;

R0 = 20*ones(1,N);
D0 = 0*ones(1,N);
G0 = (wrg*R0')';

initialvals = [R0; G0; D0];
Vprior = zeros(1,N);
subplot(2,1,1);hold on;
ldg = zeros(1,5);
vi = [0:50:200];
vo = 250;
for ci = 1:length(vi)
    Vinput = [vi(ci), vo];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(ci) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', greens{ci}, 'LineWidth',lwd);
end
xlim([-50*dt, dur]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.8]);

subplot(2,1,2);hold on;
ldg = zeros(1,8);
vi = 100;
vo = [0:50:300, 500];
for ci = 1:length(vo)
    Vinput = [vi, vo(ci)];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(ci) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', reds{ci}, 'LineWidth',lwd);
end
xlim([-50*dt, dur]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.8]);

task = 'Rprsnt2o';
h = figure;
filename = sprintf('%s',task);
hold on;
vi = [0:50:200];
vo = [0:50:300, 500];
[vi, vo] = meshgrid(vi, vo);
R1mat = vi;
R2mat = vo;
dur = 10;
for i = 1:size(vi,1)
    for j = 1:size(vi, 2)
        Vinput = [vi(i,j), vo(i,j)];
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        R1mat(i, j) = mean(R(5000:end,1));
        R2mat(i, j) = mean(R(5000:end,2));
    end
end
for j = 1:size(vi,2)
    plot(vo(:,j), R1mat(:,j), '-', 'Color', greens{j}, 'LineWidth',lwd);
end
for i = 1:size(vi,1)
    plot(vo(i,:), R1mat(i,:), 'o', 'MarkerSize', 4, 'MarkerFaceColor', reds{i}, 'MarkerEdgeColor','k');
end
xlim([-5, 500]);
ylabel('R_1 Firing rates (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.6]);

%% Example trace for feedback inhibition network
N = 2;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
task = sprintf('Trace_N2');
h = figure;
filename = sprintf('%s',task);
sgmR = 0;
sgmG = 0;
sgmInput = 0;

winput = mean(winputc(end,:))*ones(1,N);
wrg = mean(wrgc(end,:))*ones(N,N);
wgr = mean(wgrc(end,:))*ones(1,N);
a = eye(N)*0;
b = eye(N)*b0;

R0 = 20*ones(1,N);
D0 = 0*ones(1,N);
G0 = (wrg*R0')';

initialvals = [R0; G0; D0];
Vprior = zeros(1,N);
subplot(2,1,1);hold on;
ldg = zeros(1,5);
vi = [0:50:200];
vo = 250;
for ci = 1:length(vi)
    Vinput = [vi(ci), vo];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(ci) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', greens{ci}, 'LineWidth',lwd);
end
xlim([-50*dt, dur]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.8]);

subplot(2,1,2);hold on;
ldg = zeros(1,8);
vi = 100;
vo = [0:50:300, 500];
for ci = 1:length(vo)
    Vinput = [vi, vo(ci)];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(ci) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', reds{ci}, 'LineWidth',lwd);
end
xlim([-50*dt, dur]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.8]);

task = 'RprsntN2';
h = figure;
filename = sprintf('%s',task);
hold on;
vi = [0:50:200];
vo = [0:50:300, 500];
[vi, vo] = meshgrid(vi, vo);
R1mat = vi;
R2mat = vo;
dur = 10;
for i = 1:size(vi,1)
    for j = 1:size(vi, 2)
        Vinput = [vi(i,j), vo(i,j)];
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        R1mat(i, j) = mean(R(5000:end,1));
        R2mat(i, j) = mean(R(5000:end,2));
    end
end
for j = 1:size(vi,2)
    plot(vo(:,j), R1mat(:,j), '-', 'Color', greens{j}, 'LineWidth',lwd);
end
for i = 1:size(vi,1)
    plot(vo(i,:), R1mat(i,:), 'o', 'MarkerSize', 4, 'MarkerFaceColor', reds{i}, 'MarkerEdgeColor','k');
end
xlim([-5, 500]);
ylabel('R_1 Firing rates (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdir, fontsize, [1.6,1.6]);

%% Input values testing choice 
cp = [3.2 12.8, 25.6, 38.4 51.2]'/100; % percentage of coherence
vi = 100*(1+cp);
vo = 100*(1-cp);
h = figure;
filename = 'InptMtx_Choice';
hold on;
n = length(vi);
grays = flip(gray(numel(cp) + 1));
grays(1,:) = [];
% grays(end,:) = [];
for i = 1:n
    plot(vo(i), vi(i), 'o', 'MarkerSize', 4, 'MarkerFaceColor', grays(i,:), 'MarkerEdgeColor','k');
end
ylim([0, 200]);
xlim([0, 200]);
xlabel('V_2 (a.u.)');
ylabel('V_1 (a.u.)');
% axis equal;
mysavefig(h, filename, plotdir, fontsize, [1.2917, 1.1917]*1.1, .75); % [1.2917, 1.1667]*1.1

%% Example Trace for choice, without feedback inhibition
N = 2;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = dt;
dur = 1; % second
thresh = 700;
task = sprintf('ChoiceDynamics_N2o');
h = figure;
filename = sprintf('%s',task);
sgmR = 0;
sgmG = 0;
sgmInput = 0;

winput = mean(winputo(end,:))*ones(1,N);
wrg = mean(wrgo(end,:))*ones(N,N);
wgr = mean(wgro(end,:))*ones(1,N);
a = eye(N)*0;
b = eye(N)*b0;

R0 = 20*ones(1,N);
D0 = 0*ones(1,N);
G0 = (wrg*R0')';

initialvals = [R0; G0; D0];
Vprior = zeros(1,N);
hold on;
n = length(vi);
ldg = zeros(1,6);
i = 0;
for ci = 1:2:length(vi)
    i = i + 1;
    Vinput = [vi(ci), vo(ci)];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(i) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', grays(ci,:), 'LineWidth',lwd);
    ldg(i+3) = plot(x(1:5:end), R(1:5:end,2), '--', 'Color', grays(ci,:), 'LineWidth',lwd);
end
plot([.04,.1],[70,70],'-','Color', mygray(1,:));
legend(ldg, {'','','','3.2','25.6','51.2'}, 'NumColumns', 2, 'Location', 'northwest');
xlim([-5*dt, .1]);
ylim([0, 80]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.4,1.4]);
%%
task = sprintf('Vctrfld_N2o');
filename = sprintf('%s',task);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = dt;
dur = 1; % second
thresh = 115;
sgmR = 0;
sgmG = 0;
sgmInput = 0;
BR = 30;
BG = 30;
V = 100*[1, 1];
h = figure;
hold on;
winput = mean(winputo(end,:))*ones(1,N);
wrg = mean(wrgo(end,:))*ones(N,N);
wgr = mean(wgro(end,:))*ones(1,N);
R0 = 21:14:160;
% - Nullclines R1*-R2* space
R2 = (winput(2)*V(2) + BR)*ones(size(R0));

lgd2 = plot(R2, R0, 'k--','LineWidth',lwd/2); % dR2/dt = 0
Line1 = [R0' R2'];
R1 = (winput(1)*V(1) + BR)*ones(size(R0));
lgd1 = plot(R0, R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
Line2 = [R1' R0'];

[R1, R2] = meshgrid(R0, R0);
dR1 = -R1 + (winput(1).*V(1) + BR);
dR2 = -R2 + (winput(2).*V(2) + BR);
rate = 0.15;
dR1 = dR1.*rate;
dR2 = dR2.*rate;
quiver(R2,R1, dR2, dR1,'off','Color', mygray(1,:),'LineWidth',1);
% - trajectories of R1 & R2
if 1
    R0 = 20*ones(1,N);
    D0 = 0*R0;
    G0 = (wrg*R0')';

    initialvals = [R0; G0; D0];
    Vprior = zeros(1,N);
    ldgtrc = [];
    rng(42);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, V, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput+10, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, 1);
    ldgtrc(1) = plot(R(1:50:end,2), R(1:50:end,1),'-','Color','#ef476f','LineWidth',lwd/3);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, V, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput+10, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, 1);
    ldgtrc(2) = plot(R(1:50:end,2), R(1:50:end,1),'-','Color','#118ab2','LineWidth',lwd/3);
end
xlim([18, 160]);
ylim([18, 160]);
xlabel('R_2 firing rate (Hz)');
ylabel('R_1 firing rate (Hz)');
mysavefig(h, filename, plotdir, fontsize, [1.8,1.6]*1.12, 0);

%% Example Trace for choice, with feedback inhibition
N = 2;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = dt;
dur = 1; % second
thresh = 500;
task = sprintf('ChoiceDynamics_N2');
h = figure;
filename = sprintf('%s',task);
sgmR = 0;
sgmG = 0;
sgmInput = 0;

winput = mean(winputc(end,:))*ones(1,N);
wrg = mean(wrgc(end,:))*ones(N,N);
wgr = mean(wgrc(end,:))*ones(1,N);
a = eye(N)*0;
b = eye(N)*b0;

R0 = 20*ones(1,N);
D0 = 0*R0;
G0 = (wrg*R0')';

initialvals = [R0; G0; D0];
Vprior = zeros(1,N);
hold on;
n = length(vi);
ldg = zeros(1,6);
i = 0;
for ci = 1:2:length(vi)
    i = i + 1;
    Vinput = [vi(ci), vo(ci)];
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    ldg(i) = plot(x(1:5:end), R(1:5:end,1), '-', 'Color', grays(ci,:), 'LineWidth',lwd);
    ldg(i+3) = plot(x(1:5:end), R(1:5:end,2), '--', 'Color', grays(ci,:), 'LineWidth',lwd);
end
plot([.2,1],[70,70],'-','Color', mygray(1,:));
% legend(ldg, {'','','','','','3.2','12.8','25.6','38.4','51.2'}, 'NumColumns', 2, 'Location', 'northwest');
xlim([-50*dt, 1]);
ylim([0, 80]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [1.4,1.4]);
%%
h = figure;
task = sprintf('Vctrfld_N2');
filename = sprintf('%s',task);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = dt;
dur = 10; % second
thresh = 500;
sgmR = 0;
sgmG = 0;
sgmInput = 0;
BR = 30;
BG = 30;
V = 100*[1, 1];
h = figure;
hold on;
winput = mean(winputc(end,:))*ones(1,N);
wrg = mean(wrgc(end,:))*ones(N,N);
wgr = mean(wgrc(end,:))*ones(1,N);
a = eye(N)*0;
b = eye(N)*b0;
% - Nullclines R1*-R2* space
R2 = 10.^(linspace(-1,2.5,300));
R1 = ((winput(2)*V(2) + a(2,2)*R2+BR)./R2 - 1 - wgr(2)*((wrg(2,2) - b(2,2))*R2 + BG))/(wgr(2)*wrg(2,1)); % dR2/dt = 0
lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
Line1 = [R1' R2'];
R1 = 10.^(linspace(-1,2.5,300));
R2 = ((winput(1)*V(1) + a(1,1)*R1+BR)./R1 - 1 - wgr(1)*((wrg(1,1) - b(1,1))*R1 + BG))/(wgr(1)*wrg(1,2)); % dR1/dt = 0
lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
Line2 = [R1' R2'];
% set(gca, 'XScale','log');
% set(gca, 'YScale','log');
% - trajectories of R1 & R2
if 1
    R0 = 20*ones(1,N);
    D0 = 0*R0;
    G0 = (wrg*R0')';

    initialvals = [R0; G0; D0];
    Vprior = zeros(1,N);
    ldgtrc = [];
    rng(42);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, V, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput+5, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, 1);
    ldgtrc(1) = plot(R(1:50:end,2), R(1:50:end,1),'-','Color','#ef476f','LineWidth',lwd/3);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, V, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput+5, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, 1);
    ldgtrc(2) = plot(R(1:50:end,2), R(1:50:end,1),'-','Color','#118ab2','LineWidth',lwd/3);
end
xlabel('R_2 firing rate (Hz)');
ylabel('R_1 firing rate (Hz)');
R0 = 21:8:100;
[R1, R2] = meshgrid(R0, R0);
dR1 = -R1 + (winput(1).*V(1) + BR)./(1+wgr(1)*(wrg(1,1)*R1+wrg(1,2)*R2 - b(1,1)*R1 + BG));
dR2 = -R2 + (winput(2).*V(2) + BR)./(1+wgr(2)*(wrg(2,2)*R2+wrg(2,1)*R1 - b(2,2)*R2 + BG));
rate = 0.15;
dR1 = dR1.*rate;
dR2 = dR2.*rate;
quiver(R2, R1, dR2, dR1, 'off','Color', mygray(1,:),'LineWidth',1);
xlim([18, 120]);
ylim([18, 120]);
mysavefig(h, filename, plotdir, fontsize, [1.8,1.6]*1.12, 0);

%% Representation precision over before and after training
N = 2;
cp = 6.4/100;
b = eye(2)*b0;
a = eye(2)*a0;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 120*2; % second
stoprule = 1;
sgmG = 18;
task = 'Rprsnt_N2';
names = ["Early noise", "Late noise"];
for testi = 1:2
    switch testi
        case 1
            sgmInput = 66;
            sgmR = 0;
        case 2
            sgmInput = 0;
            sgmR = 26;
    end
    h = figure;
    filename = sprintf('%s_%s',task,names{testi});
    % begin...
    rng(2025);
    for i = [1,2]
        if i == 1
            winput = mean(winputo(end,:))*ones(1,N);
            wrg = mean(wrgo(end,:))*ones(N,N);
            wgr = mean(wgro(end,:))*ones(1,N);
        elseif i == 2
            winput = mean(winputc(end,:))*ones(1,N);
            wrg = mean(wrgc(end,:))*ones(N,N);
            wgr = mean(wgrc(end,:))*ones(1,N);
        end
        R0 = 20*ones(1,N);
        D0 = 0*R0;
        G0 = (wrg*R0')';
        initialvals = [R0; G0; D0];
        Vprior = zeros(1,N);
        Vinput = 100*[1+cp, 1-cp];
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        subplot(1,2,i);
        hold on;
        mylgd = [];
        [~, x] = histcounts(R(5000:end,1));
        pd1 = fitdist(R(5000:end,1),'kernel','Kernel','normal');
        y1 = pdf(pd1,x);
        mylgd(1) = fill([x, flip(x)], [y1, zeros(size(x))], [.9,.9,1], 'EdgeColor', colorpalette{end}, 'FaceAlpha', .4);
        [~, x] = histcounts(R(5000:end,2));
        pd2 = fitdist(R(5000:end,2),'kernel','Kernel','normal');
        y2 = pdf(pd2,x);
        mylgd(2) = fill([x, flip(x)], [y2, zeros(size(x))], [.9,.9,1], 'EdgeColor', colorpalette{end}, 'FaceAlpha', .4);
        ylabel('Density');
        xlabel('Activity (Hz)');
        mysavefig(h, filename, plotdir, fontsize, [4.4, 1.2]);
    end
end
%% SNR over training
%% Signal-to-Noise ratio over training 
N = 2;
task = sprintf('SNRtCourse_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 120; % second
stoprule = 1;
sgmG = 18;
names = ["Early noise", "Late noise"];
netwks = ["E", "E-I"];
myy1rng = [120,38];
for testi = 1:2
    switch testi
        case 1
            sgmInput = 66;
            sgmR = 0;
        case 2
            sgmInput = 0;
            sgmR = 26;
    end
    rate = 10000;
    x1 = 1:rate:(6*60/dt);
    x2 = (24*60/dt):rate:(25*60/dt);
    tvec = [x1, x2];
    for i = [1, 2]
        R0 = 20*ones(1,N);
        D0 = 0*R0;
        G0 = (wrg*R0')';
        initialvals = [R0; G0; D0];
        Vprior = zeros(1,N);
        Vinput = 100*[1+cp, 1-cp];
        simfile = fullfile(Simdir, sprintf('%s_%s_%i.mat',task, names{testi}, i));
        if ~exist(simfile, 'file')
            SNR = nan(length(tvec),2);
            SD = nan(length(tvec),3);
            r = nan(length(tvec), 1);
            MhnbsDstnc = nan(length(tvec), 1); % Mahalanobis distance
            for ti = 1:length(tvec)
                winput = winputc(tvec(ti));
                wrg = wrgc(tvec(ti));
                wgr = wgrc(tvec(ti));
                if i == 1
                    winput = mean(winputo(tvec(ti),:))*ones(1,N);
                    wrg = mean(wrgo(tvec(ti),:))*ones(N,N);
                    wgr = mean(wgro(tvec(ti),:))*ones(1,N);
                elseif i == 2
                    winput = mean(winputc(tvec(ti),:))*ones(1,N);
                    wrg = mean(wrgc(tvec(ti),:))*ones(N,N);
                    wgr = mean(wgrc(tvec(ti),:))*ones(1,N);
                end
                data = [];
                [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
                    sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
                data(:,1) = R(5000:end,1);
                data(:,2) = R(5000:end,2);
                mu = mean(data);
                Sigma = cov(data);
                SD(ti,:) = [Sigma(1,1), Sigma(2,2), Sigma(1,2)];
                SNR(ti,:) = mu./[Sigma(1,1), Sigma(2,2)];
                r(ti) = Sigma(1,2)/sqrt(Sigma(1,1)*Sigma(2,2));
                %MhnbsDstnc(ti) = sqrt((mu(1) - mu(2))/Sigma*(mu(1) - mu(2))');
            end
            save(simfile,'mu','Sigma','SD','SNR','r');
        else
            load(simfile);
        end
        mu1 = SNR(:,1).*SD(:,1);
        mu2 = SNR(:,2).*SD(:,2);
        dp = abs(mu1 -  mu2)./sqrt(SD(:,1) + SD(:,2) - 2*SD(:,3));

        h = figure;
        set(h, 'Units', 'inches','Position', [1,1,2.2,1.2]);
        filename = sprintf('%s_%s_%i',task, names{testi}, i);
        hold on;
        ldg = [];
        yyaxis left;
        ldg(1) = plot(x1*dt,SD(1:length(x1),2), '-', 'Color', '#ef476f', 'LineWidth',lwd);
        ldg(2) = plot(x1*dt,SD(1:length(x1),1), '-', 'Color', '#118ab2',  'LineWidth',lwd);
        plot(50+(x2-min(x2)+max(x1))*dt,SD(length(x1)+1:end,2), '-', 'Color', '#ef476f', 'LineWidth',lwd);
        plot(50+(x2-min(x2)+max(x1))*dt,SD(length(x1)+1:end,1), '-', 'Color', '#118ab2', 'LineWidth',lwd);
        ylabel('Variance');
        ylim([0, myy1rng(testi)]);
        ax = gca;
        ax.YColor = (1-.851)*[1,1,1];
        yyaxis right;
        ldg(3) = plot(x1*dt,dp(1:length(x1)), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
        plot(50+(x2-min(x2)+max(x1))*dt,dp(length(x1)+1:end), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
        ylabel('d''');
        xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
        xticklabels({'0','2','4','6','24','25'});
        xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
        ax.YColor = (1-.851)*[1,1,1];
        xlabel('Time (mins)');
        ylim([0, 2]);
        if testi == 1 && i == 1
            legend(ldg, {'R_1', 'R_2', 'd'''}, "Location","east");
        end
        % title(sprintf('%s %s', names(testi), netwks(i)));
        mysavefig(h, filename, plotdir, fontsize, [2.2,1.2]);
        % exportgraphics(h, fullfile(plotdir, [filename, '.pdf']), 'ContentType','vector');
    end
end
%% Accuracy and reaction time over training 
cp = [0 2 4 8 16 32 64 128 256 512]'/1000;
b = eye(2)*b0;
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = 0;
dur = 10; % second
stoprule = 1;
sims = 102400;
sgmG = 18;
clear Vprior Vinput;
for j = 1:2
    switch j
        case 1
            sgmInput = 66;
            sgmR = 0;
            myylim = [0,60];
            myyrng = 30;
        case 2
            sgmInput = 0;
            sgmR = 26;
            myylim = [20,40];
            myyrng = 18;
    end
    rate = 10000;
    x1 = 1:rate:(6*60/dt);
    x2 = (24*60/dt):rate:(25*60/dt);
    tvec = [x1, x2];
    filename = sprintf("PerformSimulate_N%i_%i.mat",N,j);
    simrslt = fullfile(Simdir,filename);
    if ~exist(simrslt, 'file')
        RTs = nan(length(tvec),length(cp));
        ACCs = nan(length(tvec),length(cp));
        for ti = 1:length(tvec)
            winput = winputc(tvec(ti),:);
            wrg = squeeze(wrgc(tvec(ti),:,:));
            wgr = wgrc(tvec(ti),:);
            a = squeeze(ac(tvec(ti),:,:));
%             R0 = eqlb*ones(N,1);
%             D0 = 0*ones(N,1);
%             G0 = wrg*R0+BG;

            R0 = 20*ones(1,N);
            D0 = 0*ones(1,N);
            G0 = (wrg*R0')';

            initialvals = [R0; G0; D0];
            Vprior.V1 = 100*ones(size(cp))';
            Vprior.V2 = 100*ones(size(cp))';
            V1 = 100*(1+cp)';
            V2 = 100*(1-cp)';
            Vinput.V1 = V1;
            Vinput.V2 = V2;
            [rt, choice, argmaxR] = LDDM_Rndinput_STDP_GPUrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
                sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            RTs(ti,:) = mean(squeeze(rt),2)';
            ACCs(ti,:) = mean(2-squeeze(choice),2)'*100;
        end
        save(simrslt, "RTs","ACCs");
    else
        load(simrslt);
    end
    
    h = figure;
    filename = sprintf('PerfromancetCourse%i', j);
    smpl = [6, 8, 9];
    tsmpl = [1, 3, 5, 10, 43];
    mygray = flip(gray(numel(smpl)+1));
    subplot(1,2,1); hold on;
    for i = 1:numel(smpl)
        plot(x1*dt, RTs(1:length(x1),smpl(i)), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
        plot(50+(x2-min(x2)+max(x1))*dt, RTs(length(x1)+1:end,smpl(i)), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
    end
    for ti = 1:numel(tsmpl)
        x = tvec(tsmpl(ti))*ones(1,2)*dt;
        if ti == 5
            x = (25*60+(-min(x2)+max(x1))*dt+50)*ones(1,2);
        end
        plot(x, [.3, .2], '-', 'Color', colorpalette{ti});
    end
    ylabel('Reaction time (s)');
    xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
    xticklabels({'0','2','4','6','24','25'});
    xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
    ylim([.2, 1]);
    xlabel('Time (mins)');
    mysavefig(h, filename, plotdir, fontsize, [4.8, 1.8]);

    subplot(1,2,2); hold on;
    ldg = [];
    for i = 1:numel(smpl)
        ldg(i) = plot(x1*dt,ACCs(1:length(x1),smpl(i)), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
        plot(50+(x2-min(x2)+max(x1))*dt,ACCs(length(x1)+1:end,smpl(i)), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
    end
    for ti = 1:numel(tsmpl)
        x = tvec(tsmpl(ti))*ones(1,2)*dt;
        if ti == 5
            x = (25*60+(-min(x2)+max(x1))*dt+50)*ones(1,2);
        end
        plot(x, [58, 50], '-', 'Color', colorpalette{ti});
    end
    ylabel('Accuracy (%)');
    xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
    xticklabels({'0','2','4','6','24','25'});
    xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
    xlabel('Time (mins)');
    ylim([50,100]);
    legend(ldg, {'3.2','12.8','25.6'});
    mysavefig(h, filename, plotdir, fontsize, [4.8, 1.8]);

    h = figure;
    filename = sprintf('PerformancePsychMetric%i', j);
    subplot(1,2,1); hold on;
    tsmpl = [1, 3, 5, 10, 43];
    mygray = flip(gray(numel(tsmpl) + 1));
    for ti = 1:length(tsmpl)
        plot(cp*100, RTs(tsmpl(ti),:), '-', 'Color', colorpalette{ti}, 'LineWidth',lwd);
    end
    ylim([0.2,1]);
    ylabel('Reaction time (s)');
    set(gca,'XScale','log');
    xlabel('Input coherence (c'')(%)');
    mysavefig(h, filename, plotdir, fontsize, [4.4, 1.8]);

    subplot(1,2,2); hold on;
    ldg = [];
    for ti = 1:length(tsmpl)
        ldg(ti) = plot(cp*100, ACCs(tsmpl(ti),:), '-', 'Color', colorpalette{ti}, 'LineWidth',lwd);
    end
     legend(flip(ldg), flip({'0 secs', '20 secs', '40 secs', '1.5 min', '25 mins'}), "Location","east");
    ylabel('Accuracy (%)');
    set(gca,'XScale','log');
    xlabel('Input coherence (c'')(%)');
    mysavefig(h, filename, plotdir, fontsize, [4.4, 1.8]);
end
%% Performance visualized in 2D spaces
rate = 10000;
x1 = 1:rate:(6*60/dt);
x2 = (24*60/dt):rate:(25*60/dt);
tvec = [x1, x2];
RTBoth = zeros([numel(tvec), numel(cp), 2]);
ACCBoth = zeros([numel(tvec), numel(cp), 2]);
for j = 1:2
    switch j
        case 1
            sgmInput = 66;
            sgmR = 0;
            myylim = [0,60];
            myyrng = 30;
        case 2
            sgmInput = 0;
            sgmR = 26;
            myylim = [20,40];
            myyrng = 18;
    end

    simfile = sprintf("PerformSimulate_N%i_%i.mat",N,j);
    simrslt = fullfile(Simdir,simfile);
    load(simrslt);
    RTBoth(:,:,j) = RTs;
    ACCBoth(:,:,j) = ACCs;
end
h = figure;
filename = sprintf('Perform2D');
smpl = [6, 8, 9];
tsmpl = [1, 3, 5, 10, 43];
mygray = flip(gray(numel(smpl)+1));
subplot(1,2,1); hold on;
for i = 1:numel(smpl)
    plot(RTBoth(1:length(x1), smpl(i), 2), RTBoth(1:length(x1), smpl(i), 1), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
    % plot(50+(x2-min(x2)+max(x1))*dt, RTs(length(x1)+1:end,smpl(i)), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
end
plot([.2,1],[.2,1],'--','Color',mygray(3,:));
xlim([.2,1]);
ylim([.2,1]);
xlabel('Reaction time (s) - Late Noise');
ylabel('Reaction time (s) - Early Noise');
mysavefig(h, filename, plotdir, fontsize, [4.4, 1.8]);
subplot(1,2,2); hold on;
for i = 1:numel(smpl)
    plot(ACCBoth(1:length(x1), smpl(i), 2), ACCBoth(1:length(x1), smpl(i), 1), '-', 'Color', mygray(i+1,:), 'LineWidth',lwd);
    % plot(ACCs(length(x1)+1:end,1), RTs(length(x1)+1:end,1), 'k-', 'LineWidth',lwd);
end
plot([50,100],[50,100],'--','Color',mygray(3,:));
xlim([50,100]);
ylim([50,100]);
xlabel('Accuracy (%) - Late Noise');
ylabel('Accuracy (%) - Early Noise');
mysavefig(h, filename, plotdir, fontsize, [4.4, 1.8]);

