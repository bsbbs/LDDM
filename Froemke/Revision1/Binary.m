%% Binary choices
%% Common parameters
N = 2;
w = w0*ones(N);
a = a0*eye(N);
b = b0*eye(N);
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
mygray = [.62, .62, .62; 0, 0, 0];
mygreen = [.62, 1, .62; 0, 1, 0];
myred = [1, .62, .62; 1, 0, 0];
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
    winput0 = rand(1,N)*0;
    wrg0 = rand(N,N)*0;
    wgr0 = rand(1,N)*0;
    sgmR = 0; %12;
    sgmG = 0; %12;
    sgmInput = 0; %66;
    [R, G, D, Vcourse, winputc, wrgc, wgrc,ac] = LDDM_RndInputPlastic(Vprior, Vinput, BR, BG, winput0, wrg0, wgr0, a, b,...
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
filename = sprintf('MeanFR_N2');
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
filename = sprintf('WDynamic_N2');
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

%% Plastic only on wInput
N = 1;
Vprior = 0;
Vinput = 100;
BR = 30;
BG = 30;
dur = 1500; %s
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
R0 = BR*ones(N,1);
D0 = 0*ones(N,1);
G0 = BG*ones(N,1);
initialvals = [R0'; G0'; D0'];
winput0 = rand(1,N)*0;
wrg0 = rand(N,N)*0;
wgr0 = rand(1,N)*0;
a = a0*eye(N);
b = b0*eye(N);
sgmR = 0; %12;
sgmG = 0; %12;
sgmInput = 0; %66;
[R, G, D, Vcourse, winputo, wrgo, wgro] = LDDM_RndInputPlastico(Vprior, Vinput, BR, BG, winput0, wrg0, wgr0, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
winputo(end) % 0.8612
wrgo(end) % 0.8791
wgro(end) % 0
%%
h = figure; 
filename = sprintf('MeanFR_N%io', N);
rate = 1000;
x1 = 1:rate:(6*60/dt);
x2 = (24*60/dt):rate:(25*60/dt);
x = [x1, x2];
hold on;
ldg = [];
ldg(1) = plot(x1*dt,R(x1), 'g-', 'LineWidth',lwd);
ldg(2) = plot(x1*dt,G(x1), 'r-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,R(x2), 'g-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,G(x2), 'r-', 'LineWidth',lwd);
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ylim([25,140]);
yticks([40:40:120]);
yticklabels({'40','80','12'});
legend(ldg, {'R_1','G_1'}, "Location","east");
xlabel('Time (mins)');
ylabel('Firing rates (Hz)');
mysavefig(h, filename, plotdir, fontsize, [1.45,1.15], .75);

%%
h = figure; 
filename = sprintf('WDynamic_N%io', N);
hold on;
ldg = [];
yyaxis left;
ldg(1) = plot(x1*dt, winputo(x1), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
ldg(2) = plot(x1*dt, wrgo(x1), 'g-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,winputo(x2), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wrgo(x2), 'g-', 'LineWidth',lwd);
ylabel('Excitatory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ylim([-.05, 1]);
ax = gca;
ax.YColor = (1-.851)*[1,1,1];
yyaxis right;
ldg(3) = plot(x1*dt, wgro(x1), 'r-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wgro(x2), 'r-', 'LineWidth',lwd);
ylabel('Inhibitory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ylim([-.05, 1]);
ax.YColor = (1-.851)*[1,1,1];
legend(ldg, {'w_{Input}','w_{R to G}', 'w_{G to R}'}, "Location","east");
xlabel('Time (mins)');
mysavefig(h, filename, plotdir, fontsize, [1.74,1.15], .75);

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

%% Example trace for no inhibition network
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
G0 = wrg.*R0;

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
%% Example trace for inhibition network
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
G0 = wrg.*R0;

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
cp = [.032, .128, .256, .512];
vi = 100*(1+cp);
vo = 100*(1-cp);
h = figure;
filename = 'InptMtx_Choice';
hold on;
n = length(vi);
grays = flip(gray(numel(cp) + 2));
grays(1,:) = [];
grays(end,:) = [];
for i = 1:n
    plot(vo(i), vi(i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', grays(i,:), 'MarkerEdgeColor','none');
end
ylim([0, 200]);
xlim([0, 200]);
xlabel('V_2 (a.u.)');
ylabel('V_1 (a.u.)');
axis equal;
mysavefig(h, filename, plotdir, fontsize, [1.2917, 1.1667]*1.1, .75); 
%% Accuracy and reaction time over training 
task = sprintf('PerformtCourse_N%i',N);
% cp = [.032, .064,.128, .256];
cp = [0 2 4 8 16 32 64 128 256 512]'/1000;
b = eye(2);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = 0;
dur = 10; % second
stoprule = 1;
sims = 102400;
sgmG = 18;
myxrng = 10;
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
            R0 = eqlb*ones(N,1);
            D0 = 0*ones(N,1);
            G0 = wrg*R0+BG;
            initialvals = [R0'; G0'; D0'];
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
    filename = sprintf('%s%i', task, j);
    subplot(1,2,1); hold on;
    plot(x1*dt, RTs(1:length(x1),6), 'k-', 'LineWidth',lwd);
    plot(50+(x2-min(x2)+max(x1))*dt, RTs(length(x1)+1:end,6), 'k-', 'LineWidth',lwd);
    ylabel('Reaction time (s)');
    xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
    xticklabels({'0','2','4','6','24','25'});
    xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
    xlabel('Time (mins)');
    subplot(1,2,2); hold on;
    plot(x1*dt,ACCs(1:length(x1),6), 'k-', 'LineWidth',lwd);
    plot(50+(x2-min(x2)+max(x1))*dt,ACCs(length(x1)+1:end,6), 'k-', 'LineWidth',lwd);
    ylabel('Accuracy (%)');
    xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
    xticklabels({'0','2','4','6','24','25'});
    xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
    xlabel('Time (mins)');
    mysavefig(h, filename, plotdir, fontsize, [4.8, 2.1]);

%     h = figure; hold on;
%     filename = sprintf('Perform2D%i', j);
%     plot(ACCs(1:length(x1),1), RTs(1:length(x1),1), 'k-', 'LineWidth',lwd);
%     plot(ACCs(length(x1)+1:end,1), RTs(length(x1)+1:end,1), 'k-', 'LineWidth',lwd);
%     xlabel('Accuracy (%)');
%     ylabel('Reaction time (s)');

    h = figure;
    filename = sprintf('%s%i', task, j);
    subplot(1,2,1); hold on;
    smpl = [1, 3, 5, 10, 43];
    ldg = [];
    for ti = 1:length(smpl)
        plot(cp, RTs(smpl(ti),:), '-', 'Color', mygray(ti+1,:), 'LineWidth',lwd);
    end
    %legend(ldg, {'3.2 %', '6.4 %', '12.8 %', '25.6 %'}, "Location","east");
    ylabel('Reaction time (s)');
    set(gca,'XScale','log');
    xlabel('Evidence Coherence');
    mysavefig(h, filename, plotdir, fontsize, [4.8, 2.1]);

    subplot(1,2,2); hold on;
    for ti = 1:length(smpl)
        plot(cp,ACCs(smpl(ti),:), '-', 'Color', mygray(ti+1,:), 'LineWidth',lwd);
    end
    ylabel('Accuracy (%)');
    set(gca,'XScale','log');
    xlabel('Evidence Coherence');
    mysavefig(h, filename, plotdir, fontsize, [4.8, 2.1]);
    
end

%%
