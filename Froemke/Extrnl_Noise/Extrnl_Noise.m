%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(genpath(fullfile(Homedir,'Documents','LDDM','Froemke')));
Glgdir = 'G:\My Drive';
% Glgdir = '/Volumes/GoogleDrive/My Drive';
out_dir = fullfile(Glgdir, 'LDDM/Froemke/Extrnl_Noise');
if ~exist("out_dir",'dir')
    mkdir(out_dir);
end
cd(out_dir);
plotdir = fullfile('./Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end
%% parameters for visulization
fontsize = 14;
mksz = 25;
lwd = 1.5;
cp = [.032, .064,.128, .256];
mygray = flip(gray(numel(cp) + 2));
myred = mygray; myred(:,1) = 1;
myblue = mygray; myblue(:,3) = 1;
%% parameters for simulation
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauD = .1; % second
Tau = [tauR, tauG, tauD];
w0 = 1;
a0 = 15;
b0 = 1; % 1.1;
c_rprsnt = .064;
sgmInput_rprsnt = 1/3;
c_choice = .032;
sgmInput_choice = .75;
predur = 0;
presentt = dt;
thresh = 70; % Hz
stoprule = 1;
eqlb = 32;
scale = (2*w0 - b0)*eqlb^2 + (1-a0)*eqlb; %2*w0*eqlb.^2 + (1-a0).*eqlb;
R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0)))/2/(2*w0 - b0);
D0 = b0*R0;
G0 = (2*w0-b0)*R0;
Vprior = [1,1]*scale;
%% Representation dynamic
task = 'VR_RT';
a = a0*eye(2);
b = b0*eye(2);
w = ones(2);
stimdur = Inf;
triggert = Inf;
dur = 1.8; % second
% Dynamics 
sgm = .01;
sgmInput = sgmInput_rprsnt;
Vinput = [1 + c_rprsnt, 1 - c_rprsnt]*scale;
a = a0*eye(2);
b = b0*eye(2);
w = w0*ones(2);
initialvals = [1,1;2,2;0,0]*eqlb/3;
h = figure;
filename = sprintf('timeCourse_%s_sgm%2.2f_sgmInput%.2f',task, sgm, sgmInput);
% Value representation task
rng(5);
[choice, rt, R, G, D, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
subplot(2,1,1); hold on;
x = (1:numel(R(:,1)))*dt;
plot(x, R(:,1), 'r-', 'LineWidth', lwd);
plot(x, R(:,2), 'b-', 'LineWidth', lwd);
xlim([-50*dt, dur]);
% plot(x, G(:,1), 'r--', 'LineWidth', lwd);
% plot(x, G(:,2), 'b--', 'LineWidth', lwd);
% plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
% plot(x, D(:,2), 'b-.', 'LineWidth', lwd);
yticks([10, 15, 20, 25]);
ylabel('Firing rates (Hz)');
xlabel(' ');
legend({'R_1','R_2'},'Box','off','Location','southeast',...
    'FontName','Arial', 'FontAngle','italic',...
    'FontSize',fontsize-4);
savefigs(h, filename, plotdir, fontsize, [2.9,4]);

% Reaction-time task
stimdur = Inf;
triggert = presentt;
dur = 1.8; % second
rng(8);
R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0)))/2/(2*w0 - b0);
D0 = b0*R0;
G0 = (2*w0-b0)*R0;
initialvals = [R0,R0; G0,G0; D0,D0];
subplot(2,1,2); hold on;
rtvec = [];
for vi = 1:4
    Vinput = [1+cp(vi), 1-cp(vi)]*scale;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    rtvec(vi) = rt;
    x = (1:numel(R(:,1)))*dt;
    plot(x, R(:,1), '-', 'LineWidth', lwd, 'Color', myred(vi+1,:));
    plot(x, R(:,2), '-', 'LineWidth', lwd, 'Color', myblue(vi+1,:));
    % plot(x, G(:,1), 'r--', 'LineWidth', lwd);
    % plot(x, G(:,2), 'b--', 'LineWidth', lwd);
    % plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
    % plot(x, D(:,2), 'b-.', 'LineWidth', lwd);
end
plot(rtvec, thresh*ones(size(rtvec)),'k-','LineWidth',.5);
text(mean(rtvec), thresh*1.05,'Threshold','FontName','Arial');
xlim([-50*dt, dur]);
ylim([10,thresh]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
savefigs(h, filename, plotdir, fontsize, [2.9,4]);
%% General property 1. integration of noise
predur = 0;
dur = 100;
triggert = Inf;
stimdur = Inf;
presentt = dt;
sgm = .01;
sgmInput = sgmInput_rprsnt;
Vinput = [1 + c_rprsnt, 1 - c_rprsnt]*scale;
a = a0*eye(2);
b = b0*eye(2);
w = w0*ones(2);
initialvals = [1,1;2,2;0,0]*eqlb/3;
filename = sprintf('LDDM_timeCourse_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f',a0,b0,sgm,sgmInput);
output = fullfile(Simdir,[filename, '.mat']);
if ~exist(output, 'file')
    rng(4);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    save(output,'R','G','D','Vcourse');
else
    load(output);
end
%
h = figure;  hold on;
filename = 'Input';
h1 = histogram(Vcourse(2000:end,1),100);
h1.FaceColor = 'r'; %myred(3,:);
h1.EdgeColor = 'n';
h2 = histogram(Vcourse(2000:end,2),100);
h2.FaceColor = 'b';%myblue(3,:);
h2.EdgeColor = 'n';
xlabel('Input value (a.u.)');
ylabel('Frequency');
% yticks([]);
savefigs(h, filename, plotdir, fontsize, [2, 2]);

h = figure; hold on;
filename = 'Output_baseline';
h3 = histogram(R(2000:end,1),100);
h3.FaceColor = 'r'; %myred(3,:);
h3.EdgeColor = 'n';
h4 = histogram(R(2000:end,2),100);
h4.FaceColor = 'b'; %myblue(3,:);
h4.EdgeColor = 'n';
xlabel('Firing rates (Hz)');
ylabel('Frequency');
% yticks([]);
savefigs(h, filename, plotdir, fontsize, [2, 2]);
%% Value input course
c = .064;
Vinput = [1 + c, 1 - c]*scale;
sgmInput = sgmInput_rprsnt;
filename = sprintf('LDDM_timeCourse_a%1.2f_b%1.2f_c%1.3f_sgm%1.1fsinpt%0.3f',a0,b0,c,sgm,sgmInput);
output = fullfile(Simdir,[filename, '.mat']);
if ~exist(output, 'file')
    rng(4);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    save(output,'R','G','D','Vcourse');
else
    load(output);
end

h = figure; hold on
filename = sprintf('Vcourse_c%1.3f',c);
x = (1:numel(Vcourse(1:5000,1)))*dt;
plot(x,Vcourse(1:5000,1), 'r-');
plot(x,Vcourse(1:5000,2), 'b-');
xlim([0, dur/20]);
xlabel('Time (s)');
ylabel('Input value (a.u.)');
savefigs(h, filename, plotdir, fontsize, [3.2, 2]);

%% General property 2. predicted choice behavior
task = 'RT';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
sgm = .01;
sgmInput = sgmInput_choice;
w = w0*ones(2);
a = a0*eye(2);
b = b0*eye(2);
thresh = 70;
stoprule = 1;
initialvals = [R0,R0; G0,G0; D0,D0];
filename = sprintf('LDDM_%s_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',task,a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    V1 = scale*(1+cp);
    V2 = scale*(1-cp);
    Vinput.V1 = V1; %scale*[1+cp, 1-cp];%*boost(level);
    Vinput.V2 = V2;
    Vprior.V1 = ones(size(V1))*scale;
    Vprior.V2 = ones(size(V2))*scale;
    [rt, choice, ~] = LDDM_Rndinput_GPU(Vprior, Vinput, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    ACC = gather(mean(2-squeeze(choice),2,'omitnan'));
    meanRT = gather(mean(squeeze(rt),2,'omitnan'));
    save(output,'ACC','meanRT');
else
    load(output);
end
h = figure;
filename = 'RT_ACC_baseline';
subplot(1,2,1); hold on;
plot(cp*100,ACC*100,'k.-','LineWidth',lwd,'MarkerSize',mksz);
ylim([50,100]);
xticks([1,10,100]);
xticklabels({'1','10','100'});
set(gca,'XScale','log');
ylabel('% Correct');
xlabel('Value difference (%)'); 
savefigs(h, filename, plotdir, fontsize, [4, 2]);

subplot(1,2,2); hold on;
plot(cp*100,meanRT,'k.-','LineWidth',lwd,'MarkerSize',mksz);
xticks([1,10,100]);
xticklabels({'1','10','100'});
set(gca,'XScale','log');
ylabel('RT (s)');
xlabel('Value difference (%)'); 
savefigs(h, filename, plotdir, fontsize, [4, 2]);

%% SPRT choice
% integrater = [];
% tmp = 0;
% for ti = 1:length(Vcourse)
%    tmp = tmp + (Vcourse(ti,1) - Vcourse(ti,2));
%    integrater(ti)  = tmp;
% end
% 
% h = figure;
% subplot(1,2,1); hold on;
% plot(Vcourse);
% subplot(1,2,2); hold on;
% plot(integrater);
%% E-E only, input & self-excitation
E_Eonly;
%% I-E only, G -> R
I_Eonly;
%% Both E-E & I-E
I_E_E_EBoth;

