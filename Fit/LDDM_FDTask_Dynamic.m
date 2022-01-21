% LDDM Dynamics in Fixed-duration task 
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Fit');
cd('G:\My Drive\LDDM\Fit');
out_dir = './Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
dataDynmc = load('./Data/Data.mat');
dataBhvr = LoadRoitmanData('../RoitmanDataCode');
randseed = 69094639;
rng(randseed);
% a, b, noise, scale, tauRGI, nLL
params = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132	16539.138186];
name = sprintf('LDDM_FD_a%2.2f_b%1.2f_sgm%2.1f_scale%4.1f_tau%1.2f_%1.2f_%1.2f_nLL%4.0f',params);


% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax';
sac_ax = dataDynmc.sac_ax';
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
OP = dataBhvr.OP;
% parameters to fit
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
tauR = params(5);
tauG = params(6);
tauI = params(7);
ndt = .09 + .03; % sec, 90ms after stimuli onset, resort to the saccade side,
% the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
presentt = 0; % changed for this version to move the fitting begin after the time point of recovery
scale = params(4);

% other fixed parameters
% sims = 1024;
deduction = 1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
predur = 0;
triggert = 1+1.5;%-2.5;
dur = 5;
dt =.001;
thresh = 70; %70.8399; % mean(max(m_mr1cD))+1; 
stimdur = 2.5;
stoprule = 1;
w = [1 1; 1 1];
Rstar = 32; % ~ 32 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar; 0,0];
Vprior = ones(6,2)*((1-a(1,1))*Rstar + 2*Rstar^2);

Tau = [tauR tauG tauI];
% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;
if ~exist(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),'file')
    tic;
    [rtmat, choicemat, ~, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDM_Dynmc_Trim_GPU(Vprior, Vinput, w, a, b,...
        sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
    rtmat = squeeze(rtmat)'+ndt;
    choicemat = squeeze(choicemat)';
    save(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),...
        'rtmat','choicemat','sm_mr1c','sm_mr2c','sm_mr1cD','sm_mr2cD','params');
    toc
else
    load(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)));
end

%% plot time course
load('./Data/Data.mat');
m_mr1c = m_mr1c';
m_mr2c = m_mr2c';
m_mr1cD = m_mr1cD';
m_mr2cD = m_mr2cD';
dot_ax = dot_ax';
sac_ax = sac_ax';
h = figure;
aspect = [5, 2.5];
fontsize = 10;
lwd = 1;
filename = sprintf('FittedTimeCourse_%s',name);
subplot(1,2,1);hold on;
clear flip;
colvec = flip({[218,166,109]/256,[155 110 139]/256,'#32716d','#af554d','#708d57','#3b5d64'});
for ci = 1:6
    lg(ci) = plot(dot_ax/1000, sm_mr1c(:,ci),'Color',colvec{ci},'LineWidth',lwd);
    plot(dot_ax/1000, sm_mr2c(:,ci),'--','Color',colvec{ci},'LineWidth',lwd);
end
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
% ylim([20,60]);
ylim([20,70.5]);
ylabel('Firing rate (sp/s)');
xlabel('Time (secs)');
xlim([-.05, .8]);
xticks([0:.2:.8]);
% set(gca,'FontSize',16);
savefigs(h,filename,plot_dir,fontsize,aspect);
subplot(1,2,2);hold on;
plot([0,0],[20,71],'-k');
for ci = 1:6
    lg(ci) = plot(sac_ax/1000, sm_mr1cD(:,ci),'Color',colvec{ci},'LineWidth',lwd);
    plot(sac_ax/1000, sm_mr2cD(:,ci),'--','Color',colvec{ci},'LineWidth',lwd);
end
xlim([-.8, .05]);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
yticks([]);
set(gca,'ycolor',[1 1 1]);
ylim([20,70.5]);
legend(lg,{'0','3.2','6.4','12.8','25.6','51.2'},'Location','best','FontSize',fontsize-2);
savefigs(h,filename,plot_dir,fontsize,aspect);
saveas(h,fullfile(plot_dir,[filename, '.fig']),'fig');
