%% ACC and RT as a function of GABAergic agonist level (GABAergic activation), under different levels of input
%% preparation
addpath('../');
addpath(genpath('../../../CoreFunctions'));
outdir = fullfile('./graphics','GABAActivation_ACC_RT_Values');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end

%% simulation for LDDM model
dt = .001;
w = ones(2,2);
a = eye(2,2)*10;
b = eye(2,2)*1.1;
sgm = 2;
presentt = dt;
triggert = dt;
dur = 5;
stimdur = dur;
thresh = 120;
Tau = [.1, .1, .1]*2;
initialvals = [35,35; 70,70; 0,0];
stoprule = 1;
sims = 10000;
cohr = [0 32 64 128 256 512]'/1000;
VmatDiag = 256*[1+cohr, 1-cohr];
mygray =  gray(2+length(cohr));
% Input
% h = figure; hold on;
% for vi = 1:length(cohr)
%     Vinput = VmatDiag(vi,:);
%     plot(Vinput(1),Vinput(2),'k.','MarkerSize',mksz/2,'Color',mygray(vi+1,:));
% end
% xlabel('V_1 (a.u.)');ylabel('V_2 (a.u.)');
% xlim([0,512]);ylim([0,512]);
% xticks([0,512]);yticks([0,512]);
% set(gca,'FontSize',fontsize-7);
% set(gca,'TickDir','out');
% set(gca,'LineWidth',1);
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 3.1 2.8]/3;
% saveas(h,fullfile(outdir,sprintf('ValueInputMatrix.eps')),'epsc2');
GABALevel = 10.^[linspace(0,.6021,21)];%1:10;
filename = sprintf('LcD_GABA%i_%.1f_%1.1fa%1.0fb%1.0fsgm%1.1f_values%i_sims%i',...
    length(GABALevel),min(GABALevel),max(GABALevel),a(1,1),b(1,1),sgm,numel(cohr),sims);
output = fullfile('./SimRslts',[filename, '.mat']);
if ~exist(output,'file')
    RT = [];
    Choice = [];
    for vi = 1:length(cohr)
        rng(2021);
        V1 = 256*(1 + cohr(vi));
        V2 = 256*(1 - cohr(vi));
        Vinput = [V1, V2];
        parfor gi = 1:length(GABALevel)
            Gaba = GABALevel(gi);
            [rt, choice, ~] = LcDsInhbt_Gaba_GPU(Vinput, Gaba, w, a, b,...
                sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            RT(:,gi,vi) = gather(squeeze(rt));
            Choice(:,gi,vi) = gather(squeeze(choice));
        end
    end
    save(output,'RT','Choice','GABALevel','a','b','cohr');
else
    load(output);
end
% plot
meanRT = squeeze(mean(RT,1,'omitnan'));
meanACC = squeeze(mean(2-Choice,1,'omitnan'));
h = figure;
figname = 'LcD_ACC_RT_GABA_Value';
subplot(2,1,1); hold on;
for vi = 1:length(cohr)
    plot(GABALevel,meanRT(:,vi),'.-','MarkerSize',mksz-4,'LineWidth',lwd,'Color',mygray(vi,:));
end
lgd = legend({'0','3.2','6.4','12.8','25.6','51.2'},'Box','off','FontSize',fontsize-9,'Location','Best');
title(lgd,'Coherence (%)');
%xlabel('GABAergic activation');
ylabel('RT (secs)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
subplot(2,1,2); hold on;
for vi = 1:length(cohr)
    plot(GABALevel,meanACC(:,vi)*100,'.-','MarkerSize',mksz-4,'LineWidth',lwd,'Color',mygray(vi,:));
end
xlabel('GABAergic activation');
ylabel('Accuracy (%)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);

%% RT and ACC as function of Input, under two levels of GABAergic activation
GABALevel = [1, 1.8];
c = flip([10.^linspace(-2,0,101)]');
Vinput = 256*[1+c,1-c];
sims = 10000;
filename = sprintf('LcD_GABA%i_%.1f_%1.1fa%1.0fb%1.0fsgm%1.1f_values%i_sims%i',...
    length(GABALevel),min(GABALevel),max(GABALevel),a(1,1),b(1,1),sgm,numel(c),sims);
output = fullfile('./SimRslts',[filename, '.mat']);
if ~exist(output,'file')
    RT = [];
    Choice = [];
    for gi = 1:length(GABALevel)
        GABA = GABALevel(gi);
        [rt, choice, ~] = LcDsInhbt_Gaba_GPU(Vinput, GABA, w, a, b,...
            sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        RT(:,:,gi) = gather(squeeze(rt)');
        Choice(:,:,gi) = gather(squeeze(choice)');
    end
    save(output,'RT','Choice','GABALevel','a','b','Vinput');
else
    load(output);
end

meanRT = squeeze(mean(RT,1,'omitnan'));
meanACC = squeeze(mean(2-Choice,1,'omitnan'));
h = figure;
figname = 'LcD_ACC_RT_Value_GABA';
subplot(2,1,1); hold on;
for gi = 1:length(GABALevel)
    plot(c,meanRT(:,gi),'-','MarkerSize',mksz-4,'LineWidth',lwd);
end

%xlabel('GABAergic activation');
ylabel('RT (secs)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
subplot(2,1,2); hold on;
for gi = 1:length(GABALevel)
    plot(c,meanACC(:,gi)*100,'-','MarkerSize',mksz-4,'LineWidth',lwd);
end
lgd = legend({'Placebo','Agonist'},'Box','off','FontSize',fontsize-5,'Location','Best');
xlabel('Choice difficulty (V_1 - V_2)');
ylabel('Accuracy (%)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);

%% simulation for Wong06 model
paramspecify_WongWang;
dt = .001;
sgm = .02;
sgmInput = 2/3;
presentt = dt;
dur = 5;
stimdur = dur;
thresh = 15;
initialvals = initialvals_dflt;
stoprule = 1;
sims = 10000;
cohr = [0 32 64 128 256 512]'/1000; %
VmatDiag = 256*[1+cohr, 1-cohr];
mygray =  gray(2+length(cohr));
GABALevel = 10.^[linspace(0,.6021,21)];%1:10;
filename = sprintf('XJ_GABA%i_%.1f_%1.1fsgm%0.2f_sgmInput%0.2f_values%i_sims%i',...
    length(GABALevel),min(GABALevel),max(GABALevel),sgm,sgmInput,numel(cohr),sims);
output = fullfile('./SimRslts',[filename, '.mat']);
if ~exist(output,'file')
    RT = [];
    Choice = [];
    for vi = 1:length(cohr)
        rng(2021);
        V1 = 256*(1 + cohr(vi));
        V2 = 256*(1 - cohr(vi));
        Vinput = [V1, V2];
        parfor gi = 1:length(GABALevel)
            Gaba = GABALevel(gi);
            [rt, choice, argmaxR, dR] = wong06_RndInput_Gaba_GPU(Vinput,Gaba,miu0,sgm,sgmInput,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims);
            RT(:,gi,vi) = gather(squeeze(rt));
            Choice(:,gi,vi) = gather(squeeze(choice));
        end
    end
    save(output,'RT','Choice','GABALevel','a','b','cohr');
else
    load(output);
end
% plot
meanRT = squeeze(mean(RT,1,'omitnan'));
meanACC = squeeze(mean(2-Choice,1,'omitnan'));
h = figure;
figname = 'XJ_ACC_RT_GABA_Value';
subplot(2,1,1); hold on;
for vi = 1:length(cohr)
    plot(GABALevel,meanRT(:,vi),'.-','MarkerSize',mksz-4,'LineWidth',lwd,'Color',mygray(vi,:));
end
% lgd = legend({'0','3.2','6.4','12.8','25.6','51.2'},'Box','off','FontSize',fontsize-9,'Location','Best');
% title(lgd,'Coherence (%)');
% xlabel('GABAergic activation');
ylabel('RT (secs)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
subplot(2,1,2); hold on;
for vi = 1:length(cohr)
    plot(GABALevel,meanACC(:,vi)*100,'.-','MarkerSize',mksz-4,'LineWidth',lwd,'Color',mygray(vi,:));
end
xlabel('GABAergic activation');
ylabel('Accuracy (%)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
%% RT and ACC as function of Input, under two levels of GABAergic activation
GABALevel = [1, 1.8];
c = flip([10.^linspace(-2,0,101)]');
Vinput = 256*[1+c,1-c];
sims = 10000;
filename = sprintf('XJ_GABA%i_%.1f_%1.1fsgm%1.1f_sgmInput%1.1f_values%i_sims%i',...
    length(GABALevel),min(GABALevel),max(GABALevel),sgm,sgmInput,numel(cohr),sims);
output = fullfile('./SimRslts',[filename, '.mat']);
if ~exist(output,'file')
    RT = [];
    Choice = [];
    for gi = 1:length(GABALevel)
        GABA = GABALevel(gi);
        [rt, choice, argmaxR, dR] = wong06_RndInput_Gaba_GPU(Vinput,GABA,miu0,sgm,sgmInput,I0,JN,...
            gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims);
        RT(:,:,gi) = gather(squeeze(rt)');
        Choice(:,:,gi) = gather(squeeze(choice)');
    end
    save(output,'RT','Choice','GABALevel','a','b','Vinput');
else
    load(output);
end

meanRT = squeeze(mean(RT,1,'omitnan'));
meanACC = squeeze(mean(2-Choice,1,'omitnan'));
h = figure;
figname = 'XJ_ACC_RT_Value_GABA';
subplot(2,1,1); hold on;
for gi = 1:length(GABALevel)
    plot(c,meanRT(:,gi),'-','MarkerSize',mksz-4,'LineWidth',lwd);
end

%xlabel('GABAergic activation');
ylabel('RT (secs)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
subplot(2,1,2); hold on;
for gi = 1:length(GABALevel)
    plot(c,meanACC(:,gi)*100,'-','MarkerSize',mksz-4,'LineWidth',lwd);
end
lgd = legend({'Placebo','Agonist'},'Box','off','FontSize',fontsize-5,'Location','Best');
xlabel('Choice difficulty (V_1 - V_2)');
ylabel('Accuracy (%)');
set(gca, 'XScale','log');
savefig(h, figname, outdir, fontsize, aspect2);
