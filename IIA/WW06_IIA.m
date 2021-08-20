%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = 'rslts';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
%% initalize parameters as in the paper of Wong & Wang, 2006
miu0 = 30; % Hz
sgm = .02; % nA
I0 = .3255; %nA
JN = [.2609, -.0497, -.0497;
    -.0497, .2609, -.0497;
    -.0497, -.0497, .2609]; %nA
gamma = .641;
tauS = .1; % second
tauAMPA = .002; % second
% simulation parameters
dur = 5; % second
dt = .001; % second
presentt = dt;
stimdur = dur;
thresh = 15; % Hz
initialvals = [2, 2, 2; .1, .1, .1]; % for H and S variables
stoprule = 1;
sims = 10240*2; % number of iterations
%% plot psychometric function
h = figure;
figname = 'IIA_WW06';
aspect = [7, 6];
% setting input matrix
c = [.032, .064, .128, .256, .512]';
c1 = [1 - flip(c); 1; 1 + c];
c2 = ones(size(c1));
c3 = 0:.5:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
% simulation
simrslt = fullfile(outdir,sprintf('WW06_%ic1_%ic2_%ic3_sim%i.mat',length(c1),length(unique(c2)),length(unique(c3)),sims));
if ~exist(simrslt,'file')
    [choice, rt] = wong06_GPU3(cp,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims);
    save(simrslt,'choice','rt');
else
    load(simrslt);
end
% plot
cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
x = cp.cp1 - cp.cp2;
subplot(2,2,1);
hold on;
mycol = colormap(jet(41));
for i = 1:size(x,1)
    plot(x(i,:),cndratio(i,:),'.-','color',mycol(1+(i-1)*10,:),'LineWidth',lwd-1,'MarkerSize',mksz*3);
end
xlabel('V_1 - V_2');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
lgd = legend(string(num2cell(c3)), 'FontSize', fontsize - 6, 'Location','best', 'box','off');
title(lgd, 'V_3');
savefig(h, figname, outdir, fontsize, aspect);
% plot RT as a functon of V1 - V2
subplot(2,2,3);
hold on;
mycol = colormap(jet(41));
rtmean = mean(rt,3,'omitnan');
for i = 1:size(x,1)
    plot(x(i,:),rtmean(i,:),'.-','color',mycol(1+(i-1)*10,:),'LineWidth',lwd-1,'MarkerSize',mksz*3);
end
xlabel('V_1 - V_2');
ylabel({'RT (s)'});
savefig(h, figname, outdir, fontsize, aspect);
% plot simulation settings
subplot(2,2,4); hold on;
y = 3*ones(size(c2));
plot(c2, y,'.k','MarkerSize',mksz*3);
y = 3.1*ones(size(c1));
plot(c1, y, '.k','MarkerSize',mksz*3);
% plot conditional ratio as a function of V3
c = [.128]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.05:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
% simulation
simrslt = fullfile(outdir,sprintf('WW06_%ic1_%ic2_%ic3_sim%i.mat',length(c1),length(unique(c2)),length(unique(c3)),sims));
if ~exist(simrslt,'file')
    [choice, rt] = wong06_GPU3(cp,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims);
    save(simrslt,'choice','rt');
else
    load(simrslt);
end
% plot
cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
subplot(2,2,2); hold on;
mycol = colormap(jet(length(c3)));
for i = 1:length(c3)
    plot(cp.cp3(i,1), cndratio(i,1),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
end
ylim([.75,1]);
xlabel('V_3');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
savefig(h, figname, outdir, fontsize, aspect);
% plot simulation settings
subplot(2,2,4); hold on;
y = 2.9*ones(size(c3));
for i = 1:length(c3)
    plot(c3(i),y(i),'.','color',mycol(i,:),'MarkerSize',mksz*3);
end
ylim([2.8,3.5]);
yticks([2.9,3.0,3.1]);
yticklabels({'V_3','V_2','V_1'});
set(get(gca,'YAxis'),'Visible','on');
xlabel('Option values');
savefig(h, figname, outdir, fontsize, aspect);