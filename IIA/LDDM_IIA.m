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
task = 'FDon';%'RT';%'FD';
a = eye(3)*15;
b = eye(3)*1.02;
w = ones(3);
scale = 256;
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];
% simulation parameters
dt = .001; % second
presentt = dt;
if strcmp(task,'FD')
    sgm = 5;
    dur = 7; % second
    triggert = 2.5; %presentt;
    stimdur = 2.5;
elseif strcmp(task,'FDon')
    sgm = 10;
    dur = 7; % second
    triggert = 2.5; %presentt;
    stimdur = dur;
elseif strcmp(task,'RT')
    sgm = 20;
    dur = 5; % second
    triggert = presentt;
    stimdur = dur;
end
thresh = 70; % Hz
initialvals = [4, 4, 4;12, 12, 12; 0, 0, 0]; % for R, G, and I variables
stoprule = 1;
sims = 10240*2; % number of iterations
%% plot example dynamics
h = figure; hold on;
c1 = 1.064;
c2 = 1;
for c3 = [.2 1]
    Vinput = [c1, c2, c3]*scale;
    [choice, rt, R, G, I] = LDDM(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    plot(R(:,1:2));
end
ylim([0,80]);
%% plot psychometric function
h = figure;
figname = sprintf('IIA_LDDM_%s',task);
aspect = [7, 6];
% setting input matrix
c = [.032, .064, .128, .256, .512]';
c1 = [1 - flip(c); 1; 1 + c];
c2 = ones(size(c1));
c3 = 0:.5:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
clear Vinput;
Vinput.V1 = cp.cp1*scale;
Vinput.V2 = cp.cp2*scale;
Vinput.V3 = cp.cp3*scale;
% simulation
filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_sgm%2.1f_sim%i',task,length(c1),length(unique(c2)),length(unique(c3)),sgm,sims);
simrslt = fullfile(outdir,[filename, '.mat']);
if ~exist(simrslt,'file')
    [choice, rt] = LDDM_GPU3(Vinput, w, a, b, sgm, Tau, dur,...
        dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    save(simrslt,'choice','rt');
else
    load(simrslt);
end
% plot psychometric function
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
if strcmp(task,'FD')
    c = [.512]';
elseif strcmp(task,'FDon')
    c = [.128]';
elseif strcmp(task,'RT')
    c = [.256]';
end
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.04:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
Vinput.V1 = cp.cp1*scale;
Vinput.V2 = cp.cp2*scale;
Vinput.V3 = cp.cp3*scale;
% simulation
filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_sgm%2.1f_sim%i',task,length(c1),length(unique(c2)),length(unique(c3)),sgm,sims);
simrslt = fullfile(outdir,[filename, '.mat']);
if ~exist(simrslt,'file')
    [choice, rt] = LDDM_GPU3(Vinput, w, a, b, sgm, Tau, dur,...
        dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
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
if strcmp(task,'FD')
    ylim([.8,.9]);  
elseif strcmp(task,'FDon')
    ylim([.75,1]);  
elseif strcmp(task,'RT')
    ylim([.75,.95]);
end
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