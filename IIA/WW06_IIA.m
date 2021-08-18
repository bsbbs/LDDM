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
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
colorpalettergb = [239,71,111;255,209,102;6,214,160;17,138,178;7,59,76]/255;
%% simulation begin
% inital parameters as in the paper of Wong & Wang, 2006
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
h = figure;
filename = 'IIA_WW06';
aspect = [7, 3];
subplot(1,2,1);
hold on;
mycol = colormap(jet(41));
for i = 1:size(x,1)
    plot(x(i,:),cndratio(i,:),'.-','color',mycol(1+(i-1)*10,:),'LineWidth',lwd-1,'MarkerSize',mksz*3);
end
xlabel('V_1 - V_2');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
lgd = legend(string(num2cell(c3)), 'FontSize', fontsize - 6, 'Location','best', 'box','off');
title(lgd, 'V_3');
savefig(h, filename, outdir, fontsize, aspect);
%% plot conditional ratio as a function of V3
c = [.064]';
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
subplot(1,2,2); hold on;
mycol = colormap(jet(length(c3)));
for i = 1:length(c3)
    plot(cp.cp3(i,1), cndratio(i,1),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
end
ylim([.5,.9]);
xlabel('V_3');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
savefig(h, filename, outdir, fontsize, aspect);
