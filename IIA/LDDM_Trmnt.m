%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = 'rslts/LDDM_Trmnt';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
%% initalize parameters as in the paper of Wong & Wang, 2006
task = 'WMTrmnt';
a = 118; %118;
b = eye(3)*.6;%3.1;
w = ones(3);
eqlb = 47;
scale = 3*eqlb.^2 + (1-a(1,1)).*eqlb;
tauR = .1; %0.185325; %.1; % second
tauG = 	.1; %0.224459; %.1; % second
tauI = 	.1; %.323132; %.1; % second
Tau = [tauR, tauG, tauI];
% simulation parameters
dt = .001; % second
presentt = dt;
sgm = 5;
stimdur = 2 - presentt; %max(durvec);
durvec = 2 + [.005, .01, .02, .03, .04, .07, .1, .2, .3, .4, .5, 1, 2]; % second
triggert = max(durvec); %presentt;
thresh = 70; % Hz
initialvals = [1,1,1;3,3,3; 0, 0, 0]; % for R, G, and I variables
stoprule = 1;
sims = 10240*2; % number of iterations
simname = sprintf('IIA_LDDM_%s_a%1.2f_b%1.2f_eqlb%1.2f_%iTermnt_scale%4.0f_sgm%2.1f',task,a,b(1,1),eqlb,numel(durvec),scale,sgm);
%% plot psychometric function
h = figure; hold on;
figname = simname;
aspect = [7, 6];
% plot conditional ratio as a function of V3
c = [.128];
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.04:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
% simulation
% filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_a%2.1f_w%1.1f_b%1.2f_scale%2.1f_act%1.2f_sgm%2.1f_sim%i',...
%     task,length(c1),length(unique(c2)),length(unique(c3)),a(1,1),w(1,1),b(1,1),scale,triggert,sgm,sims);
mycol2 = [colormap(winter(floor(length(durvec)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(durvec)/2))))];
lgd3 = [];
for di = 1:numel(durvec)
    dur = durvec(di);
    filename = sprintf('%s_%ic1_%ic2_%ic3_dur%1.3f',simname,length(c1),length(unique(c2)),length(unique(c3)),dur);
    simrslt = fullfile(outdir,[filename, '.mat']);
    if ~exist(simrslt,'file')
        [choice, rt, argmaxR] = LDDM_GPU3ABS(cp, eqlb, w, a, b, sgm, Tau, dur,...
            dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        save(simrslt,'choice','rt','argmaxR');
    else
        load(simrslt);
    end
    % plot
    cndratio = sum(argmaxR == 1, 3)./(sum(argmaxR == 1, 3) + sum(argmaxR == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
    mycol = colormap(jet(length(c3)));
    for i = 1:length(c3)
        plot(cp.cp3(i,1), cndratio(i,1),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
    end
    lgd3(di) = plot(cp.cp3(:,1),cndratio(:,1),'-','color',mycol2(di,:),'LineWidth',lwd/2);
end
lgd = legend(lgd3,cellstr(num2str((durvec-2)')),...
    'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
    'FontAngle','italic','NumColumns',1,'Box','off');
title(lgd,'Cutoff time (s)');
xlabel('V_3');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
savefigs(h, figname, outdir, fontsize, aspect);