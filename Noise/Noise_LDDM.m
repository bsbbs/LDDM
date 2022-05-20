%% To test different types and magnitudes of noise to value representation and irrelavant dependence of choice
% include packages
addpath('/Users/bs3667/Documents/LDDM/CoreFunctions');
addpath('/Users/bs3667/Documents/LDDM/utils/');
% set path
% Homedir = 'C:\Users\Bo';
Homedir = '~';
% Glgdir = 'G:\My Drive';
Glgdir = '/Volumes/GoogleDrive/My Drive';
out_dir = fullfile(Glgdir, 'LDDM/Noise');
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
% parameters for visulization
fontsize = 14;
mksz = 25;
lwd = 1.5;
cp = [.032, .064,.128, .256];
mygray = flip(gray(numel(cp) + 2));
myred = mygray; myred(:,1) = 1;
myblue = mygray; myblue(:,3) = 1;
%% parameters for the model
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauD = .1; % second
Tau = [tauR, tauG, tauD];
w0 = 1;
a0 = 15;
b0 = 1;
c_rprsnt = .064;
c_choice = .032;
predur = 0;
presentt = dt;
thresh = 70; % Hz
stoprule = 1;
eqlb = 32;
scale = (2*w0 - b0)*eqlb^2 + (1-a0)*eqlb; %2*w0*eqlb.^2 + (1-a0).*eqlb;
R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0)))/2/(2*w0 - b0);
D0 = b0*R0;
G0 = (2*w0-b0)*R0;
Vprior = [1, 1]*scale;
%% Step 1. check the representation under different types/magnitudes of noise
predur = 0;
dur = 100;
c = c_rprsnt;
a = a0*eye(2);
b = b0*eye(2);
w = ones(2);
triggert = Inf;
stimdur = Inf;
presentt = dt;
initialvals = [1,1;2,2;0,0]*eqlb/3;
Vprior = [1, 1]*scale;
Vinput = [1+c, 1-c]*scale;
sgm_internal_vec = [.01, 5];
sgm_external_vec = [0.01, .33];
filename = sprintf('LDDM_timeCourse_a%1.2f_b%1.2f_sgm_intrnl%1.2f_%1.2f_sgm_extrnl%1.2f_%1.2f',...
    a(1),b(1),sgm_internal_vec(1),sgm_internal_vec(2),sgm_external_vec(1),sgm_external_vec(2));
output = fullfile(Simdir,[filename, '.mat']);
if ~exist(output, 'file')
    R_conds = [];
    Vcourse_conds = [];
    for i = 1:numel(sgm_internal_vec)
        sgm_internal = sgm_internal_vec(i);
        for j = 1:numel(sgm_external_vec)
            sgm_external = sgm_external_vec(j);
            rng(5);
            [choice, rt, R, G, I, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
                sgm_internal, sgm_external*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R_conds{i,j} = R;
            Vcourse_conds{i,j} = Vcourse;
        end
    end
    save(output, 'R_conds','Vcourse_conds', 'sgm_internal_vec', 'sgm_external_vec');
else
    load(output);
end
%% Input and output dynamics
h = figure;
i = 1;
for j = 1:numel(sgm_external_vec)
    subplot(3,2,j); hold on;
    sgm_external = sgm_external_vec(j);
    x = (1:numel(R_conds{i,j}(2:2000,1)))*dt;
    plot(x, Vcourse_conds{i,j}(2:2000,1), 'r-', 'LineWidth', lwd);
    plot(x, Vcourse_conds{i,j}(2:2000,2), 'b-', 'LineWidth', lwd);
    xlim([-50*dt, 2]);
    title(['\sigma_{Exo} = ' num2str(sgm_external)]);
    xlabel('');
    if j == 1
        ylabel('Input values (a.u.)');
    else
        ylabel('');
    end
    ylim([-1000,2000]);
    savefigs(h, filename, plotdir, fontsize, [8,8]);
end
gi = 2;
for i = 1:numel(sgm_internal_vec)
    sgm_internal = sgm_internal_vec(i);
    for j = 1:numel(sgm_external_vec)
        sgm_external = sgm_external_vec(j);
        gi = gi + 1;
        subplot(3,2,gi); hold on;
        x = (1:numel(R_conds{i,j}(1:2000,1)))*dt;
        plot(x, R_conds{i,j}(1:2000,1), 'r-', 'LineWidth', lwd);
        plot(x, R_conds{i,j}(1:2000,2), 'b-', 'LineWidth', lwd);
        xlim([-50*dt, 2]);
        % plot(x, G(:,1), 'r--', 'LineWidth', lwd);
        % plot(x, G(:,2), 'b--', 'LineWidth', lwd);
        % plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
        % plot(x, D(:,2), 'b-.', 'LineWidth', lwd);
        yticks([10, 15, 20, 25]);
        xlabel('');
        ylabel('');
        title(['\sigma_{Exo} = ' num2str(sgm_external) ', \sigma_{Endo} = ' num2str(sgm_internal)]);
        if gi == 3 || gi == 5
            ylabel('Firing rates (Hz)');
        end
        if gi >= 5
            xlabel('Time (msecs)');
        end
        if gi == 3
            legend({'R_1','R_2'},'Box','off','Location','southeast',...
                'FontName','Arial', 'FontAngle','italic',...
                'FontSize',fontsize-4);
        end
        savefigs(h, filename, plotdir, fontsize, [8,8]);
    end
end


%% Input and output distribution, smoothed curve
filename = sprintf('LDDM_Dstrbtn_a%1.2f_b%1.2f_sgm_intrnl%1.2f_%1.2f_sgm_extrnl%1.2f_%1.2f',...
    a(1),b(1),sgm_internal_vec(1),sgm_internal_vec(2),sgm_external_vec(1),sgm_external_vec(2));
h = figure;
i = 1;
for j = 1:numel(sgm_external_vec)
    subplot(3,2,j); hold on;
    sgm_external = sgm_external_vec(j);
    x = 0:1000;
    pd1 = fitdist(Vcourse_conds{i,j}(2000:end,1),'kernel','Kernel','normal');
    y1 = pdf(pd1,x);
    pd2 = fitdist(Vcourse_conds{i,j}(2000:end,2),'kernel','Kernel','normal');
    y2 = pdf(pd2,x);
    plot(x,y1,'r--');
    plot(x,y2,'b--');
    xlabel('Input values (a.u.)');
    if j == 1
        ylabel('Density');
    else
        ylabel('');
    end
    savefigs(h, filename, plotdir, fontsize, [8, 8]);
end

gi = 2;
for i = 1:numel(sgm_internal_vec)
    sgm_internal = sgm_internal_vec(i);
    for j = 1:numel(sgm_external_vec)
        sgm_external = sgm_external_vec(j);
        gi = gi + 1;
        subplot(3,2,gi); hold on;
        x = 10:.1:38;
        pd3 = fitdist(R_conds{i,j}(2000:end,1),'kernel','Kernel','normal');
        y3 = pdf(pd3,x);
        pd4 = fitdist(R_conds{i,j}(2000:end,2),'kernel','Kernel','normal');
        y4 = pdf(pd4,x);
        ln3 = plot(x,y3,'r-');
        ar3 = area(x,y3,'FaceColor','r','FaceAlpha',.5,'EdgeAlpha',.3);
        ln4 = plot(x,y4,'b-');
        ar4 = area(x,y4,'FaceColor','b','FaceAlpha',.5,'EdgeAlpha',.3);
        xlim([10, 36]);
        xlabel('');
        ylabel('');
        title(['\sigma_{Exo} = ' num2str(sgm_external) ', \sigma_{Endo} = ' num2str(sgm_internal)]);
        if gi == 3 || gi == 5
            ylabel('Density');
        end
        if gi >= 5
            xlabel('Firing rates (Hz)');
        end
        if gi == 3
            legend([ar3, ar4],{'R_1','R_2'},'Box','off','Location','southeast',...
                'FontName','Arial', 'FontAngle','italic',...
                'FontSize',fontsize-4);
        end
        savefigs(h, filename, plotdir, fontsize, [8, 8]);
    end
end
%% Step 2. predicted choice behavior, relative choice ratio (V1 vs. V2) and RT as function of V3
eqlb = 32;
scale = (2*w0 - b0)*eqlb^2 + (1-a0)*eqlb; % (3*w0 - b0)*eqlb^2 + (1-a0)*eqlb;
eqlb3 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(3*w0 - b0)))/2/(3*w0 - b0);
R0 = eqlb3;
D0 = b0*R0;
G0 = (2*w0-b0)*R0;
task = 'RT';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
w = w0*ones(3);
a = a0*eye(3);
b = b0*eye(3);
thresh = 70;
stoprule = 1;
initialvals = [R0,R0,R0; G0,G0,G0; D0,D0,D0];
sgm_internal_vec = [.01, 5];
sgm_external_vec = [0.01, .33];
c = c_choice;
V3 = [0:.05:1]';
Vprior = [ones(size(V3)), ones(size(V3)), ones(size(V3))]*scale;
Vinput = [ones(size(V3))*(1+c), ones(size(V3))*(1-c), V3]*scale;
for i = 1:numel(sgm_internal_vec)
    sgm_internal = sgm_internal_vec(i);
    for j = 1:numel(sgm_external_vec)
        sgm_external = sgm_external_vec(j);
        filename = sprintf('LDDM_%s_a%1.2f_b%1.2f_sgm_intrnl%1.2f_sgm_extrnl%1.2f',...
            task, a(1), b(1), sgm_internal, sgm_external);
        output = fullfile(Simdir,[filename, '.mat']);
        ACC = [];
        meanRT = [];
        if ~exist(output,'file')
            [rt, choice, argmaxR] = LDDM3_Rndinput_GPU(Vprior, Vinput, w, a, b,...
                sgm_internal, sgm_external*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            chosenratio = ();
            ACC = (mean(2-squeeze(choice),2,'omitnan'));
            meanRT = (mean(squeeze(rt),2,'omitnan'));
            save(output,'ACC','meanRT');
        else
            load(output);
        end
    end
end
%
mygray = [.9:-.9/numel(potentiation):0];
h = figure;
filename = 'RT_ACC_I_E';
subplot(2,1,1); hold on;
for level = 1:numel(potentiation)
    plot(cp,ACC(level,:),'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz/2,'LineWidth',lwd);
end
% set(gca,'XScale','log');
ylabel('% Correct');
lgd = legend(cellstr(string(potentiation)),...
    'Location','SouthEast', 'FontSize', fontsize-4, 'Box','off');
title(lgd,'STDP','FontSize',fontsize-4);
savefigs(h, filename, plotdir,fontsize, [2, 3]) ;
subplot(2,1,2); hold on;
for level = 1:numel(potentiation)
    plot(cp,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz/2,'LineWidth',lwd);
end
% set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, [2, 3]);

%% Property 3. ACC as a function of Input noise, two STDP levels
task = 'RT_I-E';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
sgm_internal = .01;
w = w0*ones(2);
a = a0*eye(2);
b = b0*eye(2);
sims = 102400;
potentiation = [1, 2];
sgmInputvec = linspace(0,1.4,100);
filename = sprintf('LDDM_%s_STDP%.1f_%.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3fto%.3f_sims%i',task,min(potentiation),max(potentiation),a0,b0,sgm_internal,min(sgmInputvec),max(sgmInputvec),sims);
output = fullfile(Simdir,[filename, '.mat']);
STDP_v = 1;
STDP_a = 1;
c = 3.2/100;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    for pi = 1:numel(potentiation)
        STDP_G = potentiation(pi);
        fprintf('potentiation %.3f\n',potentiation(pi));
        [V1, ~] = meshgrid(scale*(1+c), sgmInputvec);
        [V2, sgm_external] = meshgrid(scale*(1-c), sgmInputvec);
        Vinput.V1 = V1;
        Vinput.V2 = V2;
        Vprior.V1 = ones(size(V1))*scale;
        Vprior.V2 = ones(size(V2))*scale;
        [rt, choice, ~] = LDDM_Rndinput_STDP_GPU(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
            sgm_internal, sgm_external*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC(pi,:) = gather(mean(2-squeeze(choice),2,'omitnan'));
        meanRT(pi,:) = gather(mean(squeeze(rt),2,'omitnan'));
    end
    save(output,'ACC','meanRT');
else
    load(output);
end
%
h = figure; hold on;
plot(sgmInputvec, ACC(1,:),'MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(1)*[1,1,1]);
plot(sgmInputvec, ACC(2,:),'MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(3)*[1,1,1]);
legend({'Baseline','STDP'},'Box','off','Location','northeast');
xlabel('Input noise (a.u.)');
ylabel('Accuracy');
savefigs(h, ['ACCoversgmInput_', filename], plotdir,fontsize, [2,1.5]);

%% largest amplitude and iSTDP
task = 'VR_I-E';
potentiation = [1:.5:2.5];
c = 3.2/100;
sgm_internal = 0;
sgm_external = 0;
triggert = Inf;
stimdur = Inf;
dur = 4.5;
STDP_v = 1;
STDP_a = 1;
initialvals = [1,1;2,2;0,0]*eqlb/3;
Rstore = [];
h = figure;
filename = sprintf('R dynamic over STDP %s',task);
for level = 1:length(potentiation)
    Vinput = [1+c, 1-c]*scale;
    a = a0*eye(2);
    b = 0*eye(2);
    w = w0*ones(2);
    STDP_G = potentiation(level);
    [choice, rt, R, G, I, Vcourse] = LDDM_RndInput_STDP(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
    sgm_internal, sgm_external*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    Rstore(level) = max(R(:,1));
    subplot(1,4,level);
    plot(R);
    ylim([10,40]);
    xlabel('Time (a.u.)');
    ylabel('Firing rates (Hz)');
end
savefigs(h, [filename], plotdir,fontsize-5, [6,2]);

h = figure; hold on;
filename = sprintf('MaxFR_%s',task);
plot(potentiation,Rstore,'k-','MarkerSize', mksz,'LineWidth',lwd);
xlabel('STDP');
ylabel('Max firing rates (Hz)');
savefigs(h, [filename], plotdir,fontsize, [2,1.5]);