%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(fullfile(Homedir,'Documents','LDDM','Froemke/SST_Inoise_v'));
cd('G:\My Drive\LDDM\Froemke\SST_Inoise_v');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Froemke');
plotdir = fullfile('./Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end

%% parameters
fontsize = 14;
mksz = 5;
lwd = 1.5;
aspect1 = [3.9,2.2]; % 16:9, for wide temporal dynamic
aspect2 = [3 3]; % for temporal dynamic
aspect3 = [2.8 2.54];
% simulation parameters
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];

w0 = ones(2);
a0 = 40;
b0 = 1.8;
task = 'VR';
%VR, stimdur = Inf, triggert = Inf;
%RT, stimdur = Inf, triggert = presentt;
%FDon, stimdur = Inf, triggert = presentt + .5;
%FDoff, stimdur = .5, triggert = presentt + .5
%WM, stimdur = .5, triggertt = presentt+.7
predur = 0;
presentt = dt;%*200;
stimdur = Inf;
triggert = Inf;%presentt + dt;
dur = 45; % second

thresh = 70; % Hz
stoprule = 1;

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0).*eqlb;

c = [.064]'; % Coherence
c1 = [1 + c];
c2 = [1 - c];
Vprior = [1,1]*scale;

boost = [1, 2];
% Dynamics
sgm = .01;
% sgmInput = 150; for OU process
sgmInput = 1/3;
h = figure;
filename = sprintf('timeCourse_%s_iSTDP%1.1fand%1.1f_sgm%2.2f_sgmInput%.2f',task, boost, sgm, sgmInput);
mycl = jet(length(boost)*10);
for level = 1:length(boost)
    Vinput = [c1, c2]*scale;%*boost(level);
    a = a0;
    b = b0;
    w = w0;%*boost(level);
    iSTDP = boost(level);
    %initialvals = (a-1)/iSTDP./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
    initialvals = (a0-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0];
    [choice, rt, R, G, I, Vcourse] = LDDM_STDP_v_RndInput(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);

    % [choice, rt, R, G, I] = LDDM_STDP(Vprior, Vinput, iSTDP, w, a, b,...
    % sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);

    subplot(3,1,1); hold on;
    plot(R(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    plot(R(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    subplot(3,1,2); hold on;
    plot(G(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    plot(G(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    subplot(3,1,3); hold on;
    plot(I(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    plot(I(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
end
celltype = {'R','G','I'};
for i = 1:3
    subplot(3,1,i); hold on;
    ylim([0,72]);
    xticks([presentt/dt, triggert/dt]);
    xticklabels({'stimuli','action'});
    ylabel('Firing Rates (Hz)');
    title(celltype{i});
end
xlabel('Time (a.u.)');
lgd = legend(' ', ' ', 'V1', 'V2',...
    'Location','Northwest','NumColumns',2, 'FontSize', fontsize-8, 'Box','off');
title(lgd, "Baseline                  iSTDP      .");
savefigs(h, filename, plotdir,fontsize-5, [3,6]);
%
if dur >= 45
    h = figure;
    filename = sprintf('Vcourse_cp%.2f_sgmInput%2.2f',c,sgmInput);
    subplot(2,1,1);
    plot(Vcourse);
    xlabel('Time (a.u.)');
    ylabel('V (a.u.)');
    legend({'V_1','V_2'},'Location','northeastoutside','Box','off','FontSize',fontsize-8);
    subplot(2,1,2); hold on;
    histogram(Vcourse(2:end,1),40);
    histogram(Vcourse(2:end,2),40);
    xlabel('Input value');
    savefigs(h, filename, plotdir,fontsize-5, [3,3]);
    %%
    h = figure;
    filename = sprintf('OutputHistogram_Staircase_cp%0.2f_sgmInput%2.2f',c,sgmInput);
    for level = 1:length(boost)
        Vinput = [c1, c2]*scale;%*boost(level);
        a = a0;
        b = b0;
        w = w0;%*boost(level);
        iSTDP = boost(level);
        %initialvals = (a-1)/iSTDP./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        initialvals = (a0-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0];
        [choice, rt, R, G, I, Vcourse] = LDDM_STDP_v_RndInput(Vprior, Vinput, iSTDP, w, a, b,...
            sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        subplot(2,1,level); hold on;
        histogram(R(1000:end,1),40);
        histogram(R(1000:end,2),40);
        xlim([12,40]);
    end
    savefigs(h, filename, plotdir,fontsize-5, [3,5]);
end

%% change iSTDP level
boost = 1:.4:4;
sims = 102400;
% dt = .001; % second
% tauR = .1; % second
% tauG = .1; % second
% tauI = .1; % second
% Tau = [tauR, tauG, tauI];

w = w0;
a = a0*eye(2);
b = b0*eye(2);

% eqlb = 32;
% scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0(1,1)).*eqlb;
task = 'RT';
predur = 0;
presentt = dt;
triggert = dt;
dur = 4.5;
stimdur = Inf;
thresh = 70;
stoprule = 1;
filename = sprintf('LDDM_%s_iSTDP%.1fto%.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',task,min(boost),max(boost),a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
        %iSTDP = boost(level);
        initialvals = (a0-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        [V1, iSTDP] = meshgrid(scale*(1+cp), boost);
        [V2, iSTDP] = meshgrid(scale*(1-cp), boost);
        Vinput.V1 = V1; %scale*[1+cp, 1-cp];%*boost(level);
        Vinput.V2 = V2;
        Vprior.V1 = ones(size(V1))*scale;
        Vprior.V2 = ones(size(V2))*scale;
        [rt, choice, ~] = LDDM_STDP_v_Rndinput_GPU(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC = gather(mean(2-squeeze(choice),3,'omitnan'));
        meanRT = gather(mean(squeeze(rt),3,'omitnan'));
    %end
    save(output,'ACC','meanRT');
else
    load(output);
end
mygray = [.9:-.9/numel(boost):0];
h = figure;
subplot(2,1,1); hold on;
for level = 1:numel(boost)
    plot(cp,ACC(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('% correct');
lgd = legend(cellstr(string(boost)),...
    'Location','SouthEast', 'FontSize', fontsize-8, 'Box','off');
title(lgd,'iSTDP','FontSize',fontsize-8);
savefigs(h, ['Choice_RT_' filename], plotdir,fontsize, flip(aspect1));
subplot(2,1,2); hold on;
for level = 1:numel(boost)
    plot(cp,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, flip(aspect1)*1.2);

%% ACC as a function of Input noise, two iSTDP levels
w = w0;
a = a0*eye(2);
b = b0*eye(2);
boost = [1, 4];
sgmInputvec = [0:.1:1];
filename = sprintf('LDDM_%s_iSTDP%.1f_%.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3fto%.3f_sims%i',task,min(boost),max(boost),a0,b0,sgm,min(sgmInput),max(sgmInput),sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [32]'/1000;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    for sgi = 1:numel(sgmInputvec)
        fprintf('sgmInput %.3f\n',sgmInput);
        sgmInput = sgmInputvec(sgi);
        initialvals = (a0-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        [V1, iSTDP] = meshgrid(scale*(1+cp), boost);
        [V2, iSTDP] = meshgrid(scale*(1-cp), boost);
        Vinput.V1 = V1; %scale*[1+cp, 1-cp];%*boost(level);
        Vinput.V2 = V2;
        Vprior.V1 = ones(size(V1))*scale;
        Vprior.V2 = ones(size(V2))*scale;
        [rt, choice, ~] = LDDM_STDP_v_Rndinput_GPU(Vprior, Vinput, iSTDP, w, a, b,...
            sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC(sgi,:) = gather(mean(2-squeeze(choice),2,'omitnan'));
        meanRT(sgi,:) = gather(mean(squeeze(rt),2,'omitnan'));
    end
    save(output,'ACC','meanRT');
else
    load(output);
end
%
h = figure; hold on;
plot(sgmInputvec, ACC(:,1));
plot(sgmInputvec, ACC(:,2));
legend({'Baseline','iSTDP'},'Box','off','Location','northeast');
xlabel('Input noise (a.u.)');
ylabel('Accuracy');
savefigs(h, ['ACCoversgmInput_', filename], plotdir,fontsize-5, [3,3]);

%% largest amplitude and iSTDP 
boost = 1:.4:4;
sgm = 0;
sgmInput = 0;
triggert = Inf;
stimdur = Inf;
dur = 4.5;
Rstore = [];
h = figure;
filename = 'R dynamic over iSTDP';
for level = 1:length(boost)
    Vinput = [c1, c2]*scale;
    a = a0;%*boost(level);
    b = b0;
    w = w0;%*boost(level);
    iSTDP = boost(level);
    %initialvals = (a-1)/iSTDP./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
    initialvals = (a-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0];
    [choice, rt, R, G, I, Vcourse] = LDDM_STDP_v_RndInput(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    Rstore(level) = R(end,1);
    subplot(2,4,level);
    plot(R);
    ylim([10,40]);
    xlabel('Time (a.u.)');
    ylabel('Firing rates (Hz)');
end
savefigs(h, [filename], plotdir,fontsize-5, [6,4]);


h = figure;
filename = 'SteadyFR_iSTDP';
plot(boost,Rstore);
xlabel('iSTDP');
ylabel('Max Firing rates (Hz)');
savefigs(h, [filename], plotdir,fontsize-5, [3,3]);

