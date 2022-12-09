%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(fullfile(Homedir,'Documents','LDDM','Froemke/SST_STDPadtv/'));
cd('G:\My Drive\LDDM\Froemke\SST_STDPadtv');
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
a0 = 15;
b0 = 2.8;
predur = 0;
presentt = 0;% dt*200;
stimdur = Inf;
triggert = presentt;% + dt;% + .8;
dur = 4.5; % second

thresh = 70; % Hz
stoprule = 1;

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0).*eqlb;

c = [.256]'; % Coherence
c1 = [1 + c];
c2 = [1 - c];
Vprior = [1,1]*scale;

boost = [0, 36];
%% Dynamics
sgm = 0;
h = figure;
hold on;
mycl = jet(length(boost)*10);
for level = 1:length(boost)
    Vinput = [c1, c2]*scale;%*boost(level);
    a = a0;%*boost(level);
    b = b0;
    w = w0;%*boost(level);
    iSTDP = boost(level);
    % initialvals = (a-1-iSTDP)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
    initialvals = (a-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0];
    [choice, rt, R, G, I] = LDDM_STDPadtv(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    plot(R(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    plot(R(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(G(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(G(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(I(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(I(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
end
ylim([0,100]);
xlabel('Time (a.u.)');
xticks([presentt/dt]);
xticklabels({'stimuli','action'});
ylabel('Firing Rates (Hz)');
lgd = legend(' ', ' ', 'R1', 'R2',...
    'Location','Northeast','NumColumns',2, 'FontSize', fontsize-8, 'Box','off');
title(lgd, "Baseline                  iSTDP      .");
filename = sprintf('timeCourse_RT_R_%1.1f_%1.1f_sgm%2.2f',boost, sgm);
savefigs(h, filename, plotdir,fontsize, aspect1);
%% behavior at two levels of boost
boost = [0, 35];
sims = 10000;
sgmInput = 0; % 1/3
sgm = 10; % 6.8/2;
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];

w0 = ones(2);
a0 = 15*eye(2);
b0 = 1.7*eye(2);

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0(1,1)).*eqlb;
predur = 0;
presentt = dt;
triggert = dt;
dur = 6;
stimdur = dur;
thresh = 70;
stoprule = 1;
filename = sprintf('LDDM_%.1f_%.1f_a%1.0f_b%1.0f_%2.0fsgm%1.1fsinpt%2.0f_sims%i',boost,a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
ACC = [];
meanRT = [];
if ~exist(output,'file')
    for level = 1:length(boost)
        a = a0;%*boost(level);
        b = b0;
        w = w0; %.*boost(level);
        iSTDP = boost(level);
        initialvals = (a(1,1)-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        Vinput = scale*[1+cp, 1-cp];%*boost(level);
        Vprior = ones(size(Vinput))*scale;
        [rt, choice, ~] = LDDM_STDP_GPU(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC(level,:,:) = gather(mean(2-squeeze(choice),2,'omitnan'));
        meanRT(level,:,:) = gather(mean(squeeze(rt),2,'omitnan'));
    end
    save(output,'ACC','meanRT');
else
    load(output);
end
h = figure;
filename = sprintf('Choice_RT_%1.1f_%1.1f_sgm%2.2f',boost, sgm);
subplot(2,1,1); hold on;
for level = 1:length(boost)
    plot(cp,ACC(level,:),'.-');
end
set(gca,'XScale','log');
ylabel('% correct');
lgd = legend('Baseline', 'iSTDP',...
    'Location','SouthEast', 'FontSize', fontsize-8, 'Box','off');
savefigs(h, filename, plotdir,fontsize, flip(aspect1));
subplot(2,1,2); hold on;
for level = 1:length(boost)
    plot(cp,meanRT(level,:),'.-');
end
set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, flip(aspect1));


%% change iSTDP level
boost = 0:5:40;
sims = 10000;
sgmInput = 0; % 1/3
sgm = 10; %6.8/2;
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];

w0 = ones(2);
a0 = 15*eye(2);
b0 = 1.7*eye(2);

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0(1,1)).*eqlb;

presentt = dt;
triggert = dt;
dur = 6;
stimdur = dur;
thresh = 70;
stoprule = 1;
filename = sprintf('LDDM_%.1fto%.1f_a%1.0f_b%1.0f_%2.0fsgm%1.1fsinpt%2.0f_sims%i',min(boost),max(boost),a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
        a = a0;%*boost(level);
        b = b0;
        w = w0; %.*boost(level);
        %iSTDP = boost(level);
        initialvals = (a(1,1)-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        [V1, iSTDP] = meshgrid(scale*(1+cp), boost);
        [V2, iSTDP] = meshgrid(scale*(1-cp), boost);
        Vinput.V1 = V1; %scale*[1+cp, 1-cp];%*boost(level);
        Vinput.V2 = V2;
        Vprior.V1 = ones(size(V1))*scale;
        Vprior.V2 = ones(size(V2))*scale;
        [rt, choice, ~] = LDDM_STDPadtv_GPU(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC = gather(mean(2-squeeze(choice),3,'omitnan'));
        meanRT = gather(mean(squeeze(rt),3,'omitnan'));
    %end
    save(output,'ACC','meanRT');
else
    load(output);
end
mygray = [.9:-.9/numel(boost):0];
h = figure;
filename = sprintf('Choice_RT_%1.1fto%1.1f_sgm%2.2f',min(boost),max(boost), sgm);
subplot(2,1,1); hold on;
for level = 1:numel(boost)
    plot(cp,ACC(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('% correct');
lgd = legend(cellstr(string(boost)),...
    'Location','SouthEast', 'FontSize', fontsize-8, 'Box','off');
title(lgd,'iSTDP','FontSize',fontsize-8);
savefigs(h, filename, plotdir,fontsize, flip(aspect1));
subplot(2,1,2); hold on;
for level = 1:numel(boost)
    plot(cp,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, flip(aspect1)*1.2);

