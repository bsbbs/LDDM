%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(fullfile(Homedir,'Documents','LDDM','Froemke/SST_STDP/'));
cd('G:\My Drive\LDDM\Froemke\SST_STDP');
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

presentt = dt*200;
stimdur = Inf;
triggert = presentt + dt + .8;
dur = 3.5; % second

thresh = 70; % Hz
stoprule = 1;

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0).*eqlb;

c = [.256]'; % Coherence
c1 = [1 + c];
c2 = [1 - c];

boost = [1, 2];
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
    initialvals = (a-1)/iSTDP./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
%     [choice, rt, R, G, I] = LDDM(Vinput, w, a, b, sgm, Tau, dur,...
%         dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    [choice, rt, R, G, I] = LDDM_SST(Vinput, iSTDP, w, a, b, sgm, Tau, dur,...
        dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    plot(R(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
    plot(R(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(G(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
% plot(G(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
%     plot(I(:,1), 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
%     plot(I(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(level-1)*10,:));
end
ylim([0,100]);
xlabel('Time (a.u.)');
xticks([presentt/dt]);
xticklabels({'stimuli','action'});
ylabel('Firing Rates (Hz)');
lgd = legend(' ', ' ', 'R1', 'R2',...
    'Location','Northeast','NumColumns',2, 'FontSize', fontsize-8, 'Box','off');
title(lgd, "Baseline                  iSTDP      .");
filename = sprintf('timeCourse_FD_R_%1.1f_%1.1f_sgm%2.2f',boost, sgm);
savefigs(h, filename, plotdir,fontsize, aspect1);
%% behavior at two levels of boost
boost = [1, 2];
sims = 10000;
sgmInput = 0; % 1/3
sgm = 10; % 6.8/2;
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];

w0 = ones(2);
a0 = 15;
b0 = 1.7;

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0).*eqlb;

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
        initialvals = (a-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        Vinput = scale*[1+cp, 1-cp];%*boost(level);
        [rt, choice, ~] = LcDsInhbt_GPU_STDP(Vinput, iSTDP, w, eye(2)*a, eye(2)*b, sgm, Tau,...
            dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
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
boost = 1:.1:2;
sims = 10000;
sgmInput = 0; % 1/3
sgm = 10; %6.8/2;
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];

w0 = ones(2);
a0 = 15;
b0 = 1.7;

eqlb = 32;
scale = 2*mean(w0,'all')*eqlb.^2 + (1-a0).*eqlb;

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
if ~exist(output,'file')
    %for level = 1:length(boost)
        a = a0;%*boost(level);
        b = b0;
        w = w0; %.*boost(level);
        %iSTDP = boost(level);
        initialvals = (a-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0]; % for R, G, and I variables, will be further scaled
        [V1, iSTDP] = meshgrid(scale*(1+cp), boost);
        [V2, iSTDP] = meshgrid(scale*(1-cp), boost);
        Vinput.V1 = V1; %scale*[1+cp, 1-cp];%*boost(level);
        Vinput.V2 = V2;
        [rt, choice, ~] = LcDsInhbt_GPU_STDP(Vinput, iSTDP, w, eye(2)*a, eye(2)*b, sgm, Tau,...
            dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC = gather(mean(2-squeeze(choice),3,'omitnan'));
        meanRT = gather(mean(squeeze(rt),3,'omitnan'));
    %end
    save(output,'ACC','meanRT');
else
    load(output);
end
mygray = [.9:-.9/numel(boost):0];
h = figure;
filename = sprintf('Choice_RT_%1.1fto%1.1f_sgm%2.2f',boost, sgm);
subplot(2,1,1); hold on;
for level = 1:numel(boost)
    plot(cp,ACC(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('% correct');
% lgd = legend('Baseline', 'iSTDP',...
%     'Location','SouthEast', 'FontSize', fontsize-8, 'Box','off');
savefigs(h, filename, plotdir,fontsize, flip(aspect1));
subplot(2,1,2); hold on;
for level = 1:numel(boost)
    plot(cp,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1]);
end
set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, flip(aspect1));

