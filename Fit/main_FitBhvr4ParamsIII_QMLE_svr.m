addpath('../');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = './Rslts/FitBhvr4ParamsIII_QMLE_SvrGPU';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
%%
% Take data from Roitman & Shadlen, 2002
dataDynmc = load('./Data/Data.mat');
dataBhvr = LoadRoitmanData('../RoitmanDataCode');
% Fix random seed for reproducibility
% rng(1);
% change random seed
t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7;
num2str(t);
rng(t);
% Define optimization starting point and bounds
%     a,    b, noise, scale
LB = [0    0.1   .1    .1*256];
UB = [70   2	100  20*256];
PLB = [15  .7	5    1*256];
PUB = [60   1.2	40   8*256];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
nLLfun = @(params) LcDFitBhvr4ParamsIII_QMLE_GPU(params, dataDynmc, dataBhvr);
[fvalbest,~,~] = nLLfun(x0)
fprintf('test succeeded\n');
% change starting points
Collect = [];
parfor i = 1:myCluster.NumWorkers*4
    !ping -c 1 www.amazon.com
    t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7 + i*10^5;
    %num2str(t);
    rng(t);
    
    % Randomize initial starting point inside plausible box
    x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;
    dlmwrite(fullfile(out_dir,'x0List.txt'),[sortNum, i, t, x0],'delimiter','\t','precision','%.6f','-append');
    % fit
    options = bads('defaults');     % Default options
    options.Display = 'iter';
    % options.UncertaintyHandling = false;    % Function is deterministic
    options.UncertaintyHandling = true;    % Function is deterministic
    [xest,fval,~,output] = bads(nLLfun,x0,LB,UB,PLB,PUB,options);
    %     xest
    %     fval
    %     output
    dlmwrite(fullfile(out_dir,'RsltList.txt'),[sortNum, i, t, xest fval],'delimiter','\t','precision','%.6f','-append');
    %     if fval < fvalbest
    %         xbest = xest;
    %         fvalbest = fval;
    %         outputbest = output;
    %     end
    % save(fullfile(out_dir,sprintf('./Rslts%i_%i.mat', sortNum, i)),'xest','fval','output');
    Collect(i).rndseed = t;
    Collect(i).x0 = x0;
    Collect(i).xest = xest;
    Collect(i).fval = fval;
    Collect(i).output = output;
    
end
t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7 + i*10^5;
save(fullfile(out_dir,sprintf('CollectRslts%i.mat',t)),'Collect');

%% hand tuning
addpath('../../CoreFunctions/');
addpath('./SvrCode/');
out_dir = './Rslts/FitBhvr4ParamsII_QMLE_SvrGPU';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
dataDynmc = load('./Data/Data.mat');
dataBhvr = LoadRoitmanData('../../RoitmanDataCode');
randseed = 22967888;
rng(randseed);
% a, b, noise, scale
params = [0.000017	0.730205	44.129237	3166.970422];
w = ones(2);
a = eye(2)*params(1); %params(3);
b = eye(2)*params(2);
sgm = params(3);%params(4);
Tau = [.1,.1,.1];
ndt = .19 + .03; % sec, 90ms after stimuli onset, activities show ~100 ms of dip and recovery;
% resort to the saccade side, the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
tgap = 0; % changed for this version to move the fitting begin after the time point of recovery
scale = params(4);

deduction = .1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
predur = 0;
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; %72.1491;
stimdur = dur;
stoprule = 1;


V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;
load('./Data/Data.mat');
m_mr1c = m_mr1c';
m_mr2c = m_mr2c';
m_mr1cD = m_mr1cD';
m_mr2cD = m_mr2cD';
dot_ax = dot_ax';
sac_ax = sac_ax';
Rstar = 42; % mean(mean([m_mr1c(1:5,:); m_mr2c(1:5,:)])); % 34.6004
initialvals = [Rstar,Rstar; 2*Rstar,2*Rstar; 0,0];
Vprior = ones(6,2)*((1-a(1,1))*Rstar + 2*Rstar^2);

% dot_ax = [0:20:1000];
%sac_ax = [-1000:20:300];
% simulation
tic;
[rtmat, choicemat, argmaxR, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LcD_Dynmc_Trim_GPU(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, tgap, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
rtmat = gather(squeeze(rtmat))' + ndt;
choicemat = gather(squeeze(choicemat))';
name = sprintf('a%2.1f_b%1.1f_tau%1.2f_%1.2f_%1.2fscl%2.1f',a(1,1),b(1,1),Tau,scale);
save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
    'rtmat','choicemat','sm_mr1c','sm_mr2c','sm_mr1cD','sm_mr2cD','a','b','w','sgm','Tau','ndt','scale');
toc
load(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)));
%% plot time course
h = figure;
subplot(1,2,1);hold on;
colvec = flip({[218,166,109]/256,[155 110 139]/256,'#32716d','#af554d','#708d57','#3b5d64'});
for ci = 1:6
    lg(ci) = plot(dot_ax/1000, sm_mr1c(:,ci),'Color',colvec{ci},'LineWidth',2);
    plot(dot_ax/1000, sm_mr2c(:,ci),'--','Color',colvec{ci},'LineWidth',2);
end
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
% ylim([20,60]);
ylim([10,70.5]);
ylabel('Firing rate (sp/s)');
xlabel('Time (secs)');
xticks(0:.2:.8);
set(gca,'FontSize',16);
legend(lg,{'0','3.2','6.4','12.8','25.6','51.2'},'Location','northwest','FontSize',11);
subplot(1,2,2);hold on;
plot([0,0],[20,71],'-k');
for ci = 1:6
    plot(sac_ax/1000, sm_mr1cD(:,ci),'Color',colvec{ci},'LineWidth',2);
    plot(sac_ax/1000, sm_mr2cD(:,ci),'--','Color',colvec{ci},'LineWidth',2);
end
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
yticks([]);
set(gca,'ycolor',[1 1 1]);
ylim([10,70.5]);
set(gca,'FontSize',16);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile(plot_dir,sprintf('FittedTimeCourse_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('FittedTimeCourse_%s.eps',name)),'epsc2');
%% plot firing rates at position a,b,c,d 
h = figure;
subplot(2,1,1);hold on;
x = Cohr*100;
y = sm_mr1c(20,:);
plot(x, y,'k.','MarkerSize',16);
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-');
y = sm_mr2c(20,:);
plot(x, y,'k.','MarkerSize',16);
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-');
ylim([10,45]);
xlim([-4,55.2]);
yticks([10:10:40]);
xticks([0:10:50]);
xticklabels({});
ylabel('Firing rates (sp/s)');
set(gca,'FontSize',14);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;

subplot(2,1,2);hold on;
y = sm_mr1cD(end-17,:);
plot(x, y,'k.','MarkerSize',16);
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-');
y = sm_mr2cD(end-17,:);
plot(x, y,'k.','MarkerSize',16);
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-');
ylim([0,60]);
xlim([-4,55.2]);
yticks([0:10:70]);
xticks([0:10:50]);
xlabel('Input strength (% coh)');
ylabel('Firing rates (sp/s)');
set(gca,'FontSize',14);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 4];
saveas(h,fullfile(plot_dir,sprintf('abcd_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('abcd_%s.eps',name)),'epsc2');

%% plot RT distribution - fitted
rate = length(rtmat)/1024;
maxrt = max(max(rtmat));
minrt = min(min(rtmat));
%segrt = maxrt - minrt;
bank1 = [];
bank2 = [];
acc = [];
meanrtc = [];
meanrtw = [];
for ii = 1:6
    gap = (dataBhvr.rtrange(ii,2) - dataBhvr.rtrange(ii,1))/dataBhvr.numbins;
    %gap = 1.757/60;
    BinEdge = [minrt:gap:(maxrt+gap)];
    hg = histogram(rtmat(choicemat(:,ii)==1,ii),'BinEdges',BinEdge);
    bank1{ii} = hg.Values/rate;
    hg = histogram(rtmat(choicemat(:,ii)==2,ii),'BinEdges',BinEdge);
    bank2{ii}= hg.Values/rate;
    BinMiddle{ii} = hg.BinEdges(1:end-1) + hg.BinWidth/2;
    acc(ii) = sum(choicemat(:,ii)==1)/(sum(choicemat(:,ii)==1) + sum(choicemat(:,ii)==2));
    meanrtc(ii) = mean(rtmat(choicemat(:,ii)==1,ii));
    meanrtw(ii) = mean(rtmat(choicemat(:,ii)==2,ii));
end
% loading Roitman's data
addpath('../../RoitmanDataCode');
ColumnNames608
load T1RT.mat;
x(:,R_RT) = x(:,R_RT)/1000;
cohlist = unique(x(:,R_COH));
maxrt = max(x(:,R_RT));
minrt = min(x(:,R_RT));
segrt = maxrt - minrt;
bins = 30;
BinEdge = [minrt:segrt/bins:maxrt];
bank1r = [];
bank2r = [];
accr = [];
meanrtcr = [];
meanrtwr = [];
for i = 1:length(cohlist)
    Lcoh = x(:,R_COH)==cohlist(i);
    if i == 1
        Dir1 = x(:,R_TRG) == 1;
        Dir2 = x(:,R_TRG) == 2;
        RT_corr = x(Lcoh & Dir1,R_RT);
        RT_wro = x(Lcoh & Dir2, R_RT);
    else
        Corr = x(:,R_DIR) == x(:,R_TRG);
        Wro = x(:,R_DIR) ~= x(:,R_TRG);
        RT_corr = x(Lcoh & Corr,R_RT);
        RT_wro = x(Lcoh & Wro, R_RT);
    end
    accr(i) = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    meanrtcr(i) = mean(RT_corr);
    meanrtwr(i) = mean(RT_wro);
    hg = histogram(RT_corr,'BinEdges',BinEdge);
    bank1r(:,i) = hg.Values;
    if ~isempty(RT_wro)
        hg = histogram(RT_wro,'BinEdges',BinEdge);
        bank2r(:,i) = hg.Values;
    else
        bank2r(:,i) = zeros(1,bins);
    end
end
BinMiddler = hg.BinEdges(1:end-1) + hg.BinWidth/2;
h = figure;
for ii = 1:6
    subplot(6,1,ii);
    %bar(BinMiddler,bank1r(:,ii),'FaceColor','#0072BD','EdgeAlpha',0);
    bar(dataBhvr.bincenter(ii,1:30),dataBhvr.histmat(ii,1:30)*1024,'FaceColor','#0072BD','EdgeAlpha',0);
    hold on;
    %bar(BinMiddler,-bank2r(:,ii),'FaceColor','#D95319','EdgeAlpha',0);
    %plot(BinMiddle{ii},bank1{ii},'c','LineWidth',1.5);
    %plot(BinMiddle{ii},-bank2{ii},'r','LineWidth',1.5);
    bar(dataBhvr.bincenter(ii,1:30),-dataBhvr.histmat(ii,31:60)*1024,'FaceColor','#D95319','EdgeAlpha',0,'EdgeColor','none');
    plot(BinMiddle{ii},bank1{ii},'c','LineWidth',2);
    plot(BinMiddle{ii},-bank2{ii},'m','LineWidth',2);
    if ii == 7
        legend({'','','Correct','Error'},'NumColumns',2,'Location','North');
        legend('boxoff');
    end
    ylim([-60,100]);
    yticks([-50:50:100]);
    yticklabels({'50','0','50','100'});
    xlim([100 1762]/1000);
    xticks([.5,1.0,1.5]);
    if ii == 6
        xticklabels({'.5','1.0','1.5'});
        xlabel('Reaction time (secs)');
    else
        xticklabels({});
    end
    if ii == 1
        ylabel('Frequency');
    end
    % title(sprintf('coherence %2.1f %%',cohlist(ii)*100));
    set(gca,'FontSize',16);
    set(gca,'TickDir','out');
    H = gca;
    H.LineWidth = 1;
    set(gca, 'box','off');
end
%set(gca,'FontSize',18);

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 3.0 10];
saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.eps',name)),'epsc2');
%% aggregated RT & ACC
cplist = Cohr*100;
h = figure;
subplot(1,2,1);
hold on;
plot(cplist, accr*100, 'xk', 'MarkerSize', 8);
plot(cplist,acc*100,'-k','LineWidth',2);
ylim([.5,1]*100);
xlim([0,100]);
ylabel('Accuracy (%)');
xlabel('Input Strength (% coh)');
set(gca, 'XScale', 'log');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
legend({'data','model'},'NumColumns',2,'Location','SouthEast','FontSize',14);
legend('boxoff');


subplot(1,2,2);
hold on;
lg1 = plot(cplist, meanrtcr, '.k', 'MarkerSize', 20);
lg2 = plot(cplist, meanrtc, '-k','LineWidth',2);
lg3 = plot(cplist, meanrtwr, 'ok', 'MarkerSize', 7);
lg4 = plot(cplist, meanrtw, '--k','LineWidth',2);
xlim([0,100]);
ylabel('Reaction time (secs)');
xlabel('Input Strength (% coh)');
set(gca, 'XScale', 'log');
set(gca,'FontSize',16);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
lgd = legend([lg3,lg1,lg4,lg2],{'','','Error','Correct'},'NumColumns',2,'Location','SouthWest','FontSize',14);
%legend({'empirical','fitted'},'NumColumns',2);
legend('boxoff');
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 9 3.0];
saveas(h,fullfile(plot_dir,sprintf('RT&ACC_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('RT&ACC_%s.eps',name)),'epsc2');



%% raw data time course
h = figure;
subplot(1,2,1);hold on;
plot(dot_ax, m_mr1c,'LineWidth',1.5);
plot(dot_ax, m_mr2c,'--','LineWidth',1.5);
set(gca,'FontSize',18);
subplot(1,2,2);hold on;
plot(sac_ax, m_mr1cD,'LineWidth',1.5);
plot(sac_ax, m_mr2cD,'--','LineWidth',1.5);
set(gca,'FontSize',18);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile(plot_dir,sprintf('Data.eps')),'epsc2');
%% disribution of fitted parameters
rslts = dlmread(fullfile(out_dir,'RsltList.txt'));
name = {'a', 'b', 'noise', 'tauR', 'tauG', 'tauI', 'ndt', 'scale', 'sigma of ll'};
h = figure;
for i = 1:9
    subplot(3,3,i);
    hist(rslts(:,i+3));
    xlabel(name{i});
end
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile(plot_dir,sprintf('FittedParamsDistribution.eps')),'epsc2');
