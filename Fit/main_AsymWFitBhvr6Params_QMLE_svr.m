addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../RecurrentModel/bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = '../../RecurrentModel/Fit/Rslts/AsymWFitBhvr6Params_QMLE_GPU';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
%%
% Take data from Roitman & Shadlen, 2002
dataDynmc = load('../../RecurrentModel/Fit/Data/Data.mat');
dataBhvr = LoadRoitmanData('../../RecurrentModel/RoitmanDataCode');
% Fix random seed for reproducibility
% rng(1);
% change random seed
t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7;
num2str(t);
rng(t);
% Define optimization starting point and bounds
%     a,    w1,  noise, scale, tauR, tauG
LB = [0    0.001   .1    .1*256, [.001,.001]];
UB = [70   1	100  20*256, [1, 1]];
PLB = [15  .1	2    .5*256, [.01 .01]];
PUB = [60   .4	20   8*256, [.2 .2]];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
nLLfun = @(params) AsymWFitBhvr6Params_QMLE_GPU(params, dataBhvr);
[fvalbest,~,~] = nLLfun(x0)
fprintf('test succeeded\n');
% change starting points
Collect = [];
parfor i = 1:myCluster.NumWorkers*8
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

if 0
%% hand tuning
% addpath('../../CoreFunctions/');
% addpath('./SvrCode/');
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Fit');
cd('G:\My Drive\LDDM\Fit');
out_dir = './Rslts/AsymWFitBhvr6Params_QMLE_GPU';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
dataDynmc = load('./Data/Data.mat');
dataBhvr = LoadRoitmanData('../RoitmanDataCode');

randseed = 105874882; % 110234755;%
rng(randseed);
% a, w1, noise, scale, tauR, tauG, nLL
% params = [25.350783	0.401201	4.272408	99.838832	0.001936	0.160071	16566.43769];
params = [36.538884	0.134803	7.401239	284.382666	0.052278	0.231087	16675.6422];
name = sprintf('a%2.1f_w%1.1f_noise%2.1f_scl%2.1f_tau%1.3f_%1.3f',params(1:6));

% simulation
if ~exist(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),'file')
    tic;
    [nLL, Chi2, BIC, AIC, rtmat, choicemat] = AsymWFitBhvr6Params_QMLE_GPU(params, dataBhvr);
    save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
        'rtmat','choicemat','params','nLL','Chi2','AIC','BIC');
    toc
else
    load(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)));
end

%% Example dynamics
lwd = 1;
mksz = 3;
fontsize = 11;
% a, w1, noise, scale, tauR, tauG, nLL
% params = [25.350783	0.401201	4.272408	99.838832	0.001936	0.160071	16566.43769];
params = [36.538884	0.134803	7.401239	284.382666	0.052278	0.231087	16675.6422];
simname = sprintf('AsymW_a%2.1f_w%1.1f_noise%2.1f_scl%2.1f_tau%1.3f_%1.3f',params(1:6));
a = params(1)*eye(2);
w1 = params(2);
sgm =  params(3)/5;
scale = params(4);
ndt = .03; % sec, resort to the saccade side, the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
tgap = 0.09; % initial dip for 90ms after stimuli onset
tauR = params(5);
tauG = params(6);
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
c1 = (1 + Cohr)';
c2 = (1 - Cohr)';
cplist = [c1, c2];
mygray = flip(gray(length(cplist)));
predur = 0;
presentt = 0;
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; %70.8399; % mean(max(m_mr1cD))+1; 
stimdur = dur;
stoprule = 1;
w = [w1 1; 1 w1];
Rstar = 32; % ~ 32 Hz according to Roitman and Shadlen's data, at the bottom of initial dip
initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar;0, 0];
eqlb = Rstar; % set equilibrium value before task as R^*
Vprior = [1, 1]*(2*mean(w,'all')*eqlb.^2 + (1-a(1)).*eqlb);
Tau = [tauR tauG];
h = figure; 
% subplot(2,1,1); 
hold on;
filename = sprintf('%s',simname);
randseed = 75245522;
rng(randseed);
for vi = 2:6
    Vinput = cplist(vi,:)*scale;
    [R, G, rt, choice] = AsymW(Vprior, Vinput, w, a, sgm, Tau,...
    predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);

    lgd2(vi-1) = plot(R(:,2), 'k-.', 'Color', mygray(vi,:), 'LineWidth',lwd);
    lgd1(vi-1) = plot(R(R(:,1)<=thresh,1), 'k-', 'Color', mygray(vi,:), 'LineWidth',lwd);
end
plot([.2, 1.2]/dt,[thresh,thresh], 'k-');
text(600,thresh*1.1,'threshold');
yticks([0,32,70]);
yticklabels({'0','32','70'});
% ylim([0,74]);
ylabel('Activity (Hz)');
xticks([0, 500, 1000, 1500]);
xticklabels({'0','.5','1.0','1.5'});
xlim([-50, 1200]);
xlabel('Time (s)');
% savefigs(h, filename, plot_dir, fontsize, [2 1.5]);

% 
% subplot(2,1,2); hold on;
% for vi = 2:6
%     Vinput = cplist(vi,:);
%     [R, G, ~, ~, ~] = AsymW(Vinput, w, a, sgm, Tau,...
%         dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%     lgd2(vi-1) = plot(diff(R(:,2)), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
%     lgd1(vi-1) = plot(diff(R(:,1)), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
% end
% ylabel('Time derivative');
% xticks([0, 500, 1000, 1500]);
% xticklabels({'0','.5','1.0','1.5'});
% xlim([-50, 1800]);
% yticks([-.1,0,.3]);
% ylim([-.1, .35]);
% xlabel('Time (s)');
% savefigs(h, filename, plot_dir, fontsize, [1.8 1.4]);

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
addpath('../RoitmanDataCode');
ColumnNames608
load T1RT.mat;
x(:,R_RT) = x(:,R_RT)/1000;
cohlist = unique(x(:,R_COH));
maxrt = max(x(:,R_RT));
minrt = min(x(:,R_RT));
segrt = maxrt - minrt;
bins = 29;
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
% saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.eps',name)),'epsc2');
%% aggregated RT & ACC
lwd = 1;
mksz = 3;
fontsize = 11;
Cohr = [0 32 64 128 256 512]/1000;
cplist = Cohr*100;
cplist(1) = 1.1;
h = figure;
filename = sprintf('RT&ACC_%s',name);
subplot(2,1,1);
hold on;
plot(cplist, accr*100, 'xk', 'MarkerSize', mksz+1);
plot(cplist, acc*100,'-k','LineWidth',lwd);
ylim([.45,1]*100);
yticks([50,100]);
xlim([1,100]);
xticks([1,10,100]);
xticklabels({'0','10','100'});
ylabel('Correct (%)');
xlabel('Input coherence (%)');
set(gca, 'XScale', 'log');
% legend({'data','model'},'NumColumns',1,'Location','SouthEastOut','FontSize',fontsize-2);
% legend('boxoff');
savefigs(h,filename,plot_dir,fontsize,[2,3.0]);

subplot(2,1,2);
hold on;
lg1 = plot(cplist, meanrtcr, '.k', 'MarkerSize', mksz*3);
lg2 = plot(cplist, meanrtc, '-k','LineWidth',lwd);
lg3 = plot(cplist, meanrtwr, 'ok', 'MarkerSize', mksz);
lg4 = plot(cplist, meanrtw, '--k','LineWidth',lwd);
xlim([1,100]);
xticks([1,10,100]);
xticklabels({'0','10','100'});
yticks([.4,1]);
ylim([.4, 1]);
ylabel('RT (secs)');
xlabel('Input coherence (%)');
set(gca, 'XScale', 'log');
% lgd = legend([lg3,lg1,lg4,lg2],{'','','error','correct'},'NumColumns',2,'Location','SoutheastOut','FontSize',fontsize-2);
% legend('boxoff');
savefigs(h,filename,plot_dir,fontsize,[2,3.0]);
%% Q-Q plot for reaction time and choice
lwd = 1.0;
mksz = 3;
fontsize = 11;
x = dataBhvr.proportionmat;
y = dataBhvr.q;
qntls = dataBhvr.qntls;
h = figure; hold on;
for vi = 1:length(x)
    xc = x(vi)*ones(size(y(:,1,vi)));
    xw = 1 - x(vi)*ones(size(y(:,2,vi)));
    plot(xc,y(:,1,vi),'gx','MarkerSize',mksz+1,'LineWidth',lwd);
    plot(xw,y(:,2,vi),'rx','MarkerSize',mksz+1,'LineWidth',lwd);
    % fitted value
    En(vi) = numel(rtmat(:,vi));
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    xr = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
    q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
    
end
for qi = 1:size(q,1)
    xq = [flip(1-x), x]';
    plot(xq,[squeeze(flip(q(qi,2,:)));squeeze(q(qi,1,:))],'k-o','MarkerSize',mksz,'LineWidth',lwd/2);
end
xlim([-.05,1.05]);
ylim([.2, 1.4]);
yticks([.2:.4:1.4]);
xlabel('Proportion');
ylabel('RT (s)');
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 4 5];
filename = sprintf('Q-QPlot_%s',name);
% saveas(h,fullfile(plot_dir,sprintf('Q-QPlot_%s.eps',name)),'epsc2');
savefigs(h, filename, plot_dir, fontsize, [2 2.5]);
%% the original space of QMLE
acc = dataBhvr.proportionmat;
ON = dataBhvr.ON;
OP = dataBhvr.OP;
Oq = dataBhvr.q;
qntls = dataBhvr.qntls;
expd_qntls = [0, qntls, 1];
P = expd_qntls(2:end) - expd_qntls(1:end-1); % proportion in each bin
En = [];
f = [];
h = figure;
for vi = 1:length(acc)
    subplot(6,1,vi); hold on;
    x = 1:10;%Oq(:,1,vi);
    y = ON(:,1,vi).*log(OP(:,1,vi));
    plot(x,y,'gx');
    x = 1:10;%Oq(:,2,vi);
    y = ON(:,2,vi).*log(OP(:,2,vi));
    plot(x,y,'rx');
    set(gca,'YScale','log');
    xlim([0, 11]);
    %ylim(-[1000,1]);
    
    En(vi) = numel(rtmat(:,vi));
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    if ~all(isnan(Oq(:,1,vi)))
        tmp = histogram(RT_corr, [0; Oq(:,1,vi); Inf], 'Visible',0);
        EN(:,1,vi) = tmp.Values;
    else
        EN(:,1,vi) = NaN(numel(Oq(:,1,vi))+1,1);
    end
    if ~all(isnan(Oq(:,2,vi)))
        tmp = histogram(RT_wro, [0; Oq(:,2,vi); Inf], 'Visible',0);
        EN(:,2,vi) = tmp.Values;
    else
        EN(:,2,vi) =  NaN(numel(Oq(:,2,vi))+1,1);
    end
    f(:,:,vi) = log((EN(:,:,vi)/En(vi)));
    plot(x,ON(:,1,vi).*f(:,1,vi),'g-');
    plot(x,ON(:,2,vi).*f(:,2,vi),'r-');
end
saveas(h,fullfile(plot_dir,sprintf('QMLE_Plot_%s.eps',name)),'epsc2');
% %% Sumed Quantile loglikelihodd over coherence
% h = figure;
% for vi = 1:length(acc)
%     Obar(vi) = sum([ON(:,1,vi).*log(OP(:,1,vi)); ON(:,2,vi).*log(OP(:,2,vi))],'omitnan');
%     Ebar(vi) = sum([ON(:,1,vi).*(f(:,1,vi)); ON(:,2,vi).*(f(:,2,vi))],'omitnan');
% end
% bar([Obar; Ebar]','grouped');
% saveas(h,fullfile(plot_dir,sprintf('SumLL_Plot_%s.eps',name)),'epsc2');
% %% proportion at each quantile
% h = figure;
% for vi = 1:length(acc)
%     subplot(6,1,vi); hold on;
%     x = 1:10;
%     plot(x, OP(:,1,vi),'gx');
%     plot(x, -OP(:,2,vi),'rx');
%     plot(x, EN(:,1,vi)/En(vi),'g-');
%     plot(x, -EN(:,2,vi)/En(vi),'r-');
%     ylim([-.1, .2]);
% end
% saveas(h,fullfile(plot_dir,sprintf('Proportion_Plot_%s.eps',name)),'epsc2');
end
