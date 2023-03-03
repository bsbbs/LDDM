addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../RecurrentModel/bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = '../../RecurrentModel/Fit/Rslts/WW06FitBhvr8ParamsVI_QMLE_GPU';
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
%    JNp,   JNn,  I0,      noise,   miu0,   gamma, H0, tauNMDA
LB = [.01  .01     0      .005       0      0      0     .01];
UB = [.6    .2     1      .1        120     1      12     1];
PLB = [.1  .02	  .2      .01       20      .5      1    .1];
PUB = [.4  .06	  .4      .03       40      .8      8    .2];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
sims = 1024*5;
nLLfun = @(params) WW06FitBhvr8ParamsVI_QMLE_GPU(params, dataBhvr, sims);
[fvalbest,~,~] = nLLfun(x0)
fprintf('test succeeded\n');
% change starting points
Collect = [];
parfor i = 1:myCluster.NumWorkers*8
    % !ping -c 1 www.amazon.com
    % !ping www.google.com
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
    % Homedir = 'C:\Users\Bo';
    Homedir = '~';
    addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
    addpath(fullfile(Homedir,'Documents','LDDM','utils'));
    addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
    % Glgdr = 'G:\My Drive\LDDM';
    Glgdr = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM';
    cd(fullfile(Glgdr, 'Fit'));
    RoitmanDataDir = fullfile(Homedir,'Documents','LDDM', 'Fit', 'RoitmanDataCode');
    dataDynmc = load(fullfile(RoitmanDataDir,'DynmcsData.mat'));
    dataBhvr = LoadRoitmanData(RoitmanDataDir);
    
    out_dir = fullfile(Glgdr, 'Fit', 'Rslts/WW06FitBhvr8ParamsVI_QMLE_GPU');
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    plot_dir = fullfile(out_dir,'graphics');
    if ~exist(plot_dir,'dir')
        mkdir(plot_dir);
    end
    randseed = 90033894;
    rng(randseed);
    %    JNp, JNn, I0, noise, miu0, gamma, H0, tauNMDA, nLL
    % params = [.2609, .0497, .3255, .02, 30, .64, .1]; % in the paper ww06
    params = [0.263217	0.022411	0.264743	0.070874	55.634434	0.588655	2.622382	0.167199	16586.81822]; % 16586.8182 ± 4.5167
    name = sprintf('JNp%2.1f_JNn%1.2f_I0%1.2f_noise%1.2f_miu0%2.2f_gamma%1.3f_H0%2.1f_tauS%0.2f_nLL%5.2f',params);
    
    
    %% simulation, behavior only
    sims = 10240;
    if  ~exist(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat] = WW06FitBhvr8ParamsVI_QMLE_GPU(params, dataBhvr, sims);
        num2str(nLL)
        num2str(AIC)
        num2str(BIC)
        save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
            'rtmat','choicemat','params','nLL','Chi2','BIC','AIC');
        toc
    else
        load(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)));
    end
    %% ditribution of RT and fitted line
    h = RT_Dstrbtn(dataBhvr, rtmat, choicemat, plot_dir, name);
    
    %% aggregated RT & ACC
    h = PlotmeanRTACC(RoitmanDataDir, rtmat, choicemat, plot_dir, name);
    
    %% Quantile probability plot for reaction time and choice
    h = QPP_Roitman(dataBhvr, rtmat, choicemat, plot_dir, name);
    
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
        f(f(:,1,vi) == -Inf,1,vi) = log(.00001/En(vi)); % set floor value of f at each point, to prevent -Inf
        f(f(:,2,vi) == -Inf,2,vi) = log(.00001/En(vi)); % set floor value of f at each point, to prevent -Inf
        plot(x,ON(:,1,vi).*f(:,1,vi),'g-');
        plot(x,ON(:,2,vi).*f(:,2,vi),'r-');
    end
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 3 10];
    saveas(h,fullfile(plot_dir,sprintf('QMLE_Plot_%s.eps',name)),'epsc2');
    %% Sumed Quantile loglikelihodd over coherence
    h = figure;
    for vi = 1:length(acc)
        Obar(vi) = sum([ON(:,1,vi).*log(OP(:,1,vi)); ON(:,2,vi).*log(OP(:,2,vi))],'omitnan');
        Ebar(vi) = sum([ON(:,1,vi).*(f(:,1,vi)); ON(:,2,vi).*(f(:,2,vi))],'omitnan');
    end
    bar([Obar; Ebar]','grouped');
    saveas(h,fullfile(plot_dir,sprintf('SumLL_Plot_%s.eps',name)),'epsc2');
    %% proportion at each quantile
    h = figure;
    for vi = 1:length(acc)
        subplot(6,1,vi); hold on;
        x = 1:10;
        plot(x, OP(:,1,vi),'gx');
        plot(x, -OP(:,2,vi),'rx');
        plot(x, EN(:,1,vi)/En(vi),'g-');
        plot(x, -EN(:,2,vi)/En(vi),'r-');
        ylim([-.1, .2]);
    end
    saveas(h,fullfile(plot_dir,sprintf('Proportion_Plot_%s.eps',name)),'epsc2');
    
    
    
    %% simulation, plot time course
    if ~exist(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat,sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = WW06Dynamic_FitBhvr8ParamsVI_QMLE_GPU(params, dataDynmc, dataBhvr);
        % [nLL, Chi2, BIC, AIC, rtmat, choicemat,sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDMDynamic_FitBhvr7ParamsIV_QMLE_GPU(params, dataDynmc, dataBhvr);
        %sm_mr1c = gather(sm_mr1c);
        save(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),...
            'rtmat','choicemat','sm_mr1c','sm_mr2c','sm_mr1cD','sm_mr2cD','params');
        toc
    else
        load(fullfile(plot_dir, sprintf('PlotDynamic_%s.mat',name)));
    end
    
    %% plot fitted dynamics and abcd values
    [h, h2, dot_tick_x] = PlotFitDynmcs(RoitmanDataDir, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD, plot_dir, name, 1);
    dot_tick_x
    
    %% Single trials dynamics
    lwd = 1;
    mksz = 3;
    fontsize = 11;
    %    JNp, JNn, I0, noise, miu0, tauS, nLL
    simname = sprintf('WW06Dynmc_JNp%2.1f_JNn%1.2f_I0%1.2f_noise%1.2f_miu0%2.2f_gamma%1.3f_tauS%0.2f_nLL%5.2f',params);
    Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
    c1 = (1 + Cohr)';
    c2 = (1 - Cohr)';
    cplist = [c1, c2];
    mygray = flip(gray(length(cplist)));
    JNp = params(1);
    JNn = params(2);
    I0 = params(3);
    sgm = params(4);
    miu0 = params(5);
    ndt = .03 + .09; % sec, initial dip for 90ms after stimuli onset, resort to the saccade side, the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
    presentt = 0;
    triggert = 0;
    dur = 5;
    dt =.001;
    thresh = 15; % to match the data in Roitman&Shadlen and most of the other evidence
    stimdur = dur;
    stoprule = 1;
    JN = [JNp -JNn; -JNn JNp];
    gamma = .641;
    tauS = params(6);   %.1; % sec
    tauAMPA = .002; % sec
    unit = 1; % secs
    H0 = 32/70*thresh;
    S0 = H0*gamma*tauS/(H0*gamma*tauS+1);
    initialvals = [H0, H0;S0, S0]; % S = H*gamma*tauS./(H*gamma*tauS+1)
    h = figure; hold on;
    filename = sprintf('%s',simname);
    rng(randseed);
    for vi = 2:6
        cp = cplist(vi,:);
        [~, ~, R, ~, ~, ~] = wong06(cp,miu0,sgm,I0,JN,...
            gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
        lgd2(vi-1) = plot(R(:,2), 'k-.', 'Color', mygray(vi,:), 'LineWidth',lwd);
        lgd1(vi-1) = plot(R(:,1), 'k-', 'Color', mygray(vi,:), 'LineWidth',lwd);
    end
    plot([.2, 1.2]/dt,[thresh,thresh], 'k-');
    text(600,thresh*1.1,'threshold');
    yticks([0,32,70]);
    yticklabels({'0','32','70'});
    ylabel('Activity (Hz)');
    % ylim([0,74]);
    xticks([0, 500, 1000, 1500]);
    xticklabels({'0','.5','1.0','1.5'});
    xlim([-50, 1200]);
    xlabel('Time (s)');
    % lgd = legend(lgd3,cellstr(num2str(c3)),...
    %     'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
    %     'FontAngle','italic','NumColumns',1,'Box','off');
    % title(lgd, 'V_3');
    savefigs(h, filename, plot_dir, fontsize, [2, 1.5]);
    
end
