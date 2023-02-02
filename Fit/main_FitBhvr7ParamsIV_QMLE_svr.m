addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../RecurrentModel/bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = '../../RecurrentModel/Fit/Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU';
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
%     a,    b, noise, scale, Tau
LB = [0    0.1   .1    .1*256 [.001,.001,.001]];
UB = [70   3	100  20*256 [1,1,1]];
PLB = [15  .9	5    1*256 [.01 .01 .01]];
PUB = [60   1.7	40   8*256 [.2 .2 .2]];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
nLLfun = @(params) LDDMFitBhvr7ParamsIV_QMLE_GPU(params, dataBhvr);
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
    % Homedir = 'C:\Users\Bo';
    Homedir = '~';
    addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
    addpath(fullfile(Homedir,'Documents','LDDM','utils'));
    addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
    % Glgdr = 'G:\My Drive\LDDM';
    Glgdr = '/Volumes/GoogleDrive/My Drive/LDDM';
    cd(fullfile(Glgdr, 'Fit'));
    RoitmanDataDir = fullfile(Glgdr, 'Fit', 'RoitmanDataCode');
    dataDynmc = load(fullfile(RoitmanDataDir,'DynmcsData.mat'));
    dataBhvr = LoadRoitmanData(RoitmanDataDir);
    
    out_dir = fullfile(Glgdr, 'Fit', 'Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU');
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    plot_dir = fullfile(out_dir,'graphics');
    if ~exist(plot_dir,'dir')
        mkdir(plot_dir);
    end
    
    randseed = 24356545;
    rng(randseed,'twister');
    % a, b, noise, scale, tauRGI, nLL
    params = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132 16539.138186]; % the old results
    name = sprintf('a%2.2f_b%1.2f_sgm%2.1f_scale%4.1f_tau%1.2f_%1.2f_%1.2f_nLL%5.2f',params);
    if ~exist(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat] = LDDMFitBhvr7ParamsIV_QMLE_GPU(params, dataBhvr);
        toc
        save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
            'rtmat','choicemat','params','nLL','Chi2','AIC','BIC');
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
        % ylim(-[1000,1]);
        
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
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 3 10];
    %saveas(h,fullfile(plot_dir,sprintf('QMLE_Plot_%s.fig',name)),'fig');
    saveas(h,fullfile(plot_dir,sprintf('QMLE_Plot_%s.eps',name)),'epsc2');
    %% Sumed Quantile loglikelihodd over coherence
    h = figure;
    for vi = 1:length(acc)
        Obar(vi) = sum([ON(:,1,vi).*log(OP(:,1,vi)); ON(:,2,vi).*log(OP(:,2,vi))],'omitnan');
        Ebar(vi) = sum([ON(:,1,vi).*(f(:,1,vi)); ON(:,2,vi).*(f(:,2,vi))],'omitnan');
    end
    bar([Obar; Ebar]','grouped');
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 4 5];
    %saveas(h,fullfile(plot_dir,sprintf('SumLL_Plot_%s.fig',name)),'fig');
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
        %ylim([-.1, .2]);
    end
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 3 10];
    %saveas(h,fullfile(plot_dir,sprintf('Proportion_Plot_%s.fig',name)),'fig');
    saveas(h,fullfile(plot_dir,sprintf('Proportion_Plot_%s.eps',name)),'epsc2');
    
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
    
    
    %% Simulate dynamics
    % if ~exist(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),'file')
    %     tic;
    %     [nLL, Chi2, BIC, AIC, rtmat, choicemat,sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDMDynamic_FitBhvr7ParamsIV_QMLE_GPU(params, dataDynmc, dataBhvr);
    %     save(fullfile(plot_dir,sprintf('PlotDynamic_%s.mat',name)),...
    %         'rtmat','choicemat','sm_mr1c','sm_mr2c','sm_mr1cD','sm_mr2cD','params');
    %     toc
    % else
    %     load(fullfile(plot_dir, sprintf('PlotDynamic_%s.mat',name)));
    % end
    destfile = fullfile(plot_dir,sprintf('PlotClassDynamic_%s.mat',name));
    if ~exist(destfile,'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD, ...
            sm_mg1c, sm_mg2c, sm_mg1cD, sm_mg2cD, sm_md1c, sm_md2c, sm_md1cD, sm_md2cD]...
            = LDDMClassDynamic_FitBhvr7ParamsIV_QMLE_GPU(params, dataDynmc, dataBhvr);
        save(destfile,...
            'rtmat','choicemat','sm_mr1c','sm_mr2c','sm_mr1cD','sm_mr2cD',...
            'sm_mg1c', 'sm_mg2c', 'sm_mg1cD', 'sm_mg2cD', 'sm_md1c', 'sm_md2c',...
            'sm_md1cD', 'sm_md2cD', 'params');
        toc
    else
        load(destfile);
    end
    %% plot fitted dynamics and abcd values
    [h, h2, dot_tick] = PlotFitDynmcs(RoitmanDataDir, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD, plot_dir, ['r_',name]);
    [h, h2, dot_tick] = PlotFitDynmcs(RoitmanDataDir, sm_mg1c, sm_mg2c, sm_mg1cD, sm_mg2cD, plot_dir, ['g_',name]);
    [h, h2, dot_tick] = PlotFitDynmcs(RoitmanDataDir, sm_md1c, sm_md2c, sm_md1cD, sm_md2cD, plot_dir, ['d_',name]);
    
    %% Single trials dynamics
    lwd = 1;
    mksz = 3;
    fontsize = 11;
    % a, b, noise, scale, tauRGI, nLL
    simname = sprintf('LDDM_Dynmc_a%2.2f_b%1.2f_sgm%2.1f_scale%4.1f_tau%1.2f_%1.2f_%1.2f_nLL%4.0f',params);
    
    a = params(1)*eye(2);
    b = params(2)*eye(2);
    sgm = params(3)/50;
    tauR = params(5);
    tauG = params(6);
    tauD = params(7);
    Tau = [tauR tauG tauD];
    ndt = .09 + .03; % sec, 90ms after stimuli onset, resort to the saccade side,
    % the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
    presentt = 0; % changed for this version to move the fitting begin after the time point of recovery
    scale = params(4);
    
    predur = 0;
    triggert = 0;
    dur = 5;
    dt =.001;
    thresh = 70; %70.8399; % mean(max(m_mr1cD))+1;
    stimdur = dur;
    stoprule = 1;
    w = [1 1; 1 1];
    Rstar = 32; % ~ 32 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
    initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar; 0,0];
    eqlb = Rstar; % set equilibrium value before task as R^*
    Vprior = [1, 1]*(2*mean(w,'all')*eqlb.^2 + (1-a(1)).*eqlb);
    
    Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
    c1 = (1 + Cohr)';
    c2 = (1 - Cohr)';
    cplist = [c1, c2];
    mygray = flip(gray(length(cplist)));
    
    h = figure;
    
    filename = sprintf('%s',simname);
    randseed = 75245522;
    rng(randseed, 'twister');
    for vi = 2:6
        Vinput = cplist(vi,:)*scale;
        [~, ~, R, G, D] = LDDM(Vprior, Vinput, w, a, b, sgm, Tau, predur, dur,...
            dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        subplot(3,1,1);
        hold on;
        lgd2(vi-1) = plot(R(:,2), 'k--', 'Color', mygray(vi,:), 'LineWidth',lwd);
        lgd1(vi-1) = plot(R(R(:,1)<=thresh,1), 'k-', 'Color', mygray(vi,:), 'LineWidth',lwd);
        subplot(3,1,2);
        hold on;
        lgd2(vi-1) = plot(G(:,2), 'k--', 'Color', mygray(vi,:), 'LineWidth',lwd);
        lgd1(vi-1) = plot(G(R(:,1)<=thresh,1), 'k-', 'Color', mygray(vi,:), 'LineWidth',lwd);
        subplot(3,1,3);
        hold on;
        lgd2(vi-1) = plot(D(:,2), 'k--', 'Color', mygray(vi,:), 'LineWidth',lwd);
        lgd1(vi-1) = plot(D(R(:,1)<=thresh,1), 'k-', 'Color', mygray(vi,:), 'LineWidth',lwd);
    end
    subplot(3,1,1);
    plot([.2, 1.2]/dt,[thresh,thresh], 'k-');
    text(600,thresh*1.1,'threshold');
    yticks([0,32,70]);
    yticklabels({'0','32','70'});
    ylim([0,74]);
    for i = 1:3
        subplot(3,1,i);
        xticks([0, 500, 1000, 1500]);
        xticklabels({'0','.5','1.0','1.5'});
        xlim([-50, 1200]);
        ylabel('Activity (Hz)');
    end
    subplot(3,1,3);
    xlabel('Time (s)');
    for i = 1:3
        subplot(3,1,i);
        savefigs(h, filename, plot_dir, fontsize, [3 6.75]);
    end
    
end