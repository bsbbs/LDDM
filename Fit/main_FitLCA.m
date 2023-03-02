
%% setup directory
addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../RecurrentModel/bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = '../../RecurrentModel/Fit/Rslts/LCAFitBhvr4Params_QMLE_SvrGPU';
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
%     k,    beta, noise, thresh
LB = [-3    0   0      .01];
UB = [3   30	5      10];
PLB = [-1  1	.1      .1];
PUB = [1   10	.5      6];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
nLLfun = @(params) LCAFitBhvr4Params_QMLE_GPU(params, dataBhvr, 10240)
[fvalbest,~,~] = nLLfun(x0)
fprintf('test succeeded\n');

% loop begin
Collect = [];
parfor i = 1:myCluster.NumWorkers*8*4
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


%% visulization
if 0
    % set directory
%     Homedir = 'C:\Users\Bo';
    Homedir = '~';
    addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
    addpath(fullfile(Homedir,'Documents','LDDM','utils'));
    addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
%     Glgdr = 'G:\My Drive\LDDM';
    Glgdr = '/Volumes/GoogleDrive/My Drive/LDDM';
%     Glgdr = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM';
    cd(fullfile(Glgdr, 'Fit'));
    RoitmanDataDir = fullfile(Homedir,'Documents','LDDM', 'Fit', 'RoitmanDataCode');
    dataDynmc = load(fullfile(RoitmanDataDir,'DynmcsData.mat'));
    dataBhvr = LoadRoitmanData(RoitmanDataDir);
    
    out_dir = fullfile(Glgdr,'Fit','Rslts','FitLCA');
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    plot_dir = fullfile(out_dir,'graphics');
    if ~exist(plot_dir,'dir')
        mkdir(plot_dir);
    end
    plot_dir = '~/Desktop';
    params = [0.477431	16.579511	0.353479	2.475246	16947.87064];
    name = sprintf('k%.3f_b%1.2f_sgm%.3f_thresh%1.2f_nLL%5.2f',params);
    %% Simulation given parameters
    sims = 10240;
    if ~exist(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat] = LCAFitBhvr4Params_QMLE_GPU(params, dataBhvr, sims);
        save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
            'rtmat','choicemat','params','nLL','Chi2','AIC','BIC');
        toc
    else
        load(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)));
    end
    
    %% ditribution of RT and fitted line
    h = RT_Dstrbtn(dataBhvr, rtmat, choicemat, plot_dir, name);
    
    %% aggregated RT & ACC
    h = PlotmeanRTACC(RoitmanDataDir, rtmat, choicemat, plot_dir, name);
    
    %% QP plot for reaction time and choice
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
    %% Simulate dynamics
    sims = 10240;
    destfile = fullfile(plot_dir,sprintf('PlotDynmcs_%s.mat',name));
    if ~exist(destfile,'file')
        tic;
        [nLL, Chi2, BIC, AIC, rtmat, choicemat, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD]...
            = LCADynmcs_FitBhvr4Params_QMLE_GPU(params, dataDynmc, dataBhvr, sims);
        save(destfile,...
            'rtmat','choicemat','params','nLL','Chi2','AIC','BIC', 'sm_mr1c', 'sm_mr2c', 'sm_mr1cD', 'sm_mr2cD');
        toc
    else
        load(destfile);
    end
    %% plot fitted dynamics and abcd values
    [h, h2, dot_tick] = PlotFitDynmcs(RoitmanDataDir, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD, plot_dir, name, 1);
    
    %% The dynamic of single trials
    mygray = gray(8);
    h = figure; hold on;
    filename = sprintf('LCA_Dynamic_%s',name);
    i = 0;
    for c = [0, .032, .064, .128, .356, .512]
        i = i + 1;
        rho = [1+c, 1-c];
        k = params(1)*eye(2);
        beta = (ones(2) - eye(2))*params(2);
        sgm = params(3);
        T0 = params(4);
        thresh = params(5);
        dur = 5;
        dt = .001;
        [choice, rt, x] = LCA(rho, k, beta, sgm, T0, thresh, dur, dt);
        t = (1:length(x(:,1)))*dt;
        plot(t, x(:,1), '-', 'Color' , mygray(i,:), 'LineWidth', lwd);
        plot(t, x(:,2), '--', 'Color' , mygray(i,:), 'LineWidth', lwd);
        
    end
    ylabel('Neural Activity (a.u.)');
    xlabel('Time (Secs)');
    savefigs(h, filename, plot_dir, fontsize, aspect);
end
