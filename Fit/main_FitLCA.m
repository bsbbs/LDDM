
%% setup directory
addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);

%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../RecurrentModel/bads/bads-master'));
addpath('../CoreFunctions/');
addpath('./SvrCode/');
out_dir = '../../RecurrentModel/Fit/Rslts/LCAFitBhvr5Params_QMLE_SvrGPU';
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
%     k,    beta, noise, T0, thresh
LB = [-3    0   0      0    .01];
UB = [3   10	5      2      10];
PLB = [-1  .1	.1      0      .1];
PUB = [1   3	.5      .2      6];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

% likelihood function
% parpool(6);
nLLfun = @(params) LCAFitBhvr5Params_QMLE_GPU(params, dataBhvr)
[fvalbest,~,~] = nLLfun(x0)
fprintf('test succeeded\n');

% loop begin
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


%% visulization
if 0
    addpath(genpath('/Users/bs3667/Documents/LDDM/utils'));
    addpath(genpath('/Users/bs3667/Documents/LDDM/CoreFunctions'));
    lwd = 1.5;
    fontsize = 14;
    aspect = [4, 3];
    outdir = '/Volumes/GoogleDrive/My Drive/LDDM/Fit/Rslts/FitLCA';
    %% The dynamic of single trials
    mygray = gray(8);
    h = figure; hold on;
    filename = 'LCA_Dynamic_sngltrials';
    i = 0;
    for c = [0, .032, .064, .128, .356, .512]
        i = i + 1;
        rho = [1+c, 1-c];
        k = .1*eye(2);
        beta = (ones(2) - eye(2))*.6;
        sgm = .2;
        T0 = .1;
        thresh = 1;
        dur = 5;
        dt = .001;
        [choice, rt, x] = LCA(rho, k, beta, sgm, T0, thresh, dur, dt);
        t = (1:length(x(:,1)))*dt;
        plot(t, x(:,1), '-', 'Color' , mygray(i,:), 'LineWidth', lwd);
        plot(t, x(:,2), '--', 'Color' , mygray(i,:), 'LineWidth', lwd);
        
    end
    ylabel('Neural Activity (a.u.)');
    xlabel('Time (Secs)');
    savefigs(h, filename, outdir, fontsize, aspect);
end
