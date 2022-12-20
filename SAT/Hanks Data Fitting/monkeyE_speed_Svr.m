%% Main file to fit Hanks data
% write your own help language here

%% setup directory
addpath('../../utils');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);


%% Model fitting with Bayesian Adaptive Direct Search (BADS) optimization algorithm
addpath(genpath('../../../RecurrentModel/bads/bads-master'));
out_dir = '../../../LDDM_Output/SAT/Hanks/monkeyE_speed_Svr';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% Take data from Hanks, et al., 2014
% monkey E
dataBhvr = load_data("behavData_eli.mat");

%% define the range of the parameters
% Define optimization starting point and bounds
%     a,    b, noise, scale, Tau
LB = [0    0.1   .1    .1*256 [.001,.001,.001]];
UB = [70   3	100  20*256 [1,1,1]];
PLB = [15  .9	5    1*256 [.01 .01 .01]];
PUB = [60   1.7 40  8*256   [.2 .2 .2]];

% Randomize initial starting point inside plausible box
x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;

%% define the negative loglikelihood function (nLLfun)
nLLfun = @(params) LDDM_fit_speed(params, dataBhvr);

%% first attempt to evaluate the nLLfun
[fvalbest, ~, ~] = nLLfun(x0);
fprintf('test succeeded\n');

%% start to fit
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

    Collect(i).rndseed = t;
    Collect(i).x0 = x0;
    Collect(i).xest = xest;
    Collect(i).fval = fval;
    Collect(i).output = output;
end
t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7 + i*10^5;
save(fullfile(out_dir,sprintf('CollectRslts%i.mat',t)),'Collect');
