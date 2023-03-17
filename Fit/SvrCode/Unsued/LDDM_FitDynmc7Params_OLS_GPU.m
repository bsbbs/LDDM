function [nll1, BIC, AIC, rtmat, choicemat, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDM_FitDynmc7Params_OLS_GPU(params, dataDynmc, dataBhvr, sims)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax;
sac_ax = dataDynmc.sac_ax;
m_mr1c = dataDynmc.m_mr1c;
m_mr2c = dataDynmc.m_mr2c;
m_mr1cD = dataDynmc.m_mr1cD;
m_mr2cD = dataDynmc.m_mr2cD;
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
%% cut the dynamic data to save computatioanl power
dot_gap = 190; % ms
sac_gap = 30; % ms
Lcut = find(dot_ax <= dot_gap, 1, 'last' )+1; % start right from the time point to fit
Rcut = find(isnan(m_mr1c(:,1)), 1 )-1; % exclude time points where m_mr*c are nans
dot_axcut = dot_ax(Lcut:Rcut); 
m_mr1ccut = m_mr1c(Lcut:Rcut,:);
m_mr2ccut = m_mr2c(Lcut:Rcut,:);
LcutD = find(isnan(m_mr1cD(:,1)), 1, 'last' )+1; % exclude time points where m_mr*c are nans
RcutD = find(sac_ax <= -sac_gap, 1, 'last' ); % stop right after the time point to fit
sac_axcut = sac_ax(LcutD:RcutD);
m_mr1cDcut = m_mr1cD(LcutD:RcutD,:);
m_mr2cDcut = m_mr2cD(LcutD:RcutD,:);
%% parameters to fit
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
tauR = params(4);
tauG = params(5);
tauD = params(6);
Tau = [tauR tauG tauD];
thresh = params(7); % mean(mean(m_mr1cD(sac_ax == -20 | sac_ax == -40,:))); % = 68.6391 % mean(max(m_mr1cD))+1 = 70.8399; 

% other fixed parameters
% ndt = .19 + .03; % sec, 190ms after stimuli onset, resort to the saccade side,
% the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
if nargin < 4
    sims = 10240;
end
w = ones(2);
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
dur = 5;
dt =.001;
stoprule = 1;
Rstar = mean(mean(m_mr1c(dot_ax == 180 | dot_ax == 200,:))); % = 43.1026 % ~ 32 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
scale = (1-params(1))*Rstar + (sum(w(1,:)) - params(2))*Rstar^2;
initialvals = [Rstar,Rstar; (sum(w(1,:)) - b(1,1))*Rstar,(sum(w(2,:)) - b(2,2))*Rstar; b(1,1)*Rstar,b(2,2)*Rstar];
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;

% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
% tic;
[rtmat, choicemat, ~, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDM_Dynmc_gap_GPU(Vinput, w, a, b,...
    sgm, Tau, dur, dt, thresh, initialvals, stoprule, sims, dot_axcut, sac_axcut, dot_gap, sac_gap);
rtmat = squeeze(rtmat)';
choicemat = squeeze(choicemat)';
% toc;

%% least squared error fitting for neural dynamics
[ll1, ~, N1] = DynamicLLFun(m_mr1ccut,m_mr2ccut,m_mr1cDcut,m_mr2cDcut,sm_mr1c,sm_mr2c,sm_mr1cD,sm_mr2cD, dot_axcut, dot_gap, sac_axcut, sac_gap);
nll1 = -ll1;
%% for reaction time and choice
[ll2, ~, N2] = BehaviorLLFun(q, On ,ON, rtmat, choicemat);

%% overall
k = numel(params);
LL = ll1 + ll2;
BIC = k*log(N1+N2) - 2*LL;
AIC = 2*k - 2*LL;
nLL = -LL;

%% fill the empty values on the dynamics
sm_mr1c = [nan(Lcut-1,6); sm_mr1c; nan(numel(dot_ax)-Rcut,6)];
sm_mr2c = [nan(Lcut-1,6); sm_mr2c; nan(numel(dot_ax)-Rcut,6)];
sm_mr1cD = [nan(LcutD-1,6); sm_mr1cD; nan(numel(sac_ax)-RcutD,6)];
sm_mr2cD = [nan(LcutD-1,6); sm_mr2cD; nan(numel(sac_ax)-RcutD,6)];