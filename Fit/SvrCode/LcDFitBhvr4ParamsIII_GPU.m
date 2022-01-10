function [nLL, R_squared, N] = LcDFitBhvr4ParamsIII_GPU(params,dataDynmc, dataBhvr)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax';
sac_ax = dataDynmc.sac_ax';
m_mr1c = dataDynmc.m_mr1c';
m_mr2c = dataDynmc.m_mr2c';
m_mr1cD = dataDynmc.m_mr1cD';
m_mr2cD = dataDynmc.m_mr2cD';
numbins = dataBhvr.numbins;
histmat = dataBhvr.histmat;
histbin = dataBhvr.histbins;
rtrange = dataBhvr.rtrange;
% parameters to fit
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
tauR = .1; %params(4);
tauG = .1; %params(5);
tauI = .1; %params(6);
ndt = .19 + .03; % sec, 90ms after stimuli onset, activities show ~100 ms of dip and recovery;
% resort to the saccade side, the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
tgap = 0; % changed for this version to move the fitting begin after the time point of recovery
scale = params(4);
% std of RSS for ll computing purpose
sigma = params(5);

% other fixed parameters
% sims = 1024;
deduction = .1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
predur = 0;
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; %70.8399; % mean(max(m_mr1cD))+1; 
stimdur = dur;
stoprule = 1;
w = [1 1; 1 1];
Rstar = 42; % ~ 42 Hz according to Roitman and Shadlen's data
initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar; 0,0];
Vprior = ones(6,2)*((1-a(1,1))*Rstar + 2*Rstar^2);

Tau = [tauR tauG tauI];
% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;
% tic;
[rtmat, choicemat, ~, ~, ~, ~, ~] = LcD_Dynmc_Trim_GPU(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, tgap, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
rtmat = squeeze(rtmat)+ndt;
choicemat = squeeze(choicemat);
% toc;

% for reaction time by histogram
SST = 0;
SSE = 0;
Chi2 = 0;
N = 0;
for vi = 1:6
    simhistvec = [];
    histmatr = histmat(vi,:);
    % histbinsr = linspace(rtrange(vi,1),rtrange(vi,2),numbins);
    % histbins = histbinsr;
    histbins = histbin(vi,:);
    rtmin = min(rtmat(vi,:));
    rtmax = max(rtmat(vi,:));
    histbins = [min(rtrange(vi,1),rtmin)*.9, histbins, max(rtrange(vi,2),rtmax)/.9];
    histmatr = [0 histmat(vi,1:(numbins)) 0, 0 histmat(vi,(numbins+1):(2*numbins)) 0];
    SlctTrial = choicemat(vi,:) == 1;
    histg_corr = histogram(rtmat(vi,SlctTrial),'BinEdges',histbins);
    tmp1 = histg_corr.Values/sims;
    SlctTrial = choicemat(vi,:) == 2;
    histg_wro = histogram(rtmat(vi,SlctTrial),'BinEdges',histbins);
    tmp2 = histg_wro.Values/sims;
    simhistvec = [tmp1 tmp2];
    Chi2 = Chi2 + sum((simhistvec - histmatr).^2);
    N = N + numel(simhistvec);
    SSE = SSE + sum((simhistvec - histmatr).^2);
    SST = SST + sum(histmatr.^2);
end
R_squared = 1 - SSE/SST;
ll = -.5*Chi2/sigma^2 - .5*N*log(2*pi*sigma^2);
nLL = -ll;