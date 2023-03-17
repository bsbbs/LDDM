function [nLL, Chi2, BIC, AIC, rtmat, choicemat,sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDMDynamic_FitBhvr_QMLE_GPUHand(params, dataDynmc, dataBhvr)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax';
sac_ax = dataDynmc.sac_ax';
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
OP = dataBhvr.OP;
% parameters to fit
w = [1 1; 1 1];
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
B0 = params(4);
tauR = params(5);
tauG = params(6);
tauI = params(7);
Rstar = params(8);
ndt = .09 + .03; % sec, 90ms after stimuli onset, resort to the saccade side,
% the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
presentt = 0; % changed for this version to move the fitting begin after the time point of recovery

% other fixed parameters
% sims = 1024;
deduction = 1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
predur = 0;
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; %70.8399; % mean(max(m_mr1cD))+1; 
stimdur = dur;
stoprule = 1;

% Rstar = 32; % ~ 32 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
I0 = params(2)*Rstar;
initialvals = [Rstar,Rstar; (sum(w(1,:)) - params(2))*Rstar,(sum(w(2,:)) - params(2))*Rstar; I0, I0];
eqlb = Rstar; % set equilibrium value before task as R^*
scale = max([5, (((2*mean(w,'all') - params(2)))*eqlb.^2 + (1-a(1)).*eqlb)/(1 + B0)]);
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = ([V1, V2]+B0)*scale;
Vprior = (ones(size(Vinput)) + B0)*scale;

Tau = [tauR tauG tauI];
% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
% tic;
[rtmat, choicemat, ~, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDM_Dynmc_Trim_GPU(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
rtmat = squeeze(rtmat)'+ndt;
choicemat = squeeze(choicemat)';
% toc;

%% for reaction time by histogram
% QMLE, quantile maximum likelihood estimation
% reference: Heathcote & Australia, and Mewhort, 2002.
% Chi-square, reference: Ratcliff & McKoon, 2007.
LL = log(0);
Chi2 = 0;
h = figure;
for vi = 1:6
    En(vi) = numel(rtmat(:,vi));
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    if ~all(isnan(q(:,1,vi)))
        tmp = histogram(RT_corr, [0; q(:,1,vi); Inf], 'Visible',0);
        EN(:,1,vi) = tmp.Values;
    else
        EN(:,1,vi) = NaN(numel(q(:,1,vi))+1,1);
    end
    if ~all(isnan(q(:,2,vi)))
        tmp = histogram(RT_wro, [0; q(:,2,vi); Inf], 'Visible',0);
        EN(:,2,vi) = tmp.Values;
    else
        EN(:,2,vi) =  NaN(numel(q(:,2,vi))+1,1);
    end
    f(:,:,vi) = log((EN(:,:,vi)/En(vi)));
    f(f(:,1,vi) == -Inf,1,vi) = log(.1/En(vi)); % set floor value of f at each point, to prevent -Inf
    f(f(:,2,vi) == -Inf,2,vi) = log(.1/En(vi)); % set floor value of f at each point, to prevent -Inf
    ON_adj(:,:,vi) = ON(:,:,vi)*En(vi)./On(vi);
end
close(h);
LL = sum(ON(:).*f(:),'omitnan');
Chi2vec = (EN - ON_adj).^2./EN;
Chi2 = sum(Chi2vec(:),'omitnan');
n = sum(~isnan(ON(:).*f(:)));
k = numel(params);
BIC = k*log(n) - 2*LL;
AIC = 2*k - 2*LL;
nLL = -LL;