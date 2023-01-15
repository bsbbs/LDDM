function [Chi2, N, nLL, BIC, AIC, rtmat, choicemat, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LDDM_FitDynmc7Params_OLS_GPU(params, dataDynmc, dataBhvr, sims)
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
OP = dataBhvr.OP;

% parameters to fit
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
tauR = params(4);
tauG = params(5);
tauD = params(6);
Tau = [tauR tauG tauD];
thresh = params(7); % mean(mean(m_mr1cD(sac_ax == -20 | sac_ax == -40,:))); % = 68.6391 % mean(max(m_mr1cD))+1 = 70.8399; 

% other fixed parameters
dot_gap = .19;
sac_gap = .03;
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
    sgm, Tau, dur, dt, thresh, initialvals, stoprule, sims, dot_ax, sac_ax, dot_gap, sac_gap);
rtmat = squeeze(rtmat)';
choicemat = squeeze(choicemat)';
% toc;

%% Weighted least square fitting for neural dynamics
Chi2 = 0;
N = 0;
for i = 1:4
    switch i
        case 1
            fit = sm_mr1c(dot_ax>=round(dot_gap/dt),:);
            dat = m_mr1c(dot_ax>=round(dot_gap/dt),:);
        case 2
            fit = sm_mr2c(dot_ax>=round(dot_gap/dt),:);
            dat = m_mr2c(dot_ax>=round(dot_gap/dt),:);
        case 3
            fit = sm_mr1cD(sac_ax<=-round(sac_gap/dt),:);
            dat = m_mr1cD(sac_ax<=-round(sac_gap/dt),:);
        case 4
            fit = sm_mr2cD(sac_ax<=-round(sac_gap/dt),:);
            dat = m_mr2cD(sac_ax<=-round(sac_gap/dt),:);
    end
    fit(isnan(fit)) = 0; % replace the missing values as zero, which will drive a large deviance from the data values
    sse = (fit - dat).^2;
    Nsse = sum(~isnan(sse(:)));
    Chi2 = Chi2 + sum(sse(:), 'omitnan');
    N = N + Nsse;
end


%% for reaction time by histogram
% QMLE, quantile maximum likelihood estimation
% reference: Heathcote & Australia, and Mewhort, 2002.
% Chi-square, reference: Ratcliff & McKoon, 2007.
LL = log(0);
%Chi2 = 0;
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
%Chi2vec = (EN - ON_adj).^2./EN;
%Chi2 = sum(Chi2vec(:),'omitnan');
n = sum(~isnan(ON(:).*f(:)));
k = numel(params);
BIC = k*log(n) - 2*LL;
AIC = 2*k - 2*LL;
nLL = -LL;