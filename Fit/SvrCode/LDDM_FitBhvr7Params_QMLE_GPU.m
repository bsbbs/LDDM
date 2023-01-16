function [nLL, BIC, AIC, rtmat, choicemat] = LDDM_FitBhvr7Params_QMLE_GPU(params, dataBhvr, sims)
% reload Roitman's data, processed
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
OP = dataBhvr.OP;
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
Rstar = 43.1026; % mean(mean(m_mr1c(dot_ax == 180 | dot_ax == 200,:))), according to Roitman and Shadlen's data
scale = (1-params(1))*Rstar + (sum(w(1,:)) - params(2))*Rstar^2;
initialvals = [Rstar,Rstar; (sum(w(1,:)) - b(1,1))*Rstar,(sum(w(2,:)) - b(2,2))*Rstar; b(1,1)*Rstar,b(2,2)*Rstar];
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;

% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
% tic;
[rtmat, choicemat, ~] = LDDM_gap_GPU(Vinput, w, a, b,...
    sgm, Tau, dur, dt, thresh, initialvals, stoprule, sims, dot_gap, sac_gap);
rtmat = squeeze(rtmat)';
choicemat = squeeze(choicemat)';
% toc;

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