function [nLL, Chi2, BIC, AIC, rtmat, choicemat, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = LCADynmcs_FitBhvr5Params_QMLE_GPU(params, dataDynmc, dataBhvr, sims)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax;
sac_ax = dataDynmc.sac_ax;
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
% parameters to fit
k = params(1)*eye(2);
beta = params(2)*eye(2);
sgm = params(3);
thresh = params(4);
T0 = 30 + 90; 
% other fixed parameters
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
dur = 5;
dt =.001;

V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
RhoMat = [V1, V2];

% tic
[rtmat, choicemat, ~, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = LCA_Dynmc_Trim_GPU(RhoMat, k, beta, sgm, thresh, dur, dt, sims, T0, dot_ax, sac_ax);
rtmat = squeeze(rtmat)';
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