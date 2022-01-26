function [nLL, Chi2, BIC, AIC, rtmat, choicemat] = WW06FitBhvr7ParamsII_QMLE_GPU(params,dataDynmc, dataBhvr)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax';
sac_ax = dataDynmc.sac_ax';
q = dataBhvr.q;
On = dataBhvr.On;
ON = dataBhvr.ON;
OP = dataBhvr.OP;
% parameters to fit
JNp = params(1);
JNn = params(2);
I0 = params(3);
sgm = params(4);
miu0 = params(5);
ndt = .03 + .09; % sec, initial dip for 90ms after stimuli onset, resort to the saccade side, the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
tgap = 0; % tgap is not assumed for the simulation of WW06 

% other fixed parameters
% sims = 1024;
deduction = .1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
presentt = 0;
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; % to match the data in Roitman&Shadlen and most of the other evidence
stimdur = dur;
stoprule = 1;
JN = [JNp -JNn; -JNn JNp];
gamma = .641;
tauS = params(6);   %.1; % sec
tauAMPA = params(7); %.002; % sec
unit = 1; % secs
H0 = 32;
S0 = H0*gamma*tauS/(H0*gamma*tauS+1);
initialvals = [H0, H0;S0, S0]; % S = H*gamma*tauS./(H*gamma*tauS+1)

% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
c1 = (1 + Cohr)';
c2 = (1 - Cohr)';
cp = [c1, c2];
% tic;
[choicemat, rtmat] = wong06_GPU(cp,miu0,sgm,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims);
rtmat = squeeze(rtmat)' + ndt;
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