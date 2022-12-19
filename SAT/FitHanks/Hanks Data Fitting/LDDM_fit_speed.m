function [nLL, Chi2, BIC, AIC, simulate_speed_rt, simulate_speed_choice] = LDDM_fit_speed(params, dataBhvr)
% reload Hanks data, processed
speed_coh = dataBhvr.speed_coh;
empirical_speed_qmat = dataBhvr.speed_qmat;
speed_empirical_num = dataBhvr.speed_observe_num;
speed_empirical_total = dataBhvr.speed_observe_total;

% parameters to fit
a = params(1)*eye(2);
b = params(2)*eye(2);
sgm = params(3);
tauR = params(5);
tauG = params(6);
tauI = params(7);
ndt = .09 + .03; % sec, 90ms after stimuli onset, resort to the saccade side,
% the activities reaches peak 30ms before initiation of saccade, according to Roitman & Shadlen
presentt = 0; % changed for this version to move the fitting begin after the time point of recovery
scale = params(4);

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
Rstar = 32; % ~ 32 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar; 0,0];
eqlb = Rstar; % set equilibrium value before task as R^*

Tau = [tauR tauG tauI];
% simulation
% fprintf('GPU Simulations %i chains ...\t', sims);
V1 = (1 + Cohr)';
V2 = (1 - Cohr)';
Vinput = [V1, V2]*scale;
Vprior = ones(size(Vinput))*(2*mean(w,'all')*eqlb.^2 + (1-a(1)).*eqlb);
% tic
[simulate_speed_rt, simulate_speed_choice, ~] = LDDM_GPU(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
simulate_speed_rt = squeeze(simulate_speed_rt)'+ndt;
simulate_speed_choice = squeeze(simulate_speed_choice)';
% toc;

%% for reaction time by histogram
% QMLE, quantile maximum likelihood estimation
% reference: Heathcote & Australia, and Mewhort, 2002.
% Chi-square, reference: Ratcliff & McKoon, 2007.

% speed
speed_h = figure;
speed_simulate_total = [];
speed_estimate_qnum = [];
speed_f = [];

for vi = 1: length(speed_coh)
    correct_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==1, vi);
    %wrong_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==0, vi);
    wrong_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==2, vi);
    speed_simulate_total(vi) = length(simulate_speed_rt(:, vi));

    if ~any(isnan(empirical_speed_qmat(:, 1, vi)))
        temp = histogram(correct_RT, [0; empirical_speed_qmat(:, 1, vi); Inf]);
        speed_estimate_qnum(:, 1, vi) = temp.Values;
    else
        speed_estimate_qnum(:, 1, vi) = NaN(length(empirical_speed_qmat(:, 1, vi))+1, 1);
    end

    if ~any(isnan(empirical_speed_qmat(:, 2, vi)))
        temp = histogram(wrong_RT, [0; empirical_speed_qmat(:, 2, vi); Inf]);
        speed_estimate_qnum(:, 2, vi) = temp.Values;
    else
        speed_estimate_qnum(:, 2, vi) = NaN(length(empirical_speed_qmat(:, 2, vi))+1, 1);
    end

    speed_f(:, :, vi) = log(speed_estimate_qnum(:, :, vi) / speed_simulate_total(vi));
    speed_f(speed_f(:, 1, vi)==-Inf, 1, vi) = log(.1/speed_simulate_total(vi));
    speed_f(speed_f(:, 2, vi)==-Inf, 2, vi) = log(.1/speed_simulate_total(vi));
    speed_empirical_adj = speed_empirical_num(:, :, vi) * speed_simulate_total(vi) ./ speed_empirical_total(vi);
end
close(speed_h);

speed_LL = sum(speed_empirical_num(:) .* speed_f(:), 'omitnan');
nLL = -speed_LL;

Chi2vec = (speed_estimate_qnum - speed_empirical_adj) .^ 2 ./ speed_estimate_qnum;
Chi2 = sum(Chi2vec(:),'omitnan');

n = sum(~isnan(speed_empirical_num(:) .* speed_f(:)));
k = numel(params);
BIC = k*log(n) - 2*speed_LL;
AIC = 2*k - 2*speed_LL;