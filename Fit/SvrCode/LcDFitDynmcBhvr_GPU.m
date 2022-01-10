function nLL = LcDFitDynmcBhvr_GPU(params,dataDynmc, dataBhvr)
% reload Roitman's data, processed
dot_ax = dataDynmc.dot_ax';
sac_ax = dataDynmc.sac_ax';
m_mr1c = dataDynmc.m_mr1c';
m_mr2c = dataDynmc.m_mr2c';
m_mr1cD = dataDynmc.m_mr1cD';
m_mr2cD = dataDynmc.m_mr2cD';
numbins = dataBhvr.numbins;
histmat = dataBhvr.histmat;
rtrange = dataBhvr.rtrange;
% parameters to fit
a = params(1);
b = params(2);
sgm = params(3);
tauR = params(4);
tauG = params(5);
tauI = params(6);
ndt = params(7);
scale = params(8);
% std of RSS for ll computing purpose
sigma = params(9);

% other fixed parameters
% sims = 1024;
deduction = 1;
sims = 1024/deduction;
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
predur = 2;
triggert = 0;
dur = 5;
dt =.001;
thresh = mean(max(m_mr1cD'))+1; % 70.8399
stimdur = dur;
stoprule = 1;
Rstar = mean(mean([m_mr1c(1:5,:); m_mr2c(1:5,:)])); % 34.6004
initialvals = [Rstar,Rstar; 2*Rstar,2*Rstar; 0,0];
Vprior = ones(6,2)*((1-a(1,1))*Rstar + 2*Rstar^2);
w = [1 1; 1 1];
Tau = [tauR tauG tauI];
% simulation
fprintf('GPU Simulations %i chains ...\t', sims);
V1 = 256*(1 + Cohr)';
V2 = 256*(1 - Cohr)';
Vinput = [V1, V2]*scale;
% tic;
[rtmat, choicemat, ~, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LcD_Dynmc_Trim_GPU(Vprior, Vinput, w, a*eye(2), b*eye(2),...
    sgm, Tau, predur, dur, dt, ndt, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
rtmat = squeeze(rtmat);
choicemat = squeeze(choicemat);
% toc;

% calculate Chi-Square, 1st for dynamic 
Chi2 = 0;
N = 0;
for i = 1:4
    switch i
        case 1
            fit = sm_mr1c(dot_ax>=0,:);
            dat = m_mr1c(dot_ax>=0,:);
        case 2
            fit = sm_mr2c(dot_ax>=0,:);
            dat = m_mr2c(dot_ax>=0,:);
        case 3
            fit = sm_mr1cD(sac_ax<=0,:);
            dat = m_mr1cD(sac_ax<=0,:);
        case 4
            fit = sm_mr2cD(sac_ax<=0,:);
            dat = m_mr2cD(sac_ax<=0,:);
    end
    sse = (fit - dat).^2;
    Ndat = sum(~isnan(dat(:)));
    Nsse = sum(~isnan(sse(:)));
    if Nsse > 10 
        Chi2 = Chi2 + sum(sse(~isnan(sse)))/Nsse*Ndat; % replace missing value by mean
    else % special case, simulation predicts weird no valid data
        Chi2 = Chi2 + sum(dat(~isnan(dat)).^2);
    end
    N = N + Ndat;
end
Chi2 = gather(Chi2);
% ll = -.5*Chi2/sigma^2 - .5*N*log(2*pi*sigma^2);
% nLL = -ll;

% 2nd for reaction time by histogram
for vi = 1:6
    simhistvec = [];
    histmatr = histmat(vi,:);
    histbinsr = linspace(rtrange(vi,1),rtrange(vi,2),numbins);
    histbins = histbinsr;
    rtmin = min(rtmat(vi,:));
    rtmax = max(rtmat(vi,:));
    if rtmin < rtrange(vi,1)
        histbins = [rtmin,histbinsr];
        histmatr = [0 histmat(vi,1:(numbins-1)) 0 histmat(vi,numbins:2*(numbins-1))];
    end
    if rtmax > rtrange(vi,2)
        histbins = [histbins,rtmax];
        n = length(histmatr);
        histmatr = [histmatr(1:n/2), 0, histmatr(n/2+1:n), 0];
    end
    
    SlctTrial = choicemat(vi,:) == 1;
    histg_corr = histogram(rtmat(vi,SlctTrial),'BinEdges',histbins, 'Visible',0);
    tmp1 = histg_corr.Values/sims;
    SlctTrial = choicemat(vi,:) == 2;
    histg_wro = histogram(rtmat(vi,SlctTrial),'BinEdges',histbins, 'Visible',0);
    tmp2 = histg_wro.Values/sims;
    simhistvec = [tmp1 tmp2];
    Chi2 = Chi2 + sum((simhistvec*1024 - histmatr*1024).^2);
    N = N + numel(simhistvec);
end
ll = -.5*Chi2/sigma^2 - .5*N*log(2*pi*sigma^2);
nLL = -ll;