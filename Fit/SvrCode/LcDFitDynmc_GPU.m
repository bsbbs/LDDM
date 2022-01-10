function nLL = LcDFitDynmc_GPU(params,data)
% reload Roitman's data, processed
dot_ax = data.dot_ax';
sac_ax = data.sac_ax';
m_mr1c = data.m_mr1c';
m_mr2c = data.m_mr2c';
m_mr1cD = data.m_mr1cD';
m_mr2cD = data.m_mr2cD';
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
% fprintf('GPU Simulations %i chains ...\t', sims);
V1 = 256*(1 + Cohr)';
V2 = 256*(1 - Cohr)';
Vinput = [V1, V2]*scale;
% tic;
[~, ~, ~, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD] = LcD_Dynmc_Trim_GPU(Vprior, Vinput, w, a*eye(2), b*eye(2),...
    sgm, Tau, predur, dur, dt, ndt, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax);
% toc;

% calculate Chi-Square, by histogram
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
Chi2 = gather(Chi2)/1000; % rescale the value a bit
ll = -.5*Chi2/sigma^2 - .5*N*log(2*pi*sigma^2);
nLL = -ll;
