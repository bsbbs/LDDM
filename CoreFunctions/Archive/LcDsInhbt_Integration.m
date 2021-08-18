function [R, G, I, rt, choice, Vtrace] = LcDsInhbt_Integration(miuV, sigmaV, w, alpha, beta, sigma, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule)
%%%%%%%%%%%%%%%%%%
% Vinput, w, alpha, beta, Tau, dt are all the same meaning as
% its own name in the model.
% thresh = threshold
% sigma is the magnitude of noise
% dur is the time for simulation in every single trial, unit in seconds.
% presentt controls the time of stimuli presenting.
% triggert controls the time of switching state, i.e. beta change from 0 to a
% positive number. 
% initialvals defines the initial values of R, G, and I. Defaultly all 
% initial values set as zeros. 
% stimdur defines the duration of stimuli presenting. After withdrawal of 
% stimuli, all input values turn into zeros. 
% stoprule is a control parameter to tell the program if need
% to stop simulating when firing rates of one of the R neurons hit the
% threshold, 1 is to stop, 0 is to continue simulating until using up the 
% rest of time for the trial.
%%%%%%%%%%%%%%%%%%%

tauN =0.002;
timeconstrain = dur - triggert;
sizeVraw = size(miuV);
if exist('initialvals','var')
    R(1,:) = initialvals(1,:);
    G(1,:) = initialvals(2,:);
    I(1,:) = initialvals(3,:);
else
    R(1,:) = zeros(sizeVraw);
    G(1,:) = zeros(sizeVraw);
    I(1,:) = zeros(sizeVraw);
end
if ~exist('stimdur', 'var')
    stimdur = .7;
end
InoiseG=zeros(sizeVraw);InoiseR=zeros(sizeVraw);InoiseI=zeros(sizeVraw);
%% 0st stage, before presenting, V = 0, alpha = 0
for ti = 2:(presentt/dt)
    dI = (-I(ti-1,:) + R(ti-1,:) * 0) * dt/Tau(3);
    dG = (-G(ti-1,:) + R(ti-1,:) * w - I(ti-1,:)) * dt/Tau(2);
    dR = (-R(ti-1,:) + 0./(1+G(ti-1,:))) * dt/Tau(1);
    dInoiseG = (-InoiseG(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseI = (-InoiseI(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseR = (-InoiseR(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    G(ti,:) = G(ti-1,:) + dG + InoiseG(ti-1,:);
    R(ti,:) = R(ti-1,:) + dR + InoiseR(ti-1,:);
    I(ti,:) = I(ti-1,:) + dI + InoiseI(ti-1,:);
    InoiseG(ti,:) = InoiseG(ti-1,:) + dInoiseG;
    InoiseI(ti,:) = InoiseI(ti-1,:) + dInoiseI;
    InoiseR(ti,:) = InoiseR(ti-1,:) + dInoiseR;
    
    G(ti,G(ti,:) < 0) = 0;
    R(ti,R(ti,:) < 0) = 0;
    I(ti,I(ti,:) < 0) = 0;
    
end
if isempty(ti)
    ti = 1;
end
%% 1st stage, after presenting, before decision, V comes in, beta still = 0, working memory on (alpha)
Vtrace = [];
ti1st = ti;
V = miuV + sigmaV * randn(sizeVraw);
for ti = (ti1st + 1):((triggert)/dt)
    if mod(ti*dt,.05) == 0
        V = miuV + sigmaV * randn(sizeVraw);
    end
    if ti > (presentt + stimdur)/dt
        V = zeros(sizeVraw);
    end
    Vtrace = [Vtrace; V];
    dI = (-I(ti-1,:) + R(ti-1,:) * 0) * dt/Tau(3);
    dG = (-G(ti-1,:) + R(ti-1,:) * w - I(ti-1,:)) * dt/Tau(2);
    dR = (-R(ti-1,:) + (V + R(ti - 1,:)*alpha)./(1+G(ti-1,:))) * dt/Tau(1);
    dInoiseG = (-InoiseG(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseI = (-InoiseI(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseR = (-InoiseR(ti-1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    G(ti,:) = G(ti-1,:) + dG + InoiseG(ti-1,:);
    R(ti,:) = R(ti-1,:) + dR + InoiseR(ti-1,:);
    I(ti,:) = I(ti-1,:) + dI + InoiseI(ti-1,:);
    InoiseG(ti,:) = InoiseG(ti-1,:) + dInoiseG;
    InoiseI(ti,:) = InoiseI(ti-1,:) + dInoiseI;
    InoiseR(ti,:) = InoiseR(ti-1,:) + dInoiseR;
    
    G(ti,G(ti,:) < 0) = 0;
    R(ti,R(ti,:) < 0) = 0;
    I(ti,I(ti,:) < 0) = 0;
end
if isempty(ti)
    ti = ti1st;
end

%% 2nd stage, after trigger, decision, beta up
dummy = 0;
rt = NaN;
choice = NaN;
tii = 0;
while (tii < timeconstrain/dt) && dummy == 0
    if mod(tii*dt,.05) == 0
        V = miuV + sigmaV * randn(sizeVraw);
    end
    if ti + tii > (presentt + stimdur)/dt
        V = zeros(sizeVraw);
    end
    Vtrace = [Vtrace; V];
    tii = tii + 1;
    dI = (-I(ti + tii -1,:) + R(ti + tii -1,:) * beta) * dt/Tau(3);
    dG = (-G(ti + tii -1,:) + R(ti + tii -1,:) * w - I(ti + tii -1,:)) * dt/Tau(2);
    dR = (-R(ti + tii -1,:) + (V + R(ti + tii - 1,:)*alpha)./(1+G(ti + tii -1,:))) * dt/Tau(1);
    dInoiseG = (-InoiseG(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseI = (-InoiseI(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    dInoiseR = (-InoiseR(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
    G(ti + tii,:) = G(ti + tii -1,:) + dG + InoiseG(ti + tii -1,:);
    R(ti + tii,:) = R(ti + tii -1,:) + dR + InoiseR(ti + tii -1,:);
    I(ti + tii,:) = I(ti + tii -1,:) + dI + InoiseI(ti + tii -1,:);
    InoiseG(ti + tii,:) = InoiseG(ti + tii -1,:) + dInoiseG;
    InoiseI(ti + tii,:) = InoiseI(ti + tii -1,:) + dInoiseI;
    InoiseR(ti + tii,:) = InoiseR(ti + tii -1,:) + dInoiseR;
    
    G(ti + tii,G(ti + tii,:) < 0) = 0;
    R(ti + tii,R(ti + tii,:) < 0) = 0;
    I(ti + tii,I(ti + tii,:) < 0) = 0;
    
    if max(R(ti + tii,:)) >= thresh
        dummy = 1;
        rt = tii*dt;
        choice = find(R(ti+tii,:) == max(R(ti+tii,:)));
    end
end

%% 3rd stage: after decision
if tii == timeconstrain/dt
    rt = NaN;
    choice = NaN;
%     rt = tii*dt;
%     choice = find(R(ti+tii,:) == max(R(ti+tii,:)));
%     if length(choice) == 2
%         choice = randi(2);
%     end
elseif stoprule == 0 && tii < timeconstrain/dt
    V = zeros(sizeVraw);
    while (tii < timeconstrain/dt)
        tii = tii + 1;
%         if ti + tii > (presentt + stimdur)/dt
%             V = zeros(sizeVraw);
%         end
        % dI = (-I(ti + tii -1,:) + R(ti + tii -1,:) * 0) * dt/Tau(3);
        dI = (-I(ti + tii -1,:) + R(ti + tii -1,:) * beta) * dt/Tau(3);
        dG = (-G(ti + tii -1,:) + R(ti + tii -1,:) * w - I(ti + tii -1,:)) * dt/Tau(2);
        % dR = (-R(ti + tii -1,:) + (V + R(ti + tii - 1,:)*0)./(1+G(ti + tii -1,:))) * dt/Tau(1);
        dR = (-R(ti + tii -1,:) + (V + R(ti + tii - 1,:)*alpha)./(1+G(ti + tii -1,:))) * dt/Tau(1);
        dInoiseG = (-InoiseG(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
        dInoiseI = (-InoiseI(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
        dInoiseR = (-InoiseR(ti + tii -1,:).*dt./tauN + randn(sizeVraw).*sqrt(dt).*sigma.*dt./tauN);
        G(ti + tii,:) = G(ti + tii -1,:) + dG + InoiseG(ti + tii -1,:);
        R(ti + tii,:) = R(ti + tii -1,:) + dR + InoiseR(ti + tii -1,:);
        I(ti + tii,:) = I(ti + tii -1,:) + dI + InoiseI(ti + tii -1,:);
        InoiseG(ti + tii,:) = InoiseG(ti + tii -1,:) + dInoiseG;
        InoiseI(ti + tii,:) = InoiseI(ti + tii -1,:) + dInoiseI;
        InoiseR(ti + tii,:) = InoiseR(ti + tii -1,:) + dInoiseR;
        
        G(ti + tii,G(ti + tii,:) < 0) = 0;
        R(ti + tii,R(ti + tii,:) < 0) = 0;
        I(ti + tii,I(ti + tii,:) < 0) = 0;
    end
end
