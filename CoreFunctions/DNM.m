function [R, G] = DNM(Vinput, w, sgm, Tau, presentt, dur, dt, initialvals, stimdur)
tauN =0.002; % time constant for Ornstein-Uhlenbeck process of noise
%% preparation
total_time_steps = round(dur/dt);
onset_of_stimuli = round(presentt/dt);
stim_duration = round(stimdur/dt);
sizeVinput = size(Vinput);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end
%% simulation begin
G(1,:) = initialvals(2,:);
R(1,:) = initialvals(1,:);
InoiseG = zeros(sizeVinput);
InoiseR = zeros(sizeVinput);
for ti = 2:total_time_steps
    % update R, G
    dG = (-G(ti-1,:)' + w*R(ti-1,:)')/Tau(2)*dt;
    dR = (-R(ti-1,:)' + Vinput'*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration)./(1+G(ti-1,:)'))/Tau(1)*dt; 
    G(ti,:) = G(ti-1,:) + dG' + InoiseG(ti-1,:);
    R(ti,:) = R(ti-1,:) + dR' + InoiseR(ti-1,:);
    % update noise
    dInoiseG = (-InoiseG(ti-1,:).*dt./tauN + randn(sizeVinput).*sqrt(dt).*sgm.*dt./tauN);
    dInoiseR = (-InoiseR(ti-1,:).*dt./tauN + randn(sizeVinput).*sqrt(dt).*sgm.*dt./tauN);
    InoiseG(ti,:) = InoiseG(ti-1,:) + dInoiseG;
    InoiseR(ti,:) = InoiseR(ti-1,:) + dInoiseR;
    % setting lower boundary, forcing neural firing rates to be non-negative
    G(ti,G(ti,:) < 0) = 0;
    R(ti,R(ti,:) < 0) = 0;
end
