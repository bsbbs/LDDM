function [nu_wind, s_wind, rt, choice, H, S] = wong06(cp,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule)
%%%%%%%%%%%%
% adapted version of the appendix function in Wong & Wang, 2006
% Created by Bo Shen, NYU, 2019
%%%%%%%%%%%%
sizeVinput = size(cp);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end

%% set parameters
% parameter in H
a = 270; %VnC^-1
b = 108; %Hz
d = .154; %s
% gamma = .641;
unit = 1; % set as 1 if the unit of time related variables is second.
% change it to .001 if the unit of time related varibales is msec. 

% connectivity strength
% JN = eye(sizeVinput(2))*(.2609+.0497) - ones(sizeVinput(2))*.0497; %nA
JAext = 5.2 * 10^-4; %nA/Hz

%% preparation
tauS = tauS/unit;
tauAMPA = tauAMPA/unit;
dt = dt/unit;
presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
% sliding window
time_wind = .050/dt/unit;  % Temporal window size for averaging, 50 msec
% not used, slided at each step of dt. slide_wind = .005/dt/unit;  % Sliding step for window, 5 msec
% time line
total_time_steps = round(dur/dt);
onset_of_stimuli = round(presentt/dt);
stim_duration = round(stimdur/dt);
rt = NaN;
choice = NaN;
%% get initial values
H = initialvals(1,:);
S = initialvals(2,:);
nu_wind = H;
s_wind = S;
%% simulation begin
for ti = 2:total_time_steps
    % update nerual firing rates and synaptic activities
    I = JAext*miu0*cp*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration);
    x = (JN*S(ti-1,:)')' + I0 + I + Inoise;
    H(ti,:) = (a*x - b)./(1 - exp(-d*(a*x - b)));
    dS = (-S(ti-1,:)/tauS + (1-S(ti-1,:)).*H(ti-1,:)*unit*gamma)*dt;
    S(ti,:) = S(ti-1,:) + dS;
    % update noise
    dInoise = (-Inoise/tauAMPA*dt + randn(sizeVinput)*sqrt(dt/tauAMPA*sgm^2)); %
    Inoise = Inoise + dInoise;
    % To ensure firing rates are always positive (noise may cause negative)
    H(ti,H(ti,:) < 0) = 0;
    loci = mod(kk+1+ti,time_wind) + (mod(kk+1+ti,time_wind) == 0)*time_wind;
    Hbuffer(loci,:) = H(ti,:);
    Sbuffer(loci,:) = S(ti,:);
    nu_wind(ti,:) = mean(Hbuffer,1);
    s_wind(ti,:) = mean(Sbuffer,1);
    % threshold detecting
    if ti >= onset_of_stimuli && isnan(rt)
        postsynapFR = nu_wind(ti,:);
        if max(postsynapFR) >= thresh
            rt = (ti - onset_of_stimuli)*dt;
            choice = find(postsynapFR == max(postsynapFR));
            if stoprule == 1
                break;
            end
        end
    end
end