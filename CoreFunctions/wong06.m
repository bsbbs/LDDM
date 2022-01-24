function [choice, rt, nu_wind, s_wind, H, S] = wong06(cp,miu0,sgm,I0,JN,...
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
% time unit
unit = 1; % set as 1 if the unit of time related variables is second.
% change it to .001 if the unit of time related varibales is msec. 

% connectivity strength
JAext = 5.2 * 10^-4; %nA/Hz

%% define parameters
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
posttask_steps = round(dur/dt);
onset_of_stimuli = round(presentt/dt);
stim_duration = round(stimdur/dt);
offset_of_stimuli = onset_of_stimuli + stim_duration;
rt = NaN;
choice = NaN;
%% stablizing noise for 200 ms
Inoise = zeros(sizeVinput);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    % update noise
    dInoise = (-Inoise/tauAMPA*dt + randn(sizeVinput)*sqrt(dt/tauAMPA*sgm^2)); %
    Inoise = Inoise + dInoise;
end
%% get initial values at ti = 0
% drop in buffer to sliding window
Hbuffer = repmat(initialvals(1,:),time_wind,1);
Sbuffer = repmat(initialvals(2,:),time_wind,1);
nu_wind(1,:) = mean(Hbuffer,1);
s_wind(1,:) = mean(Sbuffer,1);
S(1,:) = initialvals(2,:);
H(1,:) = initialvals(1,:);
%% simulation begin
t_stamp = 1;
for ti = 0:posttask_steps
    % update nerual firing rates based on inputs
    I = JAext*miu0*cp*(ti >= onset_of_stimuli & ti < offset_of_stimuli);
    x = (JN*S(ti+t_stamp,:)')' + I0 + I + Inoise;
    H(ti+t_stamp+1,:) = (a*x - b)./(1 - exp(-d*(a*x - b)));
    % update synaptic activities based on firing rates
    dS = (-S(ti+t_stamp,:)/tauS + (1-S(ti+t_stamp,:)).*H(ti+t_stamp,:)*unit*gamma)*dt;
    S(ti+t_stamp+1,:) = S(ti+t_stamp,:) + dS;
    % update noise
    dInoise = (-Inoise/tauAMPA*dt + randn(sizeVinput)*sqrt(dt/tauAMPA*sgm^2)); %
    Inoise = Inoise + dInoise;
    
    % To ensure firing rates are always positive (noise may cause negative)
    H(ti+t_stamp+1,H(ti+t_stamp+1,:) < 0) = 0;
    % drop in buffer to sliding window
    loci = mod(ti,time_wind) + (mod(ti,time_wind) == 0)*time_wind;
    Hbuffer(loci,:) = H(ti+t_stamp,:);
    Sbuffer(loci,:) = S(ti+t_stamp,:);
    nu_wind(ti+t_stamp,:) = mean(Hbuffer,1);
    s_wind(ti+t_stamp,:) = mean(Sbuffer,1);
    % threshold detecting
    if ti >= onset_of_stimuli && isnan(rt)
        postsynapFR = nu_wind(ti+t_stamp,:);
        if max(postsynapFR) >= thresh
            rt = (ti - onset_of_stimuli)*dt;
            choice = find(postsynapFR == max(postsynapFR));
            if stoprule == 1
                break;
            end
        end
    end
end