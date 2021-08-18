function [nu_wind, s_wind, rt, choice, H, S, Icourse] = wong06_RndInput_Gaba(Vinput,Gaba, miu0,sgm,sgmInput,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule)
%%%%%%%%%%%%
% adapted version of the appendix function in Wong & Wang, 2006
% Created by Bo Shen, NYU, 2019
%%%%%%%%%%%%
sizeVinput = size(Vinput);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end

cp = (Vinput/256);
% if max(cp) > 2
%     warning('The input scale exeeds the orignal settings in Wong & Wang, 2006');
% end
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
JN(1,2) = JN(1,2)*Gaba;
JN(2,1) = JN(2,1)*Gaba;
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
slide_wind = .005/dt/unit;  % Sliding step for window, 5 msec
% time line
total_time_steps = round(dur/dt);
onset_of_stimuli = round(presentt/dt);
stim_duration = round(stimdur/dt);
rt = NaN;
choice = NaN;
%% get initial values
Hi = initialvals(1,:);
Si = initialvals(2,:);
for vi = 1:sizeVinput(2)
    syms x;
    if Hi(vi) < 10^-5
        xi(vi) = 0;
    else
        s1 = solve((a*x - b)./(1 - exp(-d*(a*x - b))) == Hi(vi), x);
        xi(vi) = (double(vpa(s1)));
    end
end
%% stablizing noise
Inoise(1,:) = zeros(1,sizeVinput(2));
stablizetime = round(.2/unit/dt);
for kk = 1:stablizetime
    dInoise = (-Inoise/tauAMPA*dt + randn(sizeVinput)*sqrt(dt/tauAMPA*sgm^2)); %
    Inoise = Inoise + dInoise;
    H = (a*(xi+Inoise) - b)./(1 - exp(-d*(a*(xi+Inoise) - b)));
    H(H < 0) = 0;
    loci = mod(kk,time_wind) + (mod(kk,time_wind) == 0)*time_wind;
    Hbuffer(loci,:) = real(H);
    Sbuffer(loci,:) = real(Si);
end
H = real(mean(Hbuffer,1));
S = H*gamma*tauS./(H*gamma*tauS+1);
nu_wind = H;
s_wind = S;
ti = 1;
I = JAext*miu0*(cp*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration) + randn(size(cp))*sgmInput);
Icourse = I;
%% simulation begin
for ti = 2:total_time_steps
    % update nerual firing rates and synaptic activities
    if (mod(ti*dt, .05) == 0)
        I = JAext*miu0*(cp*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration) + randn(size(cp))*sgmInput);
    end
    Icourse(ti,:) = I;
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

%% Calculating the mean rates and gating variables with sliding window
% T_total = length(H);
% nu_wind = mean(H(1:time_wind,:),1);
% s_wind  = mean(S(1:time_wind,:),1) ;
% for t = 1:((T_total-time_wind)/slide_wind)
%     nu_wind(t+1,:) = mean(H((1:time_wind)+slide_wind*t,:),1);
%     s_wind(t+1,:)  = mean(S((1:time_wind)+slide_wind*t,:),1);
% end