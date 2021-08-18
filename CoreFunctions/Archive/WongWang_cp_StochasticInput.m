function [H, S, nu_wind, s_wind, rt, choice] = WongWang_cp_StochasticInput(Vraw, sgminput, rsmplFreq, a,b,d,gamma, JN, JAext, miu0, sigma, I0,...
    tauS, tauAMPA, unit, dur, dt, presentt, stimdur, thresh, initialvals, stoprule)
% input
cp1 = Vraw(1)/256;
cp2 = Vraw(2)/256;
sizeVraw = [1,2];

tauS = tauS/unit; % ms/unit
tauAMPA = tauAMPA/unit; % ms/unit
dt = dt/unit;
presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
timeconstrain = dur - presentt;
smoothwindow = .05/unit/dt;
if exist('initialvals','var')
    S(1,:) = initialvals(1,:);
    H(1,:) = initialvals(2,:);
    Inoise(1,:) = initialvals(3,:);
else
    S(1,:) = zeros(sizeVraw);
    H(1,:) = zeros(sizeVraw);
    Inoise = zeros(sizeVraw);
end
%% stage 0, before presenting, I = 0.
for ti = 2:(presentt/dt)
    eta = randn(1,2);
    
    I = [1,1]*0;
    x(1) = JN(1,1)*S(ti-1,1) - JN(1,2)*S(ti-1,2) + I0 + I(1) + Inoise(1);
    x(2) = JN(2,2)*S(ti-1,2) - JN(2,1)*S(ti-1,1) + I0 + I(2) + Inoise(2);
    H(ti,:) = (a.*x - b)./(1 - exp(-d.*(a.*x - b)));
    dS = (-S(ti-1,:)/tauS + (1-S(ti-1,:)).*H(ti-1,:)*unit*gamma)*dt;
    
    S(ti,:) = S(ti-1,:) + dS;
    S(ti,S(ti,:) < 0) = 0;
    %H(ti,H(ti,:) < 0) = 0;
    dInoise = (-Inoise.*dt./tauAMPA + eta.*sqrt(dt/tauAMPA*sigma^2)); %
    Inoise = Inoise + dInoise;
end
%% stage 1, presenting/decision
ti1st = ti;
dummy = 0;
rt = NaN;
choice = NaN;
tii = 0;
while tii < timeconstrain/dt && dummy == 0
    tii = tii + 1;
    if tii > stimdur/dt
        I = [0,0];
    else
        if mod(tii, 1/rsmplFreq/dt) == 1
            I = JAext*(miu0*[cp1, cp2]+randn(1,2)*sgminput);
        end
    end
    x(1) = JN(1,1)*S(ti1st+tii-1,1) - JN(1,2)*S(ti1st+tii-1,2) + I0 + I(1) + Inoise(1);
    x(2) = JN(2,2)*S(ti1st+tii-1,2) - JN(2,1)*S(ti1st+tii-1,1) + I0 + I(2) + Inoise(2);
    H(ti1st+tii,:) = (a*x - b)./(1 - exp(-d*(a*x - b)));
    dS = (-S(ti1st+tii-1,:)/tauS + (1-S(ti1st+tii-1,:)).*H(ti1st+tii-1,:)*unit*gamma)*dt;
    S(ti1st+tii,:) = S(ti1st+tii-1,:) + dS;
    S(ti1st+tii,S(ti1st+tii,:) < 0) = 0;
    %H(ti1st+tii,H(ti1st+tii,:) < 0) = 0;
    eta = randn(1,2);
    dInoise = (-Inoise.*dt/tauAMPA + eta.*sqrt(dt/tauAMPA*sigma^2)); % (-Inoise + eta*sqrt(dt*tauAMPA*sigma^2))*dt/tauAMPA;
    Inoise = Inoise + dInoise;
    postsynapFR = mean(H((ti1st + tii-smoothwindow+1):(ti1st + tii),:),1);
    if max(postsynapFR) >= thresh
        dummy = 1;
        rt = tii*dt;
        choice = find(postsynapFR == max(postsynapFR));
    end
end

%% stage 2, after decision
if tii == timeconstrain/dt
    rt = NaN;
    choice = NaN;
elseif stoprule == 0 && tii < timeconstrain/dt
    while (tii < timeconstrain/dt)
        tii = tii + 1;
        if tii > stimdur/dt
            I = [0,0];
        else
            if mod(tii, 1/rsmplFreq/dt) == 1
                I = JAext*(miu0*[cp1, cp2]+randn(1,2)*sgminput);
            end
        end
        eta = randn(1,2);
        
        x(1) = JN(1,1)*S(ti1st+tii-1,1) - JN(1,2)*S(ti1st+tii-1,2) + I0 + I(1) + Inoise(1);
        x(2) = JN(2,2)*S(ti1st+tii-1,2) - JN(2,1)*S(ti1st+tii-1,1) + I0 + I(2) + Inoise(2);
        H(ti1st+tii,:) = (a*x - b)./(1 - exp(-d*(a*x - b)));
        dS = (-S(ti1st+tii-1,:)/tauS + (1-S(ti1st+tii-1,:)).*H(ti1st+tii-1,:)*unit*gamma)*dt;
        S(ti1st+tii,:) = S(ti1st+tii-1,:) + dS;
        S(ti1st+tii,S(ti1st+tii,:) < 0) = 0;
        %H(ti1st+tii,H(ti1st+tii,:) < 0) = 0;
        dInoise = (-Inoise*dt/tauAMPA + eta*sqrt(dt/tauAMPA*sigma^2)); % (-Inoise + eta*sqrt(dt*tauAMPA*sigma^2))*dt/tauAMPA;
        Inoise = Inoise + dInoise;
        
    end
end

%% Calculating the mean rates and gating variables with sliding window
T_total = length(H);
time_wind = .050/unit/dt;  % Temporal window size for averaging
slide_wind = .005/unit/dt;  % Sliding step for window

nu_wind = [] ;
s_wind  = [] ;
nu_wind = [nu_wind; (mean(H(1:time_wind,:),1))] ;
s_wind  = [s_wind; (mean(S(1:time_wind,:),1))] ;

for t = 1:(T_total-time_wind)/slide_wind
    nu_wind = [nu_wind; (mean(H(slide_wind*t:slide_wind*t+time_wind,:),1))] ;
    s_wind  = [s_wind; (mean(S(slide_wind*t:slide_wind*t+time_wind,:),1))] ;
end