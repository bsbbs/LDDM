rng('shuffle');
% parameter in H
a = 270; %VnC^-1
b = 108; %Hz
d = .154; %s
gamma = .641;

% time constant
tauS = .1; % sec
tauAMPA = .002; % sec
unit = 1; % secs
%unit = .001; % msecs

% connectivity strength
JN = [.2609 -.0497
    -.0497   .2609]; % nA
% JN = [.1561 .0264;
%     .0264   .1561]; %nA
JAext = 5.2 * 10^-4; %nA/Hz
% JAext = .2243 * 10^-3; %nA/Hz

% input level 
miu0 = 30; % Hz
sigma_dflt = .02; % nA
sgm_dflt = sigma_dflt;
% sigma = .007; % nA
I0 = .3255; % nA
% I0 = .2346; % nA

% time simulation parameter
dt_dflt = .0001; % sec
dur_dflt = 3; % sec
presentt_dflt = .5; % sec
stimdur_dflt = 1.5; % sec

thresh_dflt = 15;
stoprule_dflt = 1;
initialvals_dflt = [2 2;.1 .1; sigma_dflt*randn, sigma_dflt*randn]; % H, S, and noise;
%Vraw = [130 40];
%cp = (Vraw(1) - Vraw(2))/mean(Vraw)/2;
cp_dflt = 6.4/100;