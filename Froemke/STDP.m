% Temporal dynamics of STDP
%% load functions
addpath('./functions');
addpath('utils');

%% simulating a temporal dyamics of connection weights over time
c = .128;
Vprior = 256 *[1, 1];
Vinput = 256*[1+c, 1-c];
STDP_V = 'dynamic';
STDP_a = 'dynamic';
STDP_G = 1;
STDP_D = 1;
w = ones(2);
a = eye(2)*15;
b = eye(2)*1.1;
sgm = .2;
sgmInput = .75;
Tau = [.1, .1, .1];
predur = 0;
dur = 10;
dt = .001;
presentt = 0;
triggert = 0;
initialvals = zeros(3,2);
stimdur = dur;
stoprule = 1;
thresh = 70;
[choice, rt, R, G, D, Vcourse, wV, wa, wG, wD] = LDDM_RndInput_STDPDynmc(Vprior, Vinput, STDP_V, STDP_a, STDP_G, STDP_D, w, a, b,...
    sgm, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%%
h = figure;
hold on;
plot(wV);
plot(wa);
plot(wG);
plot(wD);
h = figure;
hold on;
plot(R);
plot(G);
plot(D);
%% functions

