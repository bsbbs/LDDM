initialvals = zeros(3,2);
presentt = .5;
stimdur = .7;
triggert = 2.2;
dur = 3;
dt = .001;
Tau = [.1 .1 .1];
thresh = 70;
w = [1 1;
    1 1];
alpha = [1 0
    0 1]*10; % self-excitation, [a11, a12; a21, a22]
beta = [1 0
    0 1]*2; % R to I,  [b11, b12; b21, b22]
sigma = 0;
Vraw = [130 40];
stoprule = 0;
