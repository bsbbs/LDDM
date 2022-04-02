%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(genpath(fullfile(Homedir,'Documents','LDDM','Froemke')));
out_dir = 'G:\My Drive\LDDM\Froemke\Extrnl_Noise';
if ~exist("out_dir",'dir')
    mkdir(out_dir);
end
cd(out_dir);
plotdir = fullfile('./Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end
%% parameters for visulization
fontsize = 14;
mksz = 5;
lwd = 1.5;
%% parameters for simulation
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];
w0 = 1;
a0 = 15;
b0 = 1.8;
predur = 0;
presentt = dt;
thresh = 70; % Hz
stoprule = 1;
eqlb = 32;
scale = 2*w0*eqlb.^2 + (1-a0).*eqlb;
Vprior = [1,1]*scale;
%% General property 1. representation dynamic
task = 'VR_RT';
stimdur = Inf;
triggert = Inf;
dur = 1.5; % second
c = [.064]'; % Coherence
c1 = [1 + c];
c2 = [1 - c];
% Dynamics 
sgm = .01;
% sgmInput = 150; for OU process
sgmInput = 0%1/3;
Vinput = [c1, c2]*scale;
a = a0*eye(2);
b = b0*eye(2);
w = w0*ones(2);
iSTDP = 1;
% initialvals = ones(3,2)*10; %(a-1)./sum(w,2)'.*[1,1; sum(w,2)'; 0, 0];
initialvals = (a0-1)./(2*w0).*[1, 1; w0*2, w0*2; 0, 0]/.5;
initialvals = [1,1;2-b0,2-b0;b0,b0]*eqlb;
h = figure;
filename = sprintf('timeCourse_%s_sgm%2.2f_sgmInput%.2f',task, sgm, sgmInput);
% Value representation task
rng(5);
[choice, rt, R, G, I, Vcourse] = LDDM_STDP_RndInput(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
subplot(2,1,1); hold on;
x = (1:numel(R(:,1)))*dt;
plot(x, R(:,1), 'r-', 'LineWidth', lwd);
plot(x, R(:,2), 'b-', 'LineWidth', lwd);
xlim([-50,dur/dt]);
ylim([-5,thresh]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
legend({'R_1','R_2'},'Box','off');
savefigs(h, filename, plotdir, fontsize, [2.9,4]);

% Reaction-time task
stimdur = Inf;
triggert = presentt;
dur = 1.5; % second
rng(8);
[choice, rt, R, G, I, Vcourse] = LDDM_STDP_RndInput(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
subplot(2,1,2); hold on;
x = (1:numel(R(:,1)))*dt;
plot(x, R(:,1), 'r-', 'LineWidth', lwd);
plot(x, R(:,2), 'b-', 'LineWidth', lwd);
%xlim([-50,dur/dt]);
ylim([-5,thresh]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
savefigs(h, filename, plotdir, fontsize, [2.9,4]);
%% General property 2. predicted choice behavior


%% E-E only, self-excitation
 
%% E-E only, input

%% E-E only, self-excitation & input

%% I-E only, G -> R

%% Both E-E & I-E

