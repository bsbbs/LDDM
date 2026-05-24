% Simulation for the NSFC grantL: Noise in the aging brain, mouse tracking study.
%% define paths
Homedir = 'C:\Users\Bo\Documents\GitHub';
% Homedir = '~';
% addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'LDDM','NSFCgrant-Noise&Aging'));
addpath(fullfile(Homedir,'LDDM','utils'));
Drpbx = "C:\Users\Bo\Dropbox\Application - Grant\面上-2026-noise";
% Drpbx = "/Users/bs3667/NYU Langone Health Dropbox/Shen Bo/Bo Shen Working files/STDP_Project0";
out_dir = Drpbx;
plotdir = fullfile(out_dir,'Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = fullfile(out_dir,'SimRslts');
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end
%% parameters for visulization
fontsize = 10;
mksz = 25;
lwd = 1.5;
cp = [.032, .064,.128, .256];
cp = [.032, .128, .256, .512];
mygray = flip(gray(numel(cp) + 2));
myred = mygray; myred(:,1) = 1;
myblue = mygray; myblue(:,3) = 1;
mygreen = mygray; mygreen(:,2) = 1;
%% kinematic parameters
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauD = .1; % second
Tau = [tauR, tauG, tauD];
%% weighting parameters
w0 = 1;
a0 = 0;
b0 = 1.0; % 1.1;

% Simulation