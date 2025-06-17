%% %%%%%%%% 
% Created by Bo Shen Ph.D., NYU Grossman School of Medicine, May, 2025
% Simulation on STDP of E-I network when receiving noisy inputs
%%%%%%%%%%

%% define paths
Homedir = 'C:\Users\Bo\Documents\GitHub';
% Homedir = '~';
% addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'LDDM','Froemke','Revision1'));
addpath(fullfile(Homedir,'LDDM','Froemke','utils'));
Drpbx = "C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project0";
out_dir = fullfile(Drpbx, 'Revision1');
if ~exist("out_dir",'dir')
    mkdir(out_dir);
end
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
%% input parameters
c_rprsnt = .064;
sgmInput_rprsnt = 1/30;
c_choice = .032;
sgmInput_choice = .75;
BG = 30; % background inputs, Hz
BR = 30; % background inputs, Hz
%% output/activity parameters
thresh = 70; % Hz
eqlb = 32; % mean activities during equilibrium when D is off and with equal inputs


%% Simulation start...
%% Figure 3: B-F, Single input network with plasticity and two types of noise
Singleton;

%% Figures 4&5, Binary choice with E plasticity only and E-I plasticity, and two types of noise
Binary;
