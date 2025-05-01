% Demo

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
fontsize = 14;
mksz = 25;
lwd = 1.5;
cp = [.032, .064,.128, .256];
mygray = flip(gray(numel(cp) + 2));
myred = mygray; myred(:,1) = 1;
myblue = mygray; myblue(:,3) = 1;
%% parameters for simulation
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauD = .1; % second
Tau = [tauR, tauG, tauD];
w0 = 1;
a0 = 0;
b0 = 1.1; % 1.1;
c_rprsnt = .064;
sgmInput_rprsnt = 1/3;
c_choice = .032;
sgmInput_choice = .75;
predur = 0;
presentt = dt;
thresh = 70; % Hz
stoprule = 1;
eqlb = 32;
BG = 30;
BR = 30;

%%
