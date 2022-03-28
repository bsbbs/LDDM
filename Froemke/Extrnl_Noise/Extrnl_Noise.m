%% define paths
Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(fullfile(Homedir,'Documents','LDDM','Froemke/Extrnl_Noise'));
cd('G:\My Drive\LDDM\Froemke\Extrnl_Noise');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Froemke');
plotdir = fullfile('./Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end

%% E-E only, self-excitation
 
%% E-E only, input

%% E-E only, self-excitation & input

%% I-E only, G -> R

%% Both E-E & I-E

