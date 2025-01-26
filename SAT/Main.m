% Main function of generating figures for the paper: Tian & Shen, Speed
% accuracy tradeoff
[os, ~, ~] = computer;
if strcmp(os,'MACA64')
    Codedir = '/Users/bs3667/LDDM';
    rootdir = '/Users/bs3667/NYU Langone Health Dropbox/Shen Bo/Bo Shen Working files/Speed-accuracy-tradeoff';
elseif strcmp(os, 'PCWIN64')
    Codedir = 'C:\Users\Bo\Documents\LDDM';
    rootdir = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\Speed-accuracy-tradeoff';
end
addpath(genpath(fullfile(Codedir,'CoreFunctions'))); % inlcude the core functions
addpath(genpath(fullfile(Codedir,'utils'))); % inlcude utilities and functions

%% define parameters for simulation
c = [3.2 12.8, 25.6, 38.4 51.2]'/100; % percentage of coherence
scale0 = 250;
BR = 70;
a0 = 15;
BG = 0;
b0 = 1.1;
dt = .001;
Tau = ones(1,3)*.1;
thresh = 70;
%% define parameters for visualization
lwd = 2.0;
mksz = 18;
fontsize = 12;
mygray = flip(gray(length(c) + 2));
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
colorpalettergb = [239,71,111;255,209,102;6,214,160;17,138,178;7,59,76]/255;
colorpalettergb =[
    0.9373    0.2784    0.4353
    1.0000    0.8196    0.4000
    0.0235    0.8392    0.6275
    0.0667    0.5412    0.6980
    0.0275    0.2314    0.2980];
%%
aspect1 = [3.9,2.2]; % 16:9, for wide temporal dynamic
aspect2 = [3 3]; % for temporal dynamic
aspect3 = [2.8 2.54]; % 1:1 for phase plane
aspect4 = [3.302 2.54]; % 1:1 for heatmap with colorbar
aspect5 = [4.1 2.54]; % 1:1 for phase plane with legend outside
aspect6 = [2.2 4.0]; % for picked FR in fig5e
aspect7 = aspect3*.95; % for WTA surf
aspect8 = [2, 6.4]; % for the long format RT distribution fitting panels
aspect9 = [4.6 4.0]; % for fitted time course
aspect10 = [2.8 2.1]; % for fitted acc and RT
aspect11 = [2.9,.3]; % for timeline Fig7
aspect12 = [3.9,.6]; % for timeline Fig8
aspect13 = [3, 2.0]; % for dynamic Fig7
aspect14 = [2.41 3]; % for choice and RT panel
aspect15 = [3, 11]; % for combined Fig7
%% Figure 1 - Example dataset of Speed and Accuracy
% show the mean reaction time and accuracy difference, and the skewness difference 


%% Fig 2 - Example dynamics of different model
outdir = fullfile(rootdir, 'Fig2');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
plotdir = fullfile(outdir,'Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
datadir = fullfile(outdir,'SimRslts');
if ~exist(datadir,'dir')
    mkdir(datadir);
end
