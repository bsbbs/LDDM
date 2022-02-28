Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(genpath(fullfile(Homedir,'Documents','LDDM','RTSkewness')));
cd('G:\My Drive\LDDM\RTSkewness');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/RTSkewness');
out_dir = './';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
dataDynmc = load('./Data/Data.mat');
dataBhvr = LoadRoitmanData('../RoitmanDataCode');
randseed = 24356545;
rng(randseed);

params = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132	16539.138186];
name = sprintf('a%2.2f_b%1.2f_sgm%2.1f_scale%4.1f_tau%1.2f_%1.2f_%1.2f',params(1:7));
if ~exist(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),'file')
    tic;
    [nLL, Chi2, BIC, AIC, rtmat, choicemat] = LDDMFitBhvr7ParamsIV_QMLE_GPU(params, dataBhvr);
    save(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)),...
        'rtmat','choicemat','params','nLL','Chi2','AIC','BIC');
    toc
else
    load(fullfile(plot_dir,sprintf('PlotData_%s.mat',name)));
end

