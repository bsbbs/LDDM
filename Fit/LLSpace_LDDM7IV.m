%% set path
Homedir = '~/Documents';
addpath(fullfile(Homedir, 'LDDM','CoreFunctions'));
addpath(fullfile(Homedir, 'LDDM','utils'));
addpath(fullfile(Homedir, 'LDDM', 'Fit','SvrCode'));

Glgdir = '/Volumes/GoogleDrive/My Drive/LDDM';
out_dir = fullfile(Glgdir,'Fit/Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU/LLSpace');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end


%% load empirical data
dataBhvr = LoadRoitmanData(fullfile(Glgdir,'RoitmanDataCode'));

%% the best fitting parameters
% a, b, noise, scale, tauR, tauG, tauD, nLL
bestparams = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132 16539.138186];
names = {'alpha','beta','sgm','S','tauR','tauG','tauD'};
%% Check the space of alpha and beta
ivec = linspace(0,100,41); % [0, 10.^[-1:.1:3]]; %10.^[-1:.01:3];
jvec = linspace(0,4,41); %linspace(0,4,401);
idx = [1, 2];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir);
Visualization(out_dir,filename, names, idx);
h = figure;
contour(nLLmat);
%% Check the space of noise and scale
ivec = linspace(0,80,41); % noise, best fit = 25.36 
jvec = linspace(10, 8010, 41); % 10.^[-1:.1:3]*9; % scale, best fit = 3251
idx = [3, 4];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir);
Visualization(out_dir,filename, names, idx);
%% Check the space of tauR and tauG
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauR, best fit = .18
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41); % tauG, best fit = .22
idx = [5, 6];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir);
Visualization(out_dir,filename, names, idx);
% Check the space of tauR and tauD
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauR, best fit = .18
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41);% tauD, best fit = .32
idx = [5, 7];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir);
Visualization(out_dir,filename, names, idx);
%% Check the space of tauG and tauD
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauG, best fit = .22
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41); % tauD, best fit = .32
idx = [6, 7];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir);
Visualization(out_dir,filename, names, idx);
%% functions
function [filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, dataBhvr, out_dir)
iN = length(ivec);
jN = length(jvec);
filename = sprintf('LLSpace_%s%i_%s%i', names{idx(1)}, iN, names{idx(2)}, jN);
if ~exist(fullfile(out_dir,[filename, '.mat']), 'file')
    nLLmat = NaN(iN, jN);
    for i = 1:iN
        ival = ivec(i);
        fprintf('%02.1f %%\n', i/iN*100);
        parfor j = 1:jN
            params = bestparams;
            params(idx) = [ival, jvec(j)];
            [nLL, ~,~,~,~,~] = LDDMFitBhvr7ParamsIV_QMLE_GPU(params, dataBhvr);
            nLLmat(i,j) = nLL;
        end
    end
    save(fullfile(out_dir,[filename, '.mat']), 'nLLmat', 'ivec','jvec','idx');
else
    load(fullfile(out_dir,[filename, '.mat']));
end
end

function Visualization(out_dir,filename, names, idx)
% Visuliazation parameters
lwd = 2.0;
mksz = 18;
fontsize = 14;
load(fullfile(out_dir,[filename, '.mat']));
h = figure; hold on;
s = surf(nLLmat,'EdgeColor','none');
minnLL = min(nLLmat(:));
[r, c] = find(nLLmat == minnLL);
plot3(c,r,minnLL, 'rx', 'MarkerSize', mksz, 'LineWidth', lwd);
ylim([1,length(ivec)]);
xlim([1,length(jvec)]);
xticks(linspace(1,length(jvec),5));
xticklabels(jvec(linspace(1,length(jvec),5)));
yticks(linspace(1,length(ivec),5));
yticklabels(ivec(linspace(1,length(ivec),5))); % {'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel(names{idx(2)});
ylabel(names{idx(1)});
clb = colorbar;
ylabel(clb, 'Negative Log likelihood');
view(0,90);
grid on;
savefigs(h, filename, out_dir, fontsize - 2, [2.8 2.54]*.95);
end