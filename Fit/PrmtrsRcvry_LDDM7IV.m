%% set path
Homedir = 'C:\Users\Bo\Documents';
% Homedir = '~/Documents';
addpath(fullfile(Homedir, 'LDDM','CoreFunctions'));
addpath(fullfile(Homedir, 'LDDM','utils'));
addpath(fullfile(Homedir, 'LDDM', 'Fit','SvrCode'));

Glgdir = 'G:\My Drive\LDDM';
% Glgdir = '/Volumes/GoogleDrive/My Drive/LDDM';
% addpath(genpath(fullfile(Glgdir, 'Fit','bads-master')));
addpath(genpath(fullfile(Glgdir, 'Fit','bads')));
out_dir = fullfile(Glgdir,'Fit/Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU/PrmtrsRcvry');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end


%% the best fitting parameters
% a, b, noise, scale, tauR, tauG, tauD, nLL
bestparams = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132 16539.138186];
names = {'alpha','beta','sgm','S','tauR','tauG','tauD'};
latexnames = {'\alpha','\beta','\sigma','S','\tau_R','\tau_G','\tau_D'};
%% generating 'fake' data
% load empirical data - not really needed
dataBhvr = LoadRoitmanData(fullfile(Glgdir,'Fit','RoitmanDataCode'));

[~, ~, ~, ~, rtmat, choicemat] = LDDMFitBhvr7ParamsIV_QMLE_GPU(bestparams, dataBhvr, 102400);

%% calculate indicators
SimDataQp = Load_SimData(rtmat, choicemat);
%% Check the space of alpha and beta
ivec = linspace(0,100,41); % [0, 10.^[-1:.1:3]]; %10.^[-1:.01:3];
jvec = linspace(0,4,41); %linspace(0,4,401);
idx = [1, 2];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, SimDataQp, out_dir);
Visualization(out_dir,filename, names, idx, 5, latexnames);
%% Check the space of noise and scale
ivec = linspace(0,80,41); % noise, best fit = 25.36 
jvec = linspace(10, 8010, 41); % 10.^[-1:.1:3]*9; % scale, best fit = 3251
idx = [3, 4];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, SimDataQp, out_dir);
Visualization(out_dir,filename, names, idx, 5, latexnames);
%% Check the space of tauR and tauG
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauR, best fit = .18
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41); % tauG, best fit = .22
idx = [5, 6];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, SimDataQp, out_dir);
Visualization(out_dir,filename, names, idx, 4, latexnames);
%% Check the space of tauR and tauD
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauR, best fit = .18
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41);% tauD, best fit = .32
idx = [5, 7];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, SimDataQp, out_dir);
Visualization(out_dir,filename, names, idx, 4, latexnames);
%% Check the space of tauG and tauD
ivec = 10.^[-2:.1:1]; %linspace(.025,1.05,42); % tauG, best fit = .22
jvec = 10.^[-2:.1:1]; %linspace(.025,1.025,41); % tauD, best fit = .32
idx = [6, 7];
[filename, nLLmat] = SpaceCheck(bestparams, names, ivec, jvec, idx, SimDataQp, out_dir);
Visualization(out_dir,filename, names, idx, 4, latexnames);

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

%%
function Visualization(out_dir,filename, names, idx, nticks, latexnames)
% Visuliazation parameters
lwd = 2.0;
mksz = 18;
fontsize = 14;
load(fullfile(out_dir,[filename, '.mat']));
h = figure;
hold on;
minnLL = min(nLLmat(:));
[r, c] = find(nLLmat == minnLL);
contour(-nLLmat,60);
plot3(c,r,-minnLL*1.01, 'rx', 'MarkerSize', mksz/3, 'LineWidth', lwd);
xticks(linspace(1,length(jvec),nticks));
xticklabels(jvec(linspace(1,length(jvec),nticks)));
yticks(linspace(1,length(ivec),nticks));
yticklabels(ivec(linspace(1,length(ivec),nticks)));
xlabel(latexnames{idx(2)},'FontAngle', 'italic');
ylabel(latexnames{idx(1)},'FontAngle', 'italic');
clb = colorbar;
colormap('turbo');
ylabel(clb, 'Log likelihood');
savefig(h, fullfile(out_dir,['Contour',filename]));
savefigs(h, ['Contour',filename], out_dir, fontsize - 2, [4 3]);
% close(h);
end

%% refit
% % Define optimization starting point and bounds
% %     a,    b, noise, scale, Tau
% LB = [0    0.1   .1    .1*256 [.001,.001,.001]];
% UB = [70   3	100  20*256 [1,1,1]];
% PLB = [15  .9	5    1*256 [.01 .01 .01]];
% PUB = [60   1.7	40   8*256 [.2 .2 .2]];
% 
% % Randomize initial starting point inside plausible box
% x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;
% 
% % likelihood function
% % parpool(6);
% nLLfun = @(params) LDDMFitBhvr7ParamsIV_QMLE_GPU(params, SimDataQp);
% [fvalbest,~,~] = nLLfun(x0)
% fprintf('test succeeded\n');
% myCluster.NumWorkers = 6;
% sortNum = 1;
% % change starting points
% Collect = [];
% parfor i = 1:myCluster.NumWorkers*8
%     !ping -c 1 www.amazon.com
%     t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7 + i*10^5;
%     %num2str(t);
%     rng(t);
%     
%     % Randomize initial starting point inside plausible box
%     x0 = rand(1,numel(LB)) .* (PUB - PLB) + PLB;
%     dlmwrite(fullfile(out_dir,'x0List.txt'),[sortNum, i, t, x0],'delimiter','\t','precision','%.6f','-append');
%     % fit
%     options = bads('defaults');     % Default options
%     options.Display = 'iter';
%     % For this optimization, we explicitly tell BADS that the objective is
%     % noisy (it is not necessary, but it is a good habit)
%     options.UncertaintyHandling = true;    % Function is stochastic
%     % specify a rough estimate for the value of the standard deviation of the noise in a neighborhood of the solution.
%     options.NoiseSize = 2.7;  % Optional, leave empty if unknown
%     % We also limit the number of function evaluations, knowing that this is a
%     % simple example. Generally, BADS will tend to run for longer on noisy
%     % problems to better explore the noisy landscape.
%     % options.MaxFunEvals = 3000;
%     
%     % Finally, we tell BADS to re-evaluate the target at the returned solution
%     % with ** samples (10 by default). Note that this number counts towards the budget
%     % of function evaluations.
%     options.NoiseFinalSamples = 20;
%     [xest,fval,~,output] = bads(nLLfun,x0,LB,UB,PLB,PUB,[],options);
%     dlmwrite(fullfile(out_dir,'RsltList.txt'),[sortNum, i, t, xest fval],'delimiter','\t','precision','%.6f','-append');
%     
%     Collect(i).rndseed = t;
%     Collect(i).x0 = x0;
%     Collect(i).xest = xest;
%     Collect(i).fval = fval;
%     Collect(i).output = output;
%     
% end
% t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8 + sortNum*10^7 + i*10^5;
% save(fullfile(out_dir,sprintf('CollectRslts%i.mat',t)),'Collect');

%% Visualization
% load('/Volumes/GoogleDrive/My Drive/LDDM/Fit/Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU/PrmtrsRcvry/NewBADS_CollectRslts55536850.mat');
% filename = 'NewBADS_CollectRslts55536850';
% h = Vslz(Collect);
% savefig(h, fullfile(out_dir,filename));
% %% Visualization
% function h = Vslz(Collect)
% sz = length(Collect);
% Pls = nan(8,sz);
% for i = 1:sz
%     Pls(1:7,i) = Collect(i).xest;
%     Pls(8,i) = Collect(i).fval;
% end
% h = figure; hold on;
% scatter3(Pls(1,:), Pls(2,:), Pls(8,:), 190, Pls(8,:), 'Marker','.');
% loc = find(Pls(8,:) == min(Pls(8,:)));
% scatter3(Pls(1,loc), Pls(2,loc), Pls(8,loc), 190, 'r', 'Marker','*');
% xlabel('alpha');
% ylabel('beta');
% zlabel('nLL');
% grid on;
% view([20, 10]);
% end


