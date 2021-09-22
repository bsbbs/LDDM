%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath('../../RecurrentModel');
numNode = 1;
[sortNum, myCluster] = RndCtrl(numNode);
mypool = parpool(myCluster, myCluster.NumWorkers);
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = fullfile('rslts','LDDM_ABT');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
%% initalize parameters as in the paper of Wong & Wang, 2006
task = 'FDCalib';
% simulation parameters
dt = .001; % second
presentt = dt;
w = ones(3);
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];
sgm = 5;
dur = 7; % second
stimdur = 2.5;
thresh = 70; % Hz
% scale = 25;
initialvals = [4, 4, 4; 12, 12, 12; 0, 0, 0]; % for R, G, and I variables
stoprule = 1;
sims = 10240*2; % number of iterations
% plot conditional ratio as a function of V3

%%
aspect = [9, 9];
plti = 0;
panel = 0;
avec = 0:1:95; %[10,30,50,60,70,75,80,85,90]; % 0,15,20,96
bvec = [0.6:.2:4.4];
triggert = 2:.1:3;
c = .256;
c1 = 1 + c;
c2 = ones(size(c1));
c3 = 0:.1:2;
[triggertmat, cp3] = meshgrid(triggert, c3);
cp1 = repmat(c1, size(triggertmat));
cp2 = repmat(c2, size(triggertmat));
figname = sprintf('IIA_LDDM_%s_stimdur%1.1f_sgm%2.1f_%iA%iB%iT',task,stimdur,sgm,numel(avec),numel(bvec),numel(triggert));
%choicemat = nan([numel(avec),numel(bvec),size(cp1),sims]);
%rtmat = nan([numel(avec),numel(bvec),size(cp1),sims]);
cndratiomat = nan([numel(avec),numel(bvec),size(cp1)]);
rtallmat = nan([numel(avec),numel(bvec),size(cp1)]);
rtV1V2mat = nan([numel(avec),numel(bvec),size(cp1)]);
slopemat = nan([numel(avec),numel(bvec),2,numel(triggert)]);
for ai = 1:length(avec)
    scale = 3*32^2 + (1-avec(ai))*32;
    Vinput.V1 = cp1*scale;
    Vinput.V2 = cp2*scale;
    Vinput.V3 = cp3*scale;
    a = eye(3)*avec(ai);
    parfor bi = 1:length(bvec)
        b = bvec(bi);
        % plti = plti + 1;
        %         if mod(plti,9) == 1
        %             panel = panel + 1;
        %             h = figure;
        %             figname = sprintf('IIA_LDDM_%s_stimdur%1.1f_ABT_ap%i',task,stimdur,panel);
        %         end
        % simulation
        %filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_a%2.1f_w%1.1f_b%1.2f_scale%2.1f_%iact%1.2f-%1.2f_sgm%2.1f_sim%i',...
            %task,length(c1),length(unique(c2)),length(unique(c3)),a(1,1),1,b,scale,numel(triggert),triggert(1),triggert(2),sgm,sims);
        %simrslt = fullfile(outdir,[filename, '.mat']);
        %if ~exist(simrslt,'file')
            [choice, rt] = LDDM_GPU3(Vinput, w, a, eye(3)*b, sgm, Tau, dur,...
                dt, presentt, triggertmat, thresh, initialvals, stimdur, stoprule, sims);
            %save(simrslt,'choice','rt');
        %else
            %load(simrslt);
        %end
        % plot
        cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
        rtall = mean(rt,3,'omitnan');
        rtV1V2 = nan(size(cp1));
        rtV1 = nan(size(cp1));
        rtV2 = nan(size(cp1));
        for v3i = 1:numel(c3)
            for tri = 1:numel(triggert)
                rtV1V2(v3i,tri) = mean(rt(v3i,tri,choice(v3i,tri,:) == 1 | choice(v3i,tri,:) == 2));
                rtV1(v3i,tri) = mean(rt(v3i,tri,choice(v3i,tri,:) == 1));
                rtV2(v3i,tri) = mean(rt(v3i,tri,choice(v3i,tri,:) == 2));
            end
        end
        %         subplot(3,3,plti-9*floor(plti/9 - .001)); hold on;
        %         mycol = colormap(jet(length(c3)));
        %         mycol2 = [colormap(winter(floor(length(triggert)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(triggert)/2))))];
        slope = nan(2,numel(triggert));
        for ti = 1:numel(triggert)
            %             plot(cp.cp3(:,ti),cndratio(:,ti),'-','color',mycol2(ti,:),'LineWidth',lwd/2);
            %             for i = 1:length(c3)
            %                 plot(cp.cp3(i,ti), cndratio(i,ti),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
            %             end
            mdl = fitlm(cp3(c3<=c2,ti),cndratio(c3<=c2,ti),'linear');
            slope(1,ti) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp3(c3>=c2,ti),cndratio(c3>=c2,ti),'linear');
            slope(2,ti) = mdl.Coefficients.Estimate(2);
        end
        %         xlabel('V_3');
        %         ylabel({'Conditional prob.';'Choosing Opt. 1 vs. Opt. 2'});
        %         title(sprintf('a %1.1f, b %1.1f',avec(ai),bvec(bi)),'FontSize',fontsize-5,'FontWeight','normal');
        %         savefigs(h, figname, outdir, fontsize, aspect);
        cndratiomat(ai,bi,:,:) = cndratio;
        rtallmat(ai,bi,:,:) = rtall;
        rtV1V2mat(ai,bi,:,:) = rtV1V2;
        rtV1mat(ai,bi,:,:) = rtV1;
        rtV2mat(ai,bi,:,:) = rtV2;
        slopemat(ai,bi,:,:) = slope;
    end
end
save(fullfile(outdir,[figname, '.mat']),'cndratiomat','rtallmat','rtV1V2mat','rtV1mat','rtV2mat','slopemat','avec','bvec','cp1','cp2','cp3','triggert','triggertmat');