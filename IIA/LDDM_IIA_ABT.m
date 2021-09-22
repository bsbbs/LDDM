%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
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
% initalize parameters as in the paper of Wong & Wang, 2006
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
initialvals = [4, 4, 4; 12, 12, 12; 0, 0, 0]; % for R, G, and I variables
stoprule = 1;
sims = 10240*2; % number of iterations
% plot conditional ratio as a function of V3
triggert = 2:.1:3;
c = [.256]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.2:2;
[triggertmat, cp.cp3] = meshgrid(triggert, c3);
cp.cp1 = repmat(c1, size(triggertmat));
cp.cp2 = repmat(c2, size(triggertmat));
eqlb = 6;

%
aspect = [9, 9];
slope = [];
plti = 0;
panel = 0;
avec = [10,30,50,60,70,75,80,85,90]; % 0,15,20,96
bvec = [0.6:.5:5];
for bi = 1:length(bvec)
    for ai = 1:length(avec)
        plti = plti + 1;
        if mod(plti,9) == 1
            panel = panel + 1;
            h = figure;
            figname = sprintf('IIA_LDDM_%s_eqlb%2.0f_stimdur%1.1f_ABT_ap%i',task,eqlb,stimdur,panel);
        end
        scale = 3*eqlb^2 + (1-avec(ai))*eqlb;
        Vinput.V1 = cp.cp1*scale;
        Vinput.V2 = cp.cp2*scale;
        Vinput.V3 = cp.cp3*scale;
        a = eye(3)*avec(ai);
        b = eye(3)*bvec(bi);
        % simulation
        filename = sprintf('LDDM_%s_eqlb%2.0f_%ic1_%ic2_%ic3_a%2.1f_w%1.1f_b%1.2f_scale%2.1f_%iact%1.2f-%1.2f_sgm%2.1f_sim%i',...
            task,eqlb,length(c1),length(unique(c2)),length(unique(c3)),a(1,1),w(1,1),b(1,1),scale,numel(triggert),triggert(1),triggert(2),sgm,sims);
        simrslt = fullfile(outdir,[filename, '.mat']);
        if ~exist(simrslt,'file')
            [choice, rt] = LDDM_GPU3(Vinput, w, a, b, sgm, Tau, dur,...
                dt, presentt, triggertmat, thresh, initialvals, stimdur, stoprule, sims);
            save(simrslt,'choice','rt');
        else
            load(simrslt);
        end
        % plot
        cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
        subplot(3,3,plti-9*floor(plti/9 - .001)); hold on;
        mycol = colormap(jet(length(c3)));
        mycol2 = [colormap(winter(floor(length(triggert)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(triggert)/2))))];
        for ti = 1:length(triggert)
            plot(cp.cp3(:,ti),cndratio(:,ti),'-','color',mycol2(ti,:),'LineWidth',lwd/2);
            for i = 1:length(c3)
                plot(cp.cp3(i,ti), cndratio(i,ti),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
            end
            mdl = fitlm(cp.cp3(c3<=c2,ti),cndratio(c3<=c2,ti),'linear');
            slope(1,ti,bi,ai) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp.cp3(c3>=c2,ti),cndratio(c3>=c2,ti),'linear');
            slope(2,ti,bi,ai) = mdl.Coefficients.Estimate(2);
        end
        xlabel('V_3');
        ylabel({'Conditional prob.';'Choosing Opt. 1 vs. Opt. 2'});
        title(sprintf('a %1.1f, b %1.1f',avec(ai),bvec(bi)),'FontSize',fontsize-5,'FontWeight','normal');
        savefigs(h, figname, outdir, fontsize, aspect);
    end
end
figname = sprintf('IIA_LDDM_%s_eqlb%2.0f_stimdur%1.1f_sgm%2.1f_%iA%iB%iT',task,eqlb,stimdur,sgm,numel(avec),numel(bvec),numel(triggert));
save(fullfile(outdir, [figname, '.mat']), 'slope','avec','bvec','cp','triggert');

%% plotting slope
figname = sprintf('IIA_LDDM_%s_stimdur%1.1f_sgm%2.1f_%iA%iB%iT',task,stimdur,sgm,numel(avec),numel(bvec),numel(triggert));
load(fullfile(outdir, [figname, '.mat']), 'slope','avec','bvec','cp','triggert');
eqlb = 32;
[bmat, tmat, amat] = meshgrid(bvec, triggert - presentt - stimdur,avec);
segtext = {'V3<=V2','V3>=V2'};
for seg = 1:2 % V3<=V2 or V3>=V2
    h = figure; hold on;
    figname = sprintf('IIA_LDDM_%s_eqlb%2.0f_stimdur%1.1f_sgm%2.1f_%iA%iB%iT_%s',task,eqlb,stimdur,sgm,numel(avec),numel(bvec),numel(triggert),segtext{seg});
    c = squeeze(slope(seg,:,:,:));
    scatter3(tmat(:),bmat(:),amat(:),mksz*3,c(:),'filled');
    grid on;
    colormap(bluewhitered);
    colorbar;
    xlabel('action time'); ylabel('\beta'); zlabel('\alpha');
    view([79, 26]);
    set(gca,'color',[128,128,128]/256);
    savefig(h, fullfile(outdir,[figname,'.fig']));
    savefigs(h, figname, outdir, fontsize, [8 5]);
end