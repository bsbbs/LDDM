%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = fullfile('rslts','LDDM_ABS');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
% initalize parameters as in the paper of Wong & Wang, 2006
task = 'RTABS';
% simulation parameters
dt = .001; % second
presentt = dt;
triggert = presentt;
w = ones(3);
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];
sgm = 5;
dur = 4; % second
stimdur = dur;
thresh = 70; % Hz
initialvals = [1,1,1; 3,3,3; 0, 0, 0]; % for R, G, and I variables, will be further scaled in the LDDM_GPU3
stoprule = 1;
sims = 10240*2; % number of iterations
% plot conditional ratio as a function of V3

bvec = [0.6:.5:5];
eqlbvec = [5:6:70]; % 11 e
alen = 9;
amat = nan([9, numel(eqlbvec)]); % 9a * 11 e
for ei = 1:numel(eqlbvec)
    maxa = 3*eqlbvec(ei) + 1;
    afit = (3*4.7*eqlbvec(ei) + 4.7)/5.7; % based on experience, when alpha*R* and input scale in a ratio of 4.7, the model works the best
    amat(:,ei) = [flip(afit - logspace(0,log10(afit-1),ceil(alen/2))) afit + logspace(0,log10(maxa-afit-1),floor(alen/2))];
end

c = [.128]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = [0:.2:2]'; % 11 c3
cp3 = repmat(c3,1,numel(eqlbvec)); % 11 c3 * 11 e 
cp1 = repmat(c1, size(cp3));
cp2 = repmat(c2, size(cp3));
eqlbmat = repmat(eqlbvec, numel(c3), 1); % 11 c3 * 11 e

simname = sprintf('LDDM_%s_%ic1_%ic2_%ic3_%ialevel_%iblevel_w%1.1f_%ieqlbvals_sgm%2.1f_sim%i',...
            task,length(c1),length(c2),length(unique(c3)),alen,numel(bvec),w(1,1),numel(eqlbvec),sgm,sims);
%
aspect = [3,3]*ceil(sqrt(alen));
slope = [];
plti = 0;
panel = 0;
for bi = 1:length(bvec)
    for ai = 1:alen
        plti = plti + 1;
        if mod(plti,alen) == 1
            panel = panel + 1;
            h = figure;
            figname = sprintf('%s_ap%i',simname,panel);
        end
        a = repmat(amat(ai,:), numel(c3), 1); % 11 c3 * 11 e 
        b = eye(3)*bvec(bi);
        cp.cp1 = cp1;
        cp.cp2 = cp2;
        cp.cp3 = cp3;
        % simulation
        filename = sprintf('%s_alevel%i_b%1.2f',...
            simname,ai,b(1,1));
        simrslt = fullfile(outdir,[filename, '.mat']);
        if ~exist(simrslt,'file')
            [choice, rt] = LDDM_GPU3ABS(cp, eqlbmat, w, a, b, sgm, Tau, dur,...
                dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            save(simrslt,'choice','rt');
        else
            load(simrslt);
        end
        % plot
        cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
        subplot(ceil(sqrt(alen)),ceil(sqrt(alen)),plti-alen*floor(plti/alen - .001)); hold on;
        mycol = colormap(jet(length(c3)));
        mycol2 = [colormap(winter(floor(length(eqlbvec)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(eqlbvec)/2))))];
        lgd3 = [];
        for ei = 1:length(eqlbvec)
            lgd3(ei) = plot(cp.cp3(:,ei),cndratio(:,ei),'-','color',mycol2(ei,:),'LineWidth',lwd/2);
            for i = 1:length(c3)
                plot(cp.cp3(i,ei), cndratio(i,ei),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
            end
            mdl = fitlm(cp.cp3(c3<=c2,ei),cndratio(c3<=c2,ei),'linear');
            slope(1,ai,ei,bi) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp.cp3(c3>=c2,ei),cndratio(c3>=c2,ei),'linear');
            slope(2,ai,ei,bi) = mdl.Coefficients.Estimate(2);
        end
        xlabel('V_3');
        ylabel({'Conditional prob.';'Choosing Opt. 1 vs. Opt. 2'});
        title(sprintf('a level %i, b %1.1f',ai,bvec(bi)),'FontSize',fontsize-5,'FontWeight','normal');
        if mod(plti,alen) == 1
            lgd = legend(lgd3,cellstr(num2str(eqlbvec')),...
                'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
                'FontAngle','italic','NumColumns',1,'Box','off');
            title(lgd,'Eqlbvals');
        end
        savefigs(h, figname, outdir, fontsize, aspect);
    end
end
save(fullfile(outdir, [simname, '.mat']), 'slope','amat','bvec','c1','c2','c3','eqlbvec','eqlbmat');

%% to see ACC - V3 conditional on RT (controled by beta)
meanrt = [];
cndratio = [];
mycol = colormap(jet(length(c3)));
mycol3 = [colormap(winter(floor(length(bvec)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(bvec)/2))))];
h = figure; hold on;
for bi = 1:length(bvec)
    b = eye(3)*bvec(bi);
    for ai = 6
        filename = sprintf('%s_alevel%i_b%1.2f',...
            simname,ai,b(1,1));
        simrslt = fullfile(outdir,[filename, '.mat']);
        load(simrslt);
        for i = 1:6 %length(c3)
            meanrt(i,ai,bi) = mean(rt(i,8,choice(i,8,:) == 1 | choice(i,8,:) == 2),3);
            cndratio(i,ai,bi) = sum(choice(i,8,:) == 1)./(sum(choice(i,8,:) == 1) + sum(choice(i,8,:) == 2, 3));
            plot(meanrt(i,ai,bi), cndratio(i,ai,bi),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
        end
        lgd3(bi) = plot(meanrt(:,ai,bi),cndratio(:,ai,bi),'-','color',mycol3(bi,:),'LineWidth',lwd/2);
    end
end
set(gca,'XScale','log');



%% plotting slope
load(fullfile(outdir, [simname, '.mat']), 'slope','amat','bvec','c1','c2','c3','eqlbvec','eqlbmat');
[bmat, tmat, amat] = meshgrid(bvec, triggert - presentt - stimdur,avec);
[bmat, eqlbmat3] = meshgrid(bvec, eqlbmat);
[~, amat3] = meshgrid(bvec, amat);
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