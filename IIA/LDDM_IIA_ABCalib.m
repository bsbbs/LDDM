%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = fullfile('rslts','LDDM_ABC');
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
initialvals = [1,1,1; 3,3,3; 0, 0, 0]; % for R, G, and I variables, will be further scaled
stoprule = 1;
sims = 10240*2; % number of iterations
% plot conditional ratio as a function of V3
triggert = [flip(2.5 - logspace(-1,log10(2.5),5)),2.5, 2.5 + logspace(-1,log10(2.5),5)]; %11 t
c = [.256]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.2:2; % 11c
[triggertmat, cp.cp3] = meshgrid(triggert, c3); % 11 c3 * 11 t
cp.cp1 = repmat(c1, size(triggertmat));
cp.cp2 = repmat(c2, size(triggertmat));
eqlb = 47;
alen = 9;
maxa = 3*eqlb + 1;
afit = (3*4.7*eqlb + 4.7)/5.7; % based on experience, when alpha*R* and input scale in a ratio of 4.7, the model works the best
avec = [flip(afit - logspace(0,log10(afit-1),ceil(alen/2))) afit + logspace(0,log10(maxa-afit-1),floor(alen/2))];
scale = 3*eqlb.^2 + (1-avec).*eqlb;
bvec = [0.6:.5:5];
%
simname = sprintf('LDDM_%s_%ic1_%ic2_%ic3_%ialevel_%iblevel_w%1.1f_%icalibrate_sgm%2.1f_sim%i',...
            task,length(c1),length(c2),length(c3),alen,numel(bvec),w(1,1),numel(triggert),sgm,sims);
aspect = [3,3]*ceil(sqrt(alen));
slope = [];
plti = 0;
panel = 0;
cndratiomat = nan([numel(avec),numel(bvec),size(cp1)]);
rtallmat = nan([numel(avec),numel(bvec),size(cp1)]);
rtV1V2mat = nan([numel(avec),numel(bvec),size(cp1)]);
slopemat = nan([numel(avec),numel(bvec),2,numel(Dvec)]);
for bi = 1:length(bvec)
    for ai = 1:length(avec)
        plti = plti + 1;
        if mod(plti,alen) == 1
            panel = panel + 1;
            h = figure;
            figname = sprintf('%s_ap%i',simname,panel);
        end
        scale = 3*eqlb^2 + (1-avec(ai))*eqlb;
        a = avec(ai);
        b = eye(3)*bvec(bi);
        % simulation
        filename = sprintf('%s_a%2.1f_b%1.2f',...
            simname,avec(ai),b(1,1));
        simrslt = fullfile(outdir,[filename, '.mat']);
        if ~exist(simrslt,'file')
            [choice, rt] = LDDM_GPU3ABS(cp, eqlb, w, a, b, sgm, Tau, dur,...
                dt, presentt, triggertmat, thresh, initialvals, stimdur, stoprule, sims);
            save(simrslt,'choice','rt');
        else
            load(simrslt);
        end
        % plot
        cndratio = sum(choice == 2, 3)./(sum(choice == 2, 3) + sum(choice == 3, 3)); % conditional ratio of choosing 1 vs. choosing 2
        subplot(3,3,plti-9*floor(plti/9 - .001)); hold on;
        mycol = colormap(jet(length(c3)));
        mycol2 = [colormap(winter(floor(length(triggert)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(triggert)/2))))];
        lgd3 = [];
        for ti = 1:length(triggert)
            lgd3(ti) = plot(cp.cp3(:,ti),cndratio(:,ti),'-','color',mycol2(ti,:),'LineWidth',lwd/2);
            for i = 1:length(c3)
                plot(cp.cp3(i,ti), cndratio(i,ti),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
            end
            mdl = fitlm(cp.cp3(c3<c2,ti),cndratio(c3<c2,ti),'linear');
            slope(1,ti) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp.cp3(c3>c2,ti),cndratio(c3>c2,ti),'linear');
            slope(2,ti) = mdl.Coefficients.Estimate(2);
        end
        xlabel('V_3');
        ylabel({'Conditional prob.';'Choosing Opt. 1 vs. Opt. 2'});
        title(sprintf('a %1.1f, b %1.1f',avec(ai),bvec(bi)),'FontSize',fontsize-5,'FontWeight','normal');
        if mod(plti,alen) == 1
            lgd = legend(lgd3,cellstr(num2str(triggert')),...
                'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
                'FontAngle','italic','NumColumns',1,'Box','off');
            title(lgd,'Action time');
        end
        savefigs(h, figname, outdir, fontsize, aspect);
        
        rtall = mean(rt,3,'omitnan');
        rtV1V2 = nan(size(cp.cp1)); % % 11 c3 * 11 t
        rtV1 = nan(size(cp.cp1));
        rtV2 = nan(size(cp.cp1));
        for v3i = 1:numel(c3)
            for di = 1:numel(triggert)
                rtV1V2(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 2 | choice(v3i,di,:) == 3));
                rtV1(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 2));
                rtV2(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 3));
            end
        end
        cndratiomat(ai,bi,:,:) = cndratio;
        rtallmat(ai,bi,:,:) = rtall;
        rtV1V2mat(ai,bi,:,:) = rtV1V2;
        rtV1mat(ai,bi,:,:) = rtV1;
        rtV2mat(ai,bi,:,:) = rtV2;
        slopemat(ai,bi,:,:) = slope;
        
    end
end
save(fullfile(outdir, [simname, '.mat']), 'cndratiomat','rtallmat','rtV1V2mat','rtV1mat','rtV2mat','slopemat',...
    'avec','bvec','cp','triggert','triggertmat','eqlb');

%% plotting slope
load(fullfile(outdir, [simname, '.mat']), 'cndratiomat','rtallmat','rtV1V2mat','rtV1mat','rtV2mat','slopemat',...
    'avec','bvec','cp','triggert','triggertmat','eqlb');
[bmat, amat, tmat] = meshgrid(bvec,avec,triggert);
segtext = {'V3<=V2','V3>=V2'};
for seg = 1:2 % V3<=V2 or V3>=V2
    h = figure; hold on;
    figname = sprintf('%s_Slope_%s',simname,segtext{seg});
    c = squeeze(slopemat(:,:,seg,:));
    scatter3(amat(:),bmat(:),tmat(:),mksz*3,c(:),'filled');
    grid on;
    colormap(bluewhitered);
    colorbar;
    zlabel('Action time (s)'); xlabel('\alpha'); ylabel('\beta');
    view([79, 26]);
    set(gca,'color',[128,128,128]/256);
    savefig(h, fullfile(outdir,[figname,'.fig']));
    savefigs(h, figname, outdir, fontsize, [8 5]);
end