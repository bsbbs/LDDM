%% to test Wong & Wang 2006 model in predicting
% independence of irrelevant alternatives (I.I.A.)

%% including packages and dircetories
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));

%% manipulating output directories
outdir = fullfile('rslts','LDDM_ABD');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
% initalize parameters as in the paper of Wong & Wang, 2006
task = 'RTABD';
% simulation parameters
dt = .001; % second
presentt = dt;
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
eqlb = 47; 
Dvec = [logspace(-3,log10(2),11)];% 11 delays
alen = 9;
maxa = 3*eqlb + 1;
afit = (3*4.7*eqlb + 4.7)/5.7; % based on experience, when alpha*R* and input scale in a ratio of 4.7, the model works the best
avec = [flip(afit - logspace(0,log10(afit-1),ceil(alen/2))) afit + logspace(0,log10(maxa-afit-1),floor(alen/2))];
scale = 3*eqlb.^2 + (1-avec).*eqlb;
c = [.128]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = [0:.2:2]'; % 11 c3
cp3 = repmat(c3,1,numel(Dvec)); % 11 c3 * 11 d 
cp1 = repmat(c1, size(cp3));
cp2 = repmat(c2, size(cp3));
triggert = repmat(Dvec, numel(c3),1); % 11 c3 * 11 d

simname = sprintf('LDDM_%s_%ic1_%ic2_%ic3_%ialevel_%iblevel_w%1.1f_%idelays_sgm%2.1f_sim%i',...
            task,length(c1),length(c2),length(c3),alen,numel(bvec),w(1,1),numel(Dvec),sgm,sims);
%
aspect = [3,3]*ceil(sqrt(alen));
slope = [];
plti = 0;
panel = 0;
cndratiomat = nan([numel(avec),numel(bvec),size(cp1)]);
rtallmat = nan([numel(avec),numel(bvec),size(cp1)]);
rtV1V2mat = nan([numel(avec),numel(bvec),size(cp1)]);
slopemat = nan([numel(avec),numel(bvec),2,numel(Dvec)]);
for bi = 1:length(bvec)
    for ai = 1:alen
        plti = plti + 1;
        if mod(plti,alen) == 1
            panel = panel + 1;
            h = figure;
            figname = sprintf('%s_ap%i',simname,panel);
        end
        a = avec(ai);
        b = eye(3)*bvec(bi);
        cp.cp1 = cp1;
        cp.cp2 = cp2;
        cp.cp3 = cp3;
        % simulation
        filename = sprintf('%s_alevel%i_b%1.2f',...
            simname,ai,b(1,1));
        simrslt = fullfile(outdir,[filename, '.mat']);
        if ~exist(simrslt,'file')
            [choice, rt] = LDDM_GPU3ABS(cp, eqlb, w, a, b, sgm, Tau, dur,...
                dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            save(simrslt,'choice','rt');
        else
            load(simrslt);
        end
        % plot
        cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
        subplot(ceil(sqrt(alen)),ceil(sqrt(alen)),plti-alen*floor(plti/alen - .001)); hold on;
        mycol = colormap(jet(length(c3)));
        mycol2 = [colormap(winter(floor(length(Dvec)/2))); colormap(bone(1));flip(colormap(autumn(floor(length(Dvec)/2))))];
        lgd3 = [];
        for di = 1:length(Dvec)
            lgd3(di) = plot(cp.cp3(:,di),cndratio(:,di),'-','color',mycol2(di,:),'LineWidth',lwd/2);
            for i = 1:length(c3)
                plot(cp.cp3(i,di), cndratio(i,di),'.','color',mycol(i,:),'LineWidth',lwd,'MarkerSize',mksz*3);
            end
            mdl = fitlm(cp.cp3(c3<c2,di),cndratio(c3<c2,di),'linear');
            slope(1,di) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp.cp3(c3>c2,di),cndratio(c3>c2,di),'linear');
            slope(2,di) = mdl.Coefficients.Estimate(2);
        end
        xlabel('V_3');
        ylabel({'Conditional prob.';'Choosing Opt. 1 vs. Opt. 2'});
        title(sprintf('a %2.1f, b %1.1f',avec(ai),bvec(bi)),'FontSize',fontsize-5,'FontWeight','normal');
        if mod(plti,alen) == 1
            lgd = legend(lgd3,cellstr(num2str(Dvec')),...
                'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
                'FontAngle','italic','NumColumns',1,'Box','off');
            title(lgd,'beta Delays');
        end
        savefigs(h, figname, outdir, fontsize, aspect);
        
        rtall = mean(rt,3,'omitnan');
        rtV1V2 = nan(size(cp1)); % 11 c3 * 11 d 
        rtV1 = nan(size(cp1));
        rtV2 = nan(size(cp1));
        for v3i = 1:numel(c3)
            for di = 1:numel(Dvec)
                rtV1V2(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 1 | choice(v3i,di,:) == 2));
                rtV1(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 1));
                rtV2(v3i,di) = mean(rt(v3i,di,choice(v3i,di,:) == 2));
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
    'avec','bvec','cp','eqlb','c1','c2','c3','Dvec');

%% plotting slope
load(fullfile(outdir, [simname, '.mat']), 'cndratiomat','rtallmat','rtV1V2mat','rtV1mat','rtV2mat','slopemat',...
    'avec','bvec','cp','eqlb','c1','c2','c3','Dvec');
[bmat, amat, tmat] = meshgrid(bvec,avec,Dvec);
segtext = {'V3<=V2','V3>=V2'};
for seg = 1:2 % V3<=V2 or V3>=V2
    h = figure; hold on;
    figname = sprintf('%s_Slope_%s',simname,segtext{seg});
    c = squeeze(slopemat(:,:,seg,:));
    scatter3(amat(:),bmat(:),tmat(:),mksz*3,c(:),'filled');
    grid on;
    colormap(bluewhitered);
    colorbar;
    zlabel('Action time'); ylabel('\beta'); xlabel('\alpha');
    view([79, 26]);
    set(gca,'color',[128,128,128]/256);
    savefig(h, fullfile(outdir,[figname,'.fig']));
    savefigs(h, figname, outdir, fontsize, [8 5]);
end