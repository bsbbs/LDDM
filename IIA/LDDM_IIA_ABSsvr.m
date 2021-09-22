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
outdir = fullfile('rslts','LDDM_ABS');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
%% initalize parameters as in the paper of Wong & Wang, 2006
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
dur = 7; % second
stimdur = dur;
thresh = 70; % Hz
initialset = [1,1,1; 3,3,3; 0, 0, 0]; % for R, G, and I variables
stoprule = 1;
sims = 10240*2; % number of iterations
% plot conditional ratio as a function of V3

%%
avec = 0:1:95; %[10,30,50,60,70,75,80,85,90]; % 0,15,20,96
bvec = [0.6:.2:4.4];
eqlbvec = [5:5:70];
c = .128;
c1 = 1 + c;
c2 = ones(size(c1));
c3 = 0:.1:2;
[sclmat, cp3] = meshgrid(eqlbvec, c3);
cp1 = repmat(c1, size(sclmat));
cp2 = repmat(c2, size(sclmat));
initialvals = initialset.*sclmat;
figname = sprintf('IIA_LDDM_%s_stimdur%1.1f_sgm%2.1f_%iA%iB%iT',task,stimdur,sgm,numel(avec),numel(bvec),numel(eqlbvec));

cndratiomat = nan([numel(avec),numel(bvec),size(cp1)]);
rtallmat = nan([numel(avec),numel(bvec),size(cp1)]);
rtV1V2mat = nan([numel(avec),numel(bvec),size(cp1)]);
slopemat = nan([numel(avec),numel(bvec),2,numel(eqlbvec)]);
for ai = 1:length(avec)
    scale = 3*sclmat.^2 + (1-avec(ai)).*sclmat;
    Vinput.V1 = cp1.*scale;
    Vinput.V2 = cp2.*scale;
    Vinput.V3 = cp3.*scale;
    a = eye(3)*avec(ai);
    parfor bi = 1:length(bvec)
        b = bvec(bi);
        
        [choice, rt] = LDDM_GPU3(Vinput, w, a, eye(3)*b, sgm, Tau, dur,...
            dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        cndratio = sum(choice == 1, 3)./(sum(choice == 1, 3) + sum(choice == 2, 3)); % conditional ratio of choosing 1 vs. choosing 2
        rtall = mean(rt,3,'omitnan');
        rtV1V2 = nan(size(cp1));
        rtV1 = nan(size(cp1));
        rtV2 = nan(size(cp1));
        for v3i = 1:numel(c3)
            for scli = 1:numel(eqlbvec)
                rtV1V2(v3i,scli) = mean(rt(v3i,scli,choice(v3i,scli,:) == 1 | choice(v3i,scli,:) == 2));
                rtV1(v3i,scli) = mean(rt(v3i,scli,choice(v3i,scli,:) == 1));
                rtV2(v3i,scli) = mean(rt(v3i,scli,choice(v3i,scli,:) == 2));
            end
        end
        slope = nan(2,numel(eqlbvec));
        for scli = 1:numel(eqlbvec)
            mdl = fitlm(cp3(c3<=c2,scli),cndratio(c3<=c2,scli),'linear');
            slope(1,scli) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cp3(c3>=c2,scli),cndratio(c3>=c2,scli),'linear');
            slope(2,scli) = mdl.Coefficients.Estimate(2);
        end
        cndratiomat(ai,bi,:,:) = cndratio;
        rtallmat(ai,bi,:,:) = rtall;
        rtV1V2mat(ai,bi,:,:) = rtV1V2;
        rtV1mat(ai,bi,:,:) = rtV1;
        rtV2mat(ai,bi,:,:) = rtV2;
        slopemat(ai,bi,:,:) = slope;
    end
end
save(fullfile(outdir,[figname, '.mat']),'cndratiomat','rtallmat','rtV1V2mat','rtV1mat','rtV2mat','slopemat','avec','bvec','sclvec','cp1','cp2','cp3','sclmat');