%% behavior at two levels in the space of alpha and beta
addpath('../');
addpath(genpath('../../../CoreFunctions'));
numNode = 1;
% [sortNum, myCluster] = RndCtrl(numNode);
% mypool = parpool(myCluster, myCluster.NumWorkers);
%% behavior at two levels in the space of alpha and beta
outdir = fullfile('./graphics','LcD_GABA');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end
sims = 10000;
sgmInput = 0; % 1/3
dt = .001;
scale = 1;
sgm = 2; %6.8/2;
w = ones(2);
GABAtwo = [1.4, 1];
avec = 0:3:70; %[0:.5:70];
bvec = 0:.3:8;% [0:.05:8];
Tau = [.1, .1, .1]*2;
initialvals = [35,35; 70,70; 0,0];
presentt = 0;
triggert = 0;
dur = 6;
stimdur = dur;
thresh = 70;
stoprule = 1;
filename = sprintf('RNM_DNM_LD_dff_GABA%.1f_%1.1fa%1.0f_%2.0fb%1.0f_%2.0fsgm%1.1fsinpt%2.0f_mtx%i_sims%i',GABAtwo,min(avec),max(avec),min(bvec),max(bvec),sgm,sgmInput,numel(avec)*numel(bvec), sims);
output = fullfile('./SimRslts',[filename, '.mat']);
cp = [0 32 64 128 256 512]'/1000;
Vinput = 256*[1+cp, 1-cp];
%parpool(3);
if ~exist(output,'file')
    %     rtmat = [];
    %     choicemat = [];
    %     argmaxRmat = [];
    %     dRmat = [];
    ACC = [];
    meanRT = [];
    meandR = [];
    for gi = 1:length(GABAtwo)
        Gaba = GABAtwo(gi)
        for vi = 3%1:length(cp)
            Vmat = Vinput(vi,:);
            fprintf('cp = %1.2f',cp(vi));
            %for bi = 1:length(bvec)
                fprintf('.');
                [a,b] = meshgrid(avec,bvec);
                a = a';
                b = b';
                [rt, choice, ~, dR] = LcDsInhbt_dff_RndInput_Gaba_2mtrx_GPU(Vmat*scale, Gaba, w,...
                    a, b, sgm, sgmInput, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
                ACC(vi,:,:,gi) = gather(mean(2-squeeze(choice),3,'omitnan'));
                meanRT(vi,:,:,gi) = gather(mean(squeeze(rt),3,'omitnan'));
                meandR(vi,:,:,gi) = gather(mean(squeeze(dR),3,'omitnan')); % ,'omitnan'
                %                 rtmat(vi,:,bi,gi) = gather(squeeze(rt))';
                %                 choicemat(vi,:,bi,gi) = gather(squeeze(choice))';
                %                 argmaxRmat(vi,:,bi,gi) = gather(squeeze(argmaxR))';
                %                 dRmat(vi,:,bi,gi) = gather(squeeze(dR))';
            %end
            fprintf('\n');
        end
    end
    save(output,'ACC','meanRT','meandR','avec','bvec','GABAtwo','Vinput');
else
    load(output);
end

vi = 3;
x = avec;
y = bvec;
namex = '\alpha';
namey = '\beta';
H = Vslz_ParamSpc(x,y,namex,namey,filename,ACC,meanRT,meandR,thresh,outdir,vi);
