addpath('../RecurrentModel/Analysis/CoreFunctions'); % inlcude the core functions
outdir = './graphics';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
mkdir('./SimRslts');
datadir = '../RecurrentModel/Figs4Paper/SimRslts';
%% define parameters for simulation
% c = [0 3.2, 6.4, 12.8, 25.6, 51.2, 76.8]'/100; % percentage of coherence
c = [3.2 12.8, 25.6, 38.4 51.2]'/100; % percentage of coherence
VmatDiag = 30*[1+c, 1-c];
Vmat = [linspace(0,60,5)', ones(5,1)*30];
V1Iter = linspace(0,60,50);
V2Iter = 60 - linspace(0,60,50);
V1IterN = linspace(0,60,100);
V2IterN = 60 - V1IterN;
V1Iterp = linspace(0,512,129);
V2Iterp = linspace(0,512,129);
dt = .001;
Tau = ones(1,3)*.1;
lwd = 2.0;
mksz = 18;
aspect1 = [3.9,2.2]; % 16:9, for wide temporal dynamic
aspect2 = [3 3]; % for temporal dynamic
aspect3 = [2.8 2.54]; % 1:1 for phase plane
aspect4 = [3.302 2.54]; % 1:1 for heatmap with colorbar
aspect5 = [4.1 2.54]; % 1:1 for phase plane with legend outside
aspect6 = [2.2 4.0]; % for picked FR in fig5e
aspect7 = aspect3*.95; % for WTA surf
aspect8 = [2, 6.4]; % for the long format RT distribution fitting panels
aspect9 = [4.6 4.0]; % for fitted time course
aspect10 = [2.8 2.1]; % for fitted acc and RT
aspect11 = [2.9,.3]; % for timeline Fig7
aspect12 = [3.9,.6]; % for timeline Fig8
aspect13 = [3, 2.0]; % for dynamic Fig7
aspect14 = [2.41 3]; % for choice and RT panel
aspect15 = [3, 11]; % for combined Fig7
fontsize = 14;
mygray = flip(gray(length(c) + 1));
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
colorpalettergb = [239,71,111;255,209,102;6,214,160;17,138,178;7,59,76]/255;
colorpalettergb =[
    0.9373    0.2784    0.4353
    1.0000    0.8196    0.4000
    0.0235    0.8392    0.6275
    0.0667    0.5412    0.6980
    0.0275    0.2314    0.2980];
%% choice accuracy and reaction time, change b
apoint = 10;
a = eye(2)*apoint;
bpoints = [.88, 1.32];
% b = 1.02*eye(2);
w = ones(2);
initialvals = [4,4;8,8;0,0];
dur = 4;
presentt = dt;
stimdur = dur;
triggert = dt;
thresh = 70;%25; %70;
stoprule = 1;
sgm = 5;
sims = 10000;
output = fullfile(datadir,sprintf('DCM_Sim%iV%ia%2.1fb%1.2f_%1.2f_sgm%1.1fChoiceRT.mat',sims,length(V1IterN),a(1,1),bpoints,sgm));
if ~exist(output,'file')
    CRALL = []; RTALL = [];
    for bi = 1:length(bpoints)
        b = bpoints(bi)*eye(2);
        RT = [];Choice = [];
        Vinput = [V1IterN', V2IterN']; % (51:80)
        [rt, choice, ~] = LcDsInhbt_GPU(Vinput, w, a, b, sgm, Tau, dur,...
            dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        Choice = squeeze(gather(choice));
        RT = squeeze(gather(rt));
        CRALL(bi,:) = mean(Choice,2,'omitnan') - 1;
        RTALL(bi,:) = mean(RT,2,'omitnan');
    end
    save(output,'CRALL','RTALL');
else
    load(output);
end
nticks = size(CRALL,2);
h = figure;
filename = 'Fig8b';
subplot(2,1,1); hold on;
lgdtext = [];
for bi = 1:length(bpoints)
    plot(1:nticks,(1-CRALL(bi,:))*100,'LineWidth',lwd,'Color',colorpalette{9-bi*4});
    lgdtext{bi} = sprintf('\\color[rgb]{%s}\\beta = %1.1f',num2str(colorpalettergb(9-bi*4,:)),bpoints(bi));
end
plot([0,50.5],[50,50],'k--');
plot([50.5,50.5],[0,50],'k--');
xlim([1,nticks]);
xticks(linspace(1,nticks,5));
xticklabels({'0','.25','.5','.75','1.0'});
xlabel(' ');
ylim([-2,100]);
yticks([0:25:100]);
ylabel('Choice (%)');
% legend(lgdtext,'Location','NorthWest','FontSize',fontsize-5,'Box','off');
savefig(h, filename, outdir, fontsize, aspect5);
subplot(2,1,2); hold on;
for bi = 1:length(bpoints)
    plot(RTALL(bi,:),'LineWidth',lwd,'Color',colorpalette{9-bi*4});
end
plot([50.5,50.5],[0,max(RTALL(:))],'k--');
xlim([1,nticks]);
xticks(linspace(1,nticks,5));
xticklabels({'0','.25','.5','.75','1.0'});
xlabel('Input ratio'); ylabel('RT (s)');
ylim([min(RTALL(:))*.8,max(RTALL(:))*1.1]);
yticks([0:.5:max(RTALL(:))]);
% ylim([.6,1.4]);
% yticks([.6:.2:1.4]);
savefig(h, filename, outdir, fontsize, aspect14);

%% example dynamics of R1 and R2
prepresentt = .12;
predur = .8;
prestimdur = predur - prepresentt;
presentt = dt;
dur = 1.5;
stimdur = dur-presentt;
triggert = presentt;
sgm = 0;
thresh = 25; %45;
stoprule = 1;
V0 = [30,30];
h = figure; hold on;
filename = 'Fig8a';
for bi = 1:2
    vi = 2;
    b = bpoints(bi)*eye(2);
    initialvals = zeros(3,2);
    [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
    Vinput = VmatDiag(vi,:);
    initialvals = [R0(end,:);G0(end,:);I0(end,:)];
    [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh*3, initialvals, stimdur, stoprule);
    %R0(end-100:end,:) = NaN;
    R = [R0; R];
    lgd2(bi) = plot(R(:,2), 'k--', 'Color', colorpalette{9-bi*4}, 'LineWidth',lwd);
    lgd1(bi) = plot(R(R(:,1)<=thresh,1), 'k-', 'Color', colorpalette{9-bi*4}, 'LineWidth',lwd);
    lgdtext{bi} = sprintf('\\color[rgb]{%s}\\beta = %1.1f',num2str(colorpalettergb(9-bi*4,:)),bpoints(bi));
end
plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-'); % threshold
text(1200,thresh*1.1,'threshold');

% plot([1, prepresentt/dt prepresentt/dt, (prepresentt+prestimdur)/dt],...
%     [0,0,2,2]-7,'k-','LineWidth',lwd/1.5);
% plot([(predur)/dt, (predur+presentt+dur)/dt],...
%     [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs (target & stimuli)
plot([1, prepresentt/dt prepresentt/dt],...
    [0,0,2]-7,'k-','LineWidth',lwd/1.5);
plot([prepresentt/dt, (prepresentt+prestimdur)/dt],...
    [2,2]-7,'k--','LineWidth',lwd/1.5); % input, target, pre-motion
plot([(predur)/dt, (predur+presentt+dur)/dt],...
    [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion

% plot([1, (prepresentt)/dt, (prepresentt)/dt, (predur)/dt, (predur)/dt, (predur+presentt+dur)/dt],[2,2,2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-11.5,'k-','LineWidth',lwd/1.5); % disinhibition
ylim([-16.5,27]);
yticks([0]);
yticklabels({'0'});
ylabel('             Activity (a.u.)');
xticks([]);
drawaxis(gca, 'x', 0, 'movelabel', 1);
xlim([1, (prepresentt+predur+presentt+dur)/dt]);
% xlabel(' ');
%
% ylim([-1,max(ylim)*1.1]);
% yticks([0]);
% xticks([prepresentt/dt, (prepresentt+prestimdur)/dt]);
% xticklabels({'Targets','Motion'});
% yticklabels({'0'});
% ylabel('Activity (a.u.)');
%xlabel('Time (a.u.)');
legend(lgd1,lgdtext,...
    'Location','NorthWest','FontSize',fontsize-5, 'FontName','Times New Roman', 'NumColumns',1,'Box','off');
savefig(h, filename, outdir, fontsize, aspect3);

%% choice accuracy and reaction time, change alpha
% avec = [5, 20];
% b = 1.02*eye(2);
% w = ones(2);
% initialvals = [4,4;8,8;0,0];
% dur = 4;
% presentt = dt;
% stimdur = dur;
% triggert = dt;
% thresh = 25; %70;
% stoprule = 1;
% sgm = 5;
% sims = 10000;
% output = fullfile(datadir,sprintf('DCM_Sim%iV%ia%2.1f_%2.1f_b%1.2f_sgm%1.1fChoiceRT.mat',sims,length(V1IterN),avec,b(1,1),sgm));
% if ~exist(output,'file')
%     CRALL = []; RTALL = [];
%     for ai = 1:length(bvec)
%         a = avec(ai)*eye(2);
%         RT = [];Choice = [];
%         Vinput = [V1IterN', V2IterN'];
%         [rt, choice, ~] = LcDsInhbt_GPU(Vinput, w, a, b, sgm, Tau, dur,...
%             dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
%         Choice = squeeze(gather(choice));
%         RT = squeeze(gather(rt));
%         CRALL(ai,:) = mean(Choice,2,'omitnan') - 1;
%         RTALL(ai,:) = mean(RT,2,'omitnan');
%     end
%     save(output,'CRALL','RTALL');
% else
%     load(output);
% end
% nticks = size(CRALL,2);
% h = figure;
% filename = 'FigSAT_Change_a';
% subplot(2,1,1); hold on;
% for ai = 1:length(avec)
%     plot(1:nticks,1-CRALL(ai,:),'LineWidth',lwd,'Color',colorpalette{ai});
%     lgdtext{ai} = sprintf('\\color[rgb]{%s}\\alpha = %1.1f',num2str(colorpalettergb(ai,:)),avec(ai));
% end
% plot([0,50.5],[.5,.5],'k--');
% plot([50.5,50.5],[0,.5],'k--');
% xlim([1,nticks]);
% xticks(linspace(1,nticks,5));
% xticklabels({'0','.25','.5','.75','1.0'});
% xlabel(' ');
% ylim([-.02,1]);
% yticks([0:.25:1]);
% ylabel('Choice (%)');
% %legend(lgdtext,'Location','NorthWest','FontSize',fontsize-5,'Box','off');
% savefig(h, filename, outdir, fontsize, aspect5);
% subplot(2,1,2); hold on;
% for ai = 1:length(avec)
%     plot(RTALL(ai,:),'LineWidth',lwd,'Color',colorpalette{ai});
% end
% plot([50.5,50.5],[0,max(RTALL(:))],'k--');
% xlim([1,nticks]);
% xticks(linspace(1,nticks,5));
% xticklabels({'0','.25','.5','.75','1.0'});
% xlabel('Input ratio'); ylabel('RT (a.u.)');
% ylim([min(RTALL(:))*.8,max(RTALL(:))*1.1]);
% yticks([.2:.2:max(RTALL(:))]);
% % ylim([.6,1.4]);
% % yticks([.6:.2:1.4]);
% savefig(h, filename, outdir, fontsize, aspect14);

%% The full space of alpha and beta
sims = 10000;
sgmInput = 0; % 1/3
dt = .001;
scale = 1;
sgm = 5; %6.8/2;
w = 1;
v = 1;
Gaba = 1;
initialvals = [4,4;8,8;0,0];
dur = 8;
presentt = dt;
stimdur = dur;
triggert = dt;
thresh = 70;
stoprule = 1;
V = VmatDiag(2,:);
%avec = linspace(0,60,601);
avec = 10.^[-1:.01:3];
rvec = linspace(0,4,401); % ratio of beta/w
Na = numel(avec);
Nb = numel(rvec);
bvec = rvec*w(1,1);
[a,b] = meshgrid(avec,bvec);
a = a';
b = b';
filename = sprintf('RTACC_fullspace_V%1.0f_%1.0f_w%1.1f_v%1.1f_a%i_%1.2f_%1.2fr%i_%1.2f_%1.2f_thresh%i',...
    [V, w, v, length(avec),min(avec),max(avec),length(rvec),min(rvec),max(rvec),thresh]);
output = fullfile(datadir,[filename, '.mat']);
if ~exist(output, 'file')
    RT = [];
    Choice = [];
    psims = 50;
    Npiles = ceil(sims/psims);
    for pile = 1:Npiles
        fprintf('.');
        [rt, choice, ~, ~] = LcDsInhbt_RndInput_Gaba_2mtrx_GPU(V*scale, Gaba, [w,v;v,w],...
            a, b, sgm, sgmInput, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, psims);
        %         [rt, choice, ~, ~] = LcDsInhbt_dff_RndInput_Gaba_2mtrx_GPU(V*scale, Gaba, w,...
        %             a, b, sgm, sgmInput, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, psims);
        Choice(:,:,(pile-1)*psims+(1:psims)) = gather(2-squeeze(choice));
        RT(:,:,(pile-1)*psims+(1:psims)) = gather(squeeze(rt));
        %meandR(:,:,(pile-1)*psims+(1:psims)) = gather(squeeze(dR));
    end
    save(output,'Choice','RT','avec','rvec','V','w','sgm','sgmInput','-v7.3');
else
    load(output);
end
%% get an clear-cut boundary from analytical analysis
%% panel e, parameter space for choice/representation under equal inputs
V = VmatDiag(2,:);
w = 1;
v = 1;
V1 = V(1);
V2 = V(2);
%avec = linspace(0,60,601);
avec = 10.^[-1:.01:3];
rvec = linspace(0,4,401); % ratio of beta/w
Nb = length(rvec);
filename = sprintf('Stability_Numeric_V%1.0f_%1.0f_w%1.1f_v%1.1f_a%i_%1.2f_%1.2fr%i_%1.2f_%1.2f',...
    [V, w, v, length(avec),min(avec),max(avec),length(rvec),min(rvec),max(rvec)]);
output = fullfile(datadir,[filename, '.mat']);
if ~exist(output, 'file')
    Npoints = NaN(length(avec),length(rvec));
    Eigenvalue = -ones(length(avec),length(rvec),4);
    Stability = NaN(length(avec),length(rvec));
    for ai = 1:length(avec)
        a = avec(ai);
        fprintf('%2.1f.',a);
        syms R1 R2;
        vars = [R1 R2];
        parfor bi = 1:Nb
            fprintf('.');
            b = rvec(bi)*w;
            eqns = [(V1/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
                (V2/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
            [AnswR1,AnswR2] = solve(eqns, vars);
            AnswI1 = b*AnswR1;
            AnswI2 = b*AnswR2;
            AnswG1 = w*AnswR1 + v*AnswR2 - AnswI1;
            AnswG2 = v*AnswR1 + w*AnswR2 - AnswI2;
            PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
            Npoints(ai,bi) = sum(PositiveRealAnsw);
            NegEigen = -ones(1,4);
            for i = find(PositiveRealAnsw)'
                R1star = double(AnswR1(i));
                R2star = double(AnswR2(i));
                G1star = double(AnswG1(i));
                G2star = double(AnswG2(i));
                JMat = [-1 + a/(1+G1star), -(V1+a*R1star)/(1+G1star)^2, 0, 0, 0, 0
                    w, -1, -1, v, 0, 0
                    b, 0, -1, 0, 0, 0
                    0, 0, 0, -1+a/(1+G2star),   -(V2+a*R2star)/(1+G2star)^2, 0
                    v, 0, 0, w, -1, -1
                    0, 0, 0, b, 0, -1];
                A = eig(JMat);
                NegEigen(i) = all(real(A) < 0); % attractive 1, diversive 0
            end
            Stability(ai,bi)  = all(NegEigen>0);
            Eigenvalue(ai,bi,:) = NegEigen;
        end
        fprintf('\n');
    end
    visualize = sum(Eigenvalue,3);
    visualize(visualize == -2) = 2;
    save(output,'avec','rvec','V','w','v','Npoints','Eigenvalue','Stability','visualize');
else
    load(output);
end
% plot
dyadic = visualize == 2 | visualize == 0;
h = figure;colormap(colorpalettergb([3,5],:));
% filename = 'RTACC_fullspace_Stability_Numeric';
s = surf(dyadic+1,'EdgeColor','none');
ylim([1,length(avec)]);
xlim([1,length(rvec)]);
xticks(linspace(1,length(rvec),5));
xticklabels(linspace(0,4,5));
yticks(linspace(1,length(avec),5));
% yticklabels(linspace(0,60,5));
yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel('\beta/\omega');
ylabel('\alpha');
view(0,90);
savefig(h, filename, outdir, fontsize, aspect7);

%% plot - meanRT and ACC
mask = NaN(size(dyadic));
mask(~dyadic) = 1;
xmat = a;
ymat = b;
y = avec;
x = bvec;
ACC = mean(Choice,3,'omitnan'); % ,'omitnan'
h=figure;   colormap(jet);
% filename = sprintf('RTACC_fullspace_V%1.0f_%1.0f_w%1.1f_v%1.1f_a%i_%1.2f_%1.2fr%i_%1.2f_%1.2f',...
%     [V, w, v, length(avec),min(avec),max(avec),length(rvec),min(rvec),max(rvec)]);
filename = 'Fig8cd';
subplot(1,2,1); hold on;
s = surf(ACC*100); % .*mask
s.EdgeColor = 'none';
for bi = 1:numel(bpoints)
    x = find(bvec == bpoints(bi));
    y = find(avec == apoint);
    plot3(x,y,100,'.','MarkerSize',mksz,'Color',colorpalette{9-bi*4});
end
ylim([81,length(avec)]);
xlim([1,length(rvec)]);
xticks(linspace(1,length(rvec),5));
xticklabels(linspace(0,4,5));
yticks(linspace(1,find(avec==1000),5));
yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel('\beta');
ylabel('\alpha');
view(0,90);
c = colorbar;
ylabel(c, 'Accuracy [%]');
savefig(h, filename, outdir, fontsize, aspect4.*[2.4,1]);

meanRT = mean(RT,3, 'omitnan'); % , 'omitnan'
subplot(1,2,2); hold on;
s = surf(meanRT); %.*mask
s.EdgeColor = 'none';

for bi = 1:numel(bpoints)
    x = find(bvec == bpoints(bi));
    y = find(avec == apoint);
    plot3(x,y,10,'.','MarkerSize',mksz,'Color',colorpalette{9-bi*4});
end

ylim([81,length(avec)]);
xlim([1,length(rvec)]);
xticks(linspace(1,length(rvec),5));
xticklabels(linspace(0,4,5));
yticks(linspace(1,find(avec==1000),5));
yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel('\beta');
ylabel('\alpha');
view(0,90);
c = colorbar;
set(gca, 'ColorScale','log');
ylabel(c, 'RT [s]');
savefig(h, filename, outdir, fontsize, aspect4.*[2.4,1]);
%% See change of ACC as a function of RT
ACC = mean(Choice,3,'omitnan'); % ,'omitnan'
meanRT = mean(RT,3, 'omitnan'); % , 'omitnan'
h = figure; hold on;
filename = 'WrkingMem';
plot(meanRT(90,:),ACC(90,:)*100,'r.'); % low alpha
plot(meanRT(160,:),ACC(160,:)*100,'b.'); % high alpha
legend({'low \alpha','high \alpha'},'Location','Best');
xlabel('RT(secs)');
ylabel('Accuracy rate (%)');
savefig(h, filename, outdir, fontsize, aspect4);


%% NMDA Antagonism
cp = 10.^[-2:.1:-.5]';%[5 32:32:512]'/1000;
sims = 20000;
dt = .001;
sgm = 5; %6.8/2;
NMDAlevel = [1, .6];
w = ones(2);
a = eye(2)*30;
b = eye(2)*1.02;
Tau = [.1, .1, .1];
initialvals = [4,4;8,8;0,0];
presentt = dt;
triggert = presentt;
dur = 5;
stimdur = dur;
thresh = 70;
stoprule = 1;
output = fullfile('./SimRslts',sprintf('LDDM_NMDA_a%2.0f_b%.2f_sgm%1.1f_%i.mat',a(1,1),b(1,1),sgm,sims));
Vinput = 256*[1+cp, 1-cp];
if ~exist(output,'file')
    ACC = [];
    meanRT = [];
    for Ni = 1:numel(NMDAlevel)
        NMDA = NMDAlevel(Ni);
        [rt, choice, ~] = LcDsInhbt_NMDA_GPU(Vinput, NMDA, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC(:,Ni) = gather(mean(2-squeeze(choice),2,'omitnan'));
        meanRT(:,Ni) = gather(mean(squeeze(rt),2,'omitnan')); %,'omitnan'
    end
    save(output,'ACC','meanRT','Vinput','a','b','w','Tau','initialvals','NMDAlevel');
else
    load(output);
end
h = figure; hold on;
filename = 'NMDA';
lg = [];
for Ni = 1:numel(NMDAlevel)
    subplot(2,1,1); hold on;
    lg(Ni) = plot(cp,meanRT(:,Ni),'.-');
    set(gca, 'XScale','log');
    ylabel('RT (secs)');
    subplot(2,1,2); hold on;
    plot(cp,ACC(:,Ni),'.-');
    set(gca, 'XScale','log');
    ylabel('Accuracy');
    xlabel('Coherence');
end
subplot(2,1,1); hold on;
legend(lg,{'control','Antagonism'},'Location','southwest');
savefig(h, filename, outdir, fontsize-5, aspect14);
subplot(2,1,2);
savefig(h, filename, outdir, fontsize-5, aspect14);
