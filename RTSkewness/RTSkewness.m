Homedir = 'C:\Users\Bo';
% Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
addpath(genpath(fullfile(Homedir,'Documents','LDDM','Fit')));
cd('G:\My Drive\LDDM\RTSkewness');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/RTSkewness');
out_dir = './';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
sim_dir = fullfile(out_dir,'SimRslts');
if ~exist(sim_dir,'dir')
    mkdir(sim_dir);
end
dataDynmc = load('../Fit/Data/Data.mat');
dataBhvr = LoadRoitmanData('../RoitmanDataCode');
randseed = 24356545;
rng(randseed);
%%
% avec = 10.^[-1:.01:1.9];
avec = 0.5:.5:80;
% bvec = linspace(0,4,401);
bvec = 0.01:.01:2;
[Amat,Bmat] = meshgrid(avec,bvec);
params = [0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132	16539.138186];
sgm = params(3);
tauR = params(5);
tauG = params(6);
tauI = params(7);

predur = 0;
presentt = 0; 
triggert = 0;
dur = 5;
dt =.001;
thresh = 70; %70.8399; % mean(max(m_mr1cD))+1; 
stimdur = dur;
stoprule = 1;
w = [1 1; 1 1];
Rstar = 42; % ~ 42 Hz at the bottom of initial fip, according to Roitman and Shadlen's data
initialvals = [Rstar,Rstar; sum(w(1,:))*Rstar,sum(w(2,:))*Rstar; 0,0];
eqlb = Rstar; % set equilibrium value before task as R^*
Tau = [tauR tauG tauI];
sims = 10240;
Npiles = 1;
Bdim = numel(bvec);
Nbs = 20;
Bsec = ceil(Bdim/Nbs);
cp = .128;

effN = [];
SK = [];
KT = [];
M = [];
SD = [];
ACC = [];
name = sprintf('Amat%i_Bmat%i_sgm%2.1f_tau%1.2f_%1.2f_%1.2f_%i',numel(avec),numel(bvec),params([3,5:7]),sims*Npiles);
if ~exist(fullfile(sim_dir,sprintf('PlotData_%s.mat',name)),'file')
    for pi = 1:Npiles
        for bi = 1:Nbs
            startp = ((bi-1)*Bsec + 1);
            endp = min((bi*Bsec), Bdim);
            bbvec = bvec(startp:endp);
            [amat,bmat] = meshgrid(avec,bbvec);
            scale = 2*1*eqlb.^2 + (1-amat).*eqlb; % a = 1-(scale/eqlb - 2*eqlb)
            Vinput.V1 = (1+cp)*scale;
            Vinput.V2 = (1-cp)*scale;
            Vprior.V1 = scale;
            Vprior.V2 = scale;
            tic;
            [rt, choice, ~] = LDDM_GPU_ABMtrx(Vprior, Vinput, w, amat, bmat,...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            toc
            rtmat(startp:endp,:,(sims*(pi-1)+1):(sims*pi)) = rt;
            choicemat(startp:endp,:,(sims*(pi-1)+1):(sims*pi)) = choice;
        end
    end

    for ai = 1:numel(avec)
        for bi = 1:numel(bvec)
            % mode: mean, deviance, skewness, turkosis
            simdata = rtmat(bi,ai,~isnan(choicemat(bi,ai,:)));
            ACC(bi,ai) = mean(2-choicemat(bi,ai,:),'omitnan');
            effN(bi,ai) = numel(simdata);
            SK(bi,ai) = skewness(simdata);
            KT(bi,ai) = kurtosis(simdata);
            M(bi,ai) = mean(simdata);
            SD(bi,ai) = std(simdata);
        end
    end

    save(fullfile(sim_dir,sprintf('PlotData_%s.mat',name)),...
                'rtmat','choicemat','-v7.3');
    save(fullfile(sim_dir,sprintf('CalculatedData_%s.mat',name)),...
                'ACC','effN','SK','KT','M','SD');
    
else
    % load(fullfile(sim_dir,sprintf('PlotData_%s.mat',name)));
    load(fullfile(sim_dir,sprintf('CalculatedData_%s.mat',name)));
end



%% 
h = figure;
filename = sprintf('RT Distrib over a%i and b%i',numel(avec),numel(bvec));
subplot(2,3,1); hold on;
s=surf(Bmat,Amat,SK);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Skewness');
view(0,90);
c = colorbar;
% ylabel(c, 'Skewness');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,2); hold on;
s=surf(Bmat,Amat,KT);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Kurtosis');
view(0,90);
c = colorbar;
% ylabel(c, 'Kurtosis');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,4); hold on;
s=surf(Bmat,Amat,M);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Mean');
view(0,90);
c = colorbar;
% ylabel(c, 'Mean');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,5); hold on;
s=surf(Bmat,Amat,SD);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Deviance');
view(0,90);
c = colorbar;
% ylabel(c, 'Deviance');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,6); hold on;
s=surf(Bmat,Amat,ACC);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Accuracy');
view(0,90);
c = colorbar;
% ylabel(c, 'Accuracy');
xlabel('\beta');
ylabel('\alpha');
savefig(h,fullfile(plot_dir,[filename, '.fig']));
%% get an clear-cut boundary from analytical analysis
% mypool = parpool(7);
%% panel e, parameter space for choice/representation under equal inputs
w = 1;
v = 1;
scale = 2*1*eqlb.^2 + (1-Amat).*eqlb;
V1 = (1+cp)*scale;
V2 = (1-cp)*scale;
Nb = length(bvec);
filename = sprintf('Stability_Numeric_cp%1.3f_w%1.1f_v%1.1f_Amat%i_%1.2f_%1.2fBmat%i_%1.2f_%1.2f',...
    [cp, w, v, numel(avec),min(avec),max(avec),numel(bvec),min(bvec),max(bvec)]);
output = fullfile(sim_dir,[filename, '.mat']);
if ~exist(output, 'file')
    Npoints = NaN(size(Amat));
    Eigenvalue = -ones([size(Amat),4]);
    Stability = NaN(size(Amat));
    for ai = 1:length(avec)
        a = avec(ai);
        fprintf('%2.1f.',a);
        syms R1 R2;
        vars = [R1 R2];
        parfor bi = 1:Nb
            fprintf('.');
            b = bvec(bi);
            eqns = [(V1(bi,ai)/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
                (V2(bi,ai)/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
            [AnswR1,AnswR2] = solve(eqns, vars);
            AnswI1 = b*AnswR1;
            AnswI2 = b*AnswR2;
            AnswG1 = w*AnswR1 + v*AnswR2 - AnswI1;
            AnswG2 = v*AnswR1 + w*AnswR2 - AnswI2;
            PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
            Npoints(bi,ai) = sum(PositiveRealAnsw);
            NegEigen = -ones(1,4);
            for i = find(PositiveRealAnsw)'
                R1star = double(AnswR1(i));
                R2star = double(AnswR2(i));
                G1star = double(AnswG1(i));
                G2star = double(AnswG2(i));
                JMat = [-1 + a/(1+G1star), -(V1(bi,ai)+a*R1star)/(1+G1star)^2, 0, 0, 0, 0
                    w, -1, -1, v, 0, 0
                    b, 0, -1, 0, 0, 0
                    0, 0, 0, -1+a/(1+G2star),   -(V2(bi,ai)+a*R2star)/(1+G2star)^2, 0
                    v, 0, 0, w, -1, -1
                    0, 0, 0, b, 0, -1];
                A = eig(JMat);
                NegEigen(i) = all(real(A) < 0); % attractive 1, diversive 0
            end
            Stability(bi,ai)  = all(NegEigen>0);
            Eigenvalue(bi,ai,:) = NegEigen;
        end
        fprintf('\n');
    end
    visualize = sum(Eigenvalue,3);
    visualize(visualize == -2) = 2;
    save(output,'Amat','Bmat','avec','bvec','V1','V2','cp','w','v','Npoints','Eigenvalue','Stability','visualize');
else
    load(output);
end
% plot
dyadic = visualize == 2 | visualize == 0;
colorpalettergb = [239,71,111;255,209,102;6,214,160;17,138,178;7,59,76]/255;
h = figure; colormap(colorpalettergb([3,5],:));
% filename = 'RTACC_fullspace_Stability_Numeric';
s = surf(Bmat,Amat,dyadic+1,'EdgeColor','none');
% colorbar;
s.EdgeColor = 'none';
%set(gca,'YScale','log');
%ylim([1,length(avec)]);
%xlim([1,length(rvec)]);
%xticks(linspace(1,length(rvec),5));
%xticklabels(linspace(0,4,5));
%yticks(linspace(1,length(avec),5));
% yticklabels(linspace(0,60,5));
%yticklabels({'10^{-1}','10^0','10^1','10^2'});
xlabel('\beta');
ylabel('\alpha');
view(0,90);
savefigs(h, filename, plot_dir, 12, [2.8, 2.54]);

%% plot with mask
mask = NaN(size(dyadic));
mask(~dyadic) = 1;
h = figure;
filename = sprintf('RT Distrib over a%i and b%i [mask]',numel(avec),numel(bvec));
subplot(2,3,1); hold on;
s=surf(Bmat,Amat,SK.*mask);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Skewness');
view(0,90);
c = colorbar;
% ylabel(c, 'Skewness');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,2); hold on;
s=surf(Bmat,Amat,KT.*mask);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Kurtosis');
view(0,90);
c = colorbar;
% ylabel(c, 'Kurtosis');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,4); hold on;
s=surf(Bmat,Amat,M.*mask);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Mean');
view(0,90);
c = colorbar;
% ylabel(c, 'Mean');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,5); hold on;
s=surf(Bmat,Amat,SD.*mask);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Deviance');
view(0,90);
c = colorbar;
% ylabel(c, 'Deviance');
xlabel('\beta');
ylabel('\alpha');

subplot(2,3,6); hold on;
s=surf(Bmat,Amat,ACC.*mask);
s.EdgeColor = 'none';
%set(gca,'YScale','log');
xlim([0,1.8]);
title('Accuracy');
view(0,90);
c = colorbar;
% ylabel(c, 'Accuracy');
xlabel('\beta');
ylabel('\alpha');
savefig(h,fullfile(plot_dir,[filename, '.fig']));
%% plot RT distribution - fitted
rate = length(rtmat)/1024;
maxrt = max(max(rtmat));
minrt = min(min(rtmat));
%segrt = maxrt - minrt;
bank1 = [];
bank2 = [];
acc = [];
meanrtc = [];
meanrtw = [];
for ii = 1:6
    gap = (dataBhvr.rtrange(ii,2) - dataBhvr.rtrange(ii,1))/dataBhvr.numbins;
    %gap = 1.757/60;
    BinEdge = [minrt:gap:(maxrt+gap)];
    hg = histogram(rtmat(choicemat(:,ii)==1,ii),'BinEdges',BinEdge);
    bank1{ii} = hg.Values/rate;
    hg = histogram(rtmat(choicemat(:,ii)==2,ii),'BinEdges',BinEdge);
    bank2{ii}= hg.Values/rate;
    BinMiddle{ii} = hg.BinEdges(1:end-1) + hg.BinWidth/2;
    acc(ii) = sum(choicemat(:,ii)==1)/(sum(choicemat(:,ii)==1) + sum(choicemat(:,ii)==2));
    meanrtc(ii) = mean(rtmat(choicemat(:,ii)==1,ii));
    meanrtw(ii) = mean(rtmat(choicemat(:,ii)==2,ii));
end
% loading Roitman's data
addpath('../RoitmanDataCode');
ColumnNames608
load T1RT.mat;
x(:,R_RT) = x(:,R_RT)/1000;
cohlist = unique(x(:,R_COH));
maxrt = max(x(:,R_RT));
minrt = min(x(:,R_RT));
segrt = maxrt - minrt;
bins = 30;
BinEdge = [minrt:segrt/bins:maxrt];
bank1r = [];
bank2r = [];
accr = [];
meanrtcr = [];
meanrtwr = [];
for i = 1:length(cohlist)
    Lcoh = x(:,R_COH)==cohlist(i);
    if i == 1
        Dir1 = x(:,R_TRG) == 1;
        Dir2 = x(:,R_TRG) == 2;
        RT_corr = x(Lcoh & Dir1,R_RT);
        RT_wro = x(Lcoh & Dir2, R_RT);
    else
        Corr = x(:,R_DIR) == x(:,R_TRG);
        Wro = x(:,R_DIR) ~= x(:,R_TRG);
        RT_corr = x(Lcoh & Corr,R_RT);
        RT_wro = x(Lcoh & Wro, R_RT);
    end
    accr(i) = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    meanrtcr(i) = mean(RT_corr);
    meanrtwr(i) = mean(RT_wro);
    hg = histogram(RT_corr,'BinEdges',BinEdge);
    bank1r(:,i) = hg.Values;
    if ~isempty(RT_wro)
        hg = histogram(RT_wro,'BinEdges',BinEdge);
        bank2r(:,i) = hg.Values;
    else
        bank2r(:,i) = zeros(1,bins);
    end
end
BinMiddler = hg.BinEdges(1:end-1) + hg.BinWidth/2;
h = figure;
for ii = 1:6
    subplot(6,1,ii);
    % bar(BinMiddler,bank1r(:,ii),'FaceColor','#0072BD','EdgeAlpha',0);
    bar(dataBhvr.bincenter(ii,1:30),dataBhvr.histmat(ii,1:30)*1024,'FaceColor','#0072BD','EdgeAlpha',0);
    hold on;
    % bar(BinMiddler,-bank2r(:,ii),'FaceColor','#D95319','EdgeAlpha',0);
    %plot(BinMiddle{ii},bank1{ii},'c','LineWidth',1.5);
    %plot(BinMiddle{ii},-bank2{ii},'r','LineWidth',1.5);
    bar(dataBhvr.bincenter(ii,1:30),-dataBhvr.histmat(ii,31:60)*1024,'FaceColor','#D95319','EdgeAlpha',0,'EdgeColor','none');
    plot(BinMiddle{ii},bank1{ii},'c','LineWidth',2);
    plot(BinMiddle{ii},-bank2{ii},'m','LineWidth',2);
    if ii == 7
        legend({'','','Correct','Error'},'NumColumns',2,'Location','North');
        legend('boxoff');
    end
    ylim([-60,100]);
    yticks([-50:50:100]);
    yticklabels({'50','0','50','100'});
    xlim([100 1762]/1000);
    xticks([.5,1.0,1.5]);
    if ii == 6
        xticklabels({'.5','1.0','1.5'});
        xlabel('Reaction time (secs)');
    else
        xticklabels({});
    end
    if ii == 1
        ylabel('Frequency');
    end
    % title(sprintf('coherence %2.1f %%',cohlist(ii)*100));
    set(gca,'FontSize',16);
    set(gca,'TickDir','out');
    H = gca;
    H.LineWidth = 1;
    set(gca, 'box','off');
end
%set(gca,'FontSize',18);

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 3.0 10];
%saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.fig',name)),'fig');
saveas(h,fullfile(plot_dir,sprintf('RTDistrb_%s.eps',name)),'epsc2');
% Q-Q plot for reaction time and choice
lwd = 1.0;
mksz = 3;
fontsize = 11;
x = dataBhvr.proportionmat;
y = dataBhvr.q;
qntls = dataBhvr.qntls;
h = figure; hold on;
for vi = 1:length(x)
    xc = x(vi)*ones(size(y(:,1,vi)));
    xw = 1 - x(vi)*ones(size(y(:,2,vi)));
    lgc = plot(xc,y(:,1,vi),'gx','MarkerSize',mksz+1,'LineWidth',lwd);
    lge = plot(xw,y(:,2,vi),'rx','MarkerSize',mksz+1,'LineWidth',lwd);
    % fitted value
    En(vi) = numel(rtmat(:,vi));
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    xr = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
    q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
end
for qi = 1:size(q,1)
    xq = [flip(1-x), x]';
    plot(xq,[squeeze(flip(q(qi,2,:)));squeeze(q(qi,1,:))],'k-o','MarkerSize',mksz,'LineWidth',lwd/2);
end
legend([lge,lgc],{'error','correct'},"NumColumns",2,'Location','northeast','FontSize',fontsize-2);
legend('box','off');
xlim([-.05 1.05]);
ylim([.2, 1.4]);
yticks([.2:.4:1.4]);
xlabel('Proportion');
ylabel('RT (s)');
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 4 5];
% saveas(h,fullfile(plot_dir,sprintf('Q-QPlot_%s.fig',name)),'fig');
filename = sprintf('Q-QPlot_%s',name);
% saveas(h,fullfile(plot_dir,filename),'epsc2');
savefigs(h, filename, plot_dir, fontsize, [2.5 2.5]);
