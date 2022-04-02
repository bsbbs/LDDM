%% Phase-plane analysis
%% define paths
Homedir = '~';
addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'Documents','LDDM','utils'));
% cd('G:\My Drive\LDDM\Froemke');
cd('/Volumes/GoogleDrive/My Drive/LDDM/Froemke/SSTVIP_NMDA/');
plotdir = fullfile('./Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = './SimRslts';
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end

% parameters
fontsize = 14;
mksz = 18;
lwd = 1.5;
aspect1 = [3.9,2.2]; % 16:9, for wide temporal dynamic
aspect2 = [3 3]; % for temporal dynamic
aspect3 = [2.8 2.54];

%% nullclines for R1 and R2 under unequal inputs
rng('default'); rng(2);
c = .8; % .872;
a0 = 10;
b0 = 0;%.7;
w0 = 1;
v0 = 1;
% - time couse R1 & R2
dt = .001;Tau = [.1,.1,.1];
sgm = .1; dur = 2;
presentt = dt; triggert = presentt; thresh = Inf;
initialvals = [4,4;8,8;0,0]/15; stimdur = dur; stoprule = 0;
% - Nullclines R1*-R2* space
h = figure; hold on;
filename = 'PhasePlane_Represent';
boost = [1,2];
mycl = jet(length(boost)*10);
for level = 1:2
    V = 30*[1+c, 1-c]; %*boost(level);
    a = a0; %*boost(level);
    b = b0*boost(level);
    w = w0*boost(level);
    v = v0*boost(level);
    R2 = linspace(.1,35,200);
    R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
    lgd2(level) = plot(R1,R2,'--', 'Color',mycl(5+(level-1)*10,:),'LineWidth',lwd); % dR2/dt = 0
    plot([V(1)/(1-a), V(1)/(1-a)],[min(R2), V(2)/(1-a)],'--k','LineWidth',1);
    Line1 = [R1' R2'];
    R1 = linspace(.1,35,200);
    R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
    lgd1(level) = plot(R1,R2,'-', 'Color',mycl(5+(level-1)*10,:),'LineWidth',lwd); % dR1/dt = 0
    plot([min(R1), V(1)/(1-a)],[V(2)/(1-a), V(2)/(1-a)],'--k','LineWidth',1);
    Line2 = [R1' R2'];
    if 1
        [R, G, I, ~, ~] = LcDsInhbt(V, [w,v;v,w], a*eye(2), b*eye(2),...
            sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        ldgtrc(level) = plot(R(round(presentt/dt):end,1), R(round(presentt/dt):end,2),'-','Color',mycl(5+(level-1)*10,:),'LineWidth',lwd/2); % '#ef476f'
%         [R, G, I, ~, ~] = LcDsInhbt(V, [w,v;v,w], a*eye(2), b*eye(2),...
%             sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%         ldgtrc(2) = plot(R(round(presentt/dt):end,1), R(round(presentt/dt):end,2),'-','Color','#118ab2','LineWidth',lwd/2);
    end
    if 1 %b == 0 || V(1) ~= V(2)
        plot([0.1,1000],[0.1,1000],'--k');
    end
    syms R1 R2
    eqns = [(V(1)/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
        (V(2)/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
    vars = [R1 R2];
    [AnswR1,AnswR2] = solve(eqns, vars);
    AnswI1 = b*AnswR1;
    AnswI2 = b*AnswR2;
    AnswG1 = w*AnswR1 + v*AnswR2 - AnswI1;
    AnswG2 = v*AnswR1 + w*AnswR2 - AnswI2;
    PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
    Npoints = sum(PositiveRealAnsw);
    lgdvec = [];
    for i = [find(PositiveRealAnsw)]'
        R1star = double(AnswR1(i));
        R2star = double(AnswR2(i));
        G1star = double(AnswG1(i));
        G2star = double(AnswG2(i));
        JMat = [-1 + a/(1+G1star), -(V(1)+a*R1star)/(1+G1star)^2, 0, 0, 0, 0
            w, -1, -1, v, 0, 0
            b, 0, -1, 0, 0, 0
            0, 0, 0, -1+a/(1+G2star),   -(V(2)+a*R2star)/(1+G2star)^2, 0
            v, 0, 0, w, -1, -1
            0, 0, 0, b, 0, -1];
        A = eig(JMat);
        Stability(i) = all(real(A) < 0); % attractive 1, diversive 0
        if Stability(i)
            lgdvec(end+1)=plot(R1star, R2star, '.', 'Color','#06d6a0', 'MarkerSize',mksz);
            plot(R1star, R2star, 'o', 'Color','k', 'MarkerSize',mksz/3);
        else
            lgdvec(end+1)=plot(R1star, R2star, '.', 'Color', '#ffd166', 'MarkerSize',mksz);
            plot(R1star, R2star, 'o', 'Color','k', 'MarkerSize',mksz/3);
        end
    end
end
set(gca, 'XScale','log');
set(gca, 'YScale','log');
xlim([.1,10^2]);ylim([.1,10^2]);
xticks([1,10,100]);yticks([1,10,100]);
xlabel('R_1 activity (a.u.)');
ylabel('R_2 activity (a.u.)');
lgd = legend([lgd1(1),lgd2(1),ldgtrc(1), lgd1(2),lgd2(2),ldgtrc(2)],{'','','','dR_1/dt = 0','dR_2/dt = 0','trace'},...
    'NumColumns',2,'FontName','Times New Roman', ...
    'FontAngle','Italic','FontSize',fontsize-8, 'Location','northeastoutside','Box','off');
title(lgd, "Baseline                      iSTDP                 .");
savefigs(h, filename, plotdir, fontsize, aspect3 + [1.9,0]);
