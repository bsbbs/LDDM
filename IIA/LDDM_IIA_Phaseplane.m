% Phase plane analysis for triplet system

%%
addpath(genpath('../CoreFunctions'));
addpath(genpath('../utils'));
outdir = './rslts/LDDM_IIA_Phaseplane';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
dt = .001;
Tau = ones(1,3)*.1;
lwd = 2.0;
mksz = 18;
%% panel c_middle, nullclines for R1 and R2 under unequal inputs
rng('default'); rng(5);
c = .256;
V = 256*[1+c, 1, 1.1];  
a = 10;
b = 1.8;
w = 1;
v = 1;
% - time couse R1, R2, and R3
sgm = 1; dur = 4;
presentt = dt; triggert = presentt; thresh = Inf;
initialvals = [4,4,4;12,12,12;0,0,0]/15; stimdur = dur; stoprule = 0;
% - Nullclines R1*-R2*-R3* space
h = figure; hold on;
filename = sprintf('LDDMTriplet_a%1.1f_b%1.1f_w%1.1f_v%1.1f_V1%3.0f_V2%3.0f_V3%3.0f',a,b,w,v,V);
[R3, R2] = meshgrid(linspace(.1,155,200),linspace(.1,155,200));
R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v; % dR3/dt = 0
lgd3 = surf(R1,R2,R3,'FaceAlpha',.3); % dR3/dt = 0
lgd3.EdgeColor = 'none';
lgd3.FaceColor = '#ef476f';

[R2, R1] = meshgrid(linspace(.1,155,200),linspace(.1,155,200));
R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v; % dR2/dt = 0
lgd2 = surf(R1,R2,R3,'FaceAlpha',.3); % dR2/dt = 0
lgd2.EdgeColor = 'none';
lgd2.FaceColor = '#118ab2';


[R1, R3] = meshgrid(linspace(.1,155,200),linspace(.1,155,200));
R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v; % dR1/dt = 0
lgd1 = surf(R1,R2,R3,'FaceAlpha',.3); % dR1/dt = 0
lgd1.EdgeColor = 'none';
lgd1.FaceColor = '#ffd166';

if 1
    [R, G, I, ~, ~] = LcDsInhbt(V, [w,v,v;v,w,v;v,v,w], a*eye(3), b*eye(3),...
        sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    ldgtrc(1) = plot3(R(round(presentt/dt):end,1), R(round(presentt/dt):end,2),R(round(presentt/dt):end,3),'-','Color','#ef476f','LineWidth',lwd/2);
    [R, G, I, ~, ~] = LcDsInhbt(V, [w,v,v;v,w,v;v,v,w], a*eye(3), b*eye(3),...
        sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    ldgtrc(2) = plot3(R(round(presentt/dt):end,1), R(round(presentt/dt):end,2),R(round(presentt/dt):end,3),'-','Color','#118ab2','LineWidth',lwd/2);
    [R, G, I, ~, ~] = LcDsInhbt(V, [w,v,v;v,w,v;v,v,w], a*eye(3), b*eye(3),...
        sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    ldgtrc(3) = plot3(R(round(presentt/dt):end,1), R(round(presentt/dt):end,2),R(round(presentt/dt):end,3),'-','Color','#073b4c','LineWidth',lwd/2);
end
syms R1 R2 R3
eqns = [(V(1)/R1 - (w - b)*R1 - (1-a) - v*R3)/v == R2, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
    (V(2)/R2 - (w - b)*R2 - (1-a) - v*R1)/v == R3, ... % dR2/dt = 0: R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v
    (V(3)/R3 - (w - b)*R3 - (1-a) - v*R2)/v == R1]; % dR3/dt = 0: R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v;
vars = [R1 R2 R3];
[AnswR1,AnswR2,AnswR3] = solve(eqns, vars);
AnswI1 = b*AnswR1;
AnswI2 = b*AnswR2;
AnswI3 = b*AnswR3;
AnswG1 = w*AnswR1 + v*AnswR2 + v*AnswR3  - AnswI1;
AnswG2 = v*AnswR1 + w*AnswR2 + v*AnswR3  - AnswI2;
AnswG3 = v*AnswR1 + v*AnswR2 + w*AnswR3 - AnswI3;
PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0 & double(AnswR2) >0 & double(imag(AnswR2)) == 0 & double(AnswR3) >0 & double(imag(AnswR3)) == 0;
Npoints = sum(PositiveRealAnsw);
lgdvec = [];
for i = [find(PositiveRealAnsw)]'
    R1star = double(AnswR1(i));
    R2star = double(AnswR2(i));
    R3star = double(AnswR3(i));
    G1star = double(AnswG1(i));
    G2star = double(AnswG2(i));
    G3star = double(AnswG3(i));
    JMat = [-1 + a/(1+G1star), -(V(1)+a*R1star)/(1+G1star)^2, 0,    0, 0, 0,       0, 0, 0
            w, -1, -1,                                              v, 0, 0,       v, 0, 0
            b, 0, -1,                                               0, 0, 0,       0, 0, 0
            0, 0, 0,    -1+a/(1+G2star),   -(V(2)+a*R2star)/(1+G2star)^2, 0,       0, 0, 0
            v, 0, 0,    w, -1, -1,                                                 v, 0, 0
            0, 0, 0,    b, 0, -1,                                                  0, 0, 0
            0, 0, 0,    0, 0, 0,       -1 + a/(1+G3star), -(V(3)+a*R3star)/(1+G3star)^2, 0
            v, 0, 0,    v, 0, 0,                                                  w, -1, -1
            0, 0, 0,    0, 0, 0,                                                  b, 0, -1];
    A = eig(JMat);
    Stability(i) = all(real(A) < 0); % attractive 1, diversive 0
    if Stability(i)
        lgdvec(end +1) = plot3(R1star, R2star, R3star, '.', 'Color','#06d6a0', 'MarkerSize',mksz);
        plot3(R1star, R2star, R3star, 'o', 'Color','k', 'MarkerSize',mksz/3);
    else
        lgdvec(end +1) = plot3(R1star, R2star, R3star, '.', 'Color', '#ffd166', 'MarkerSize',mksz);
        plot3(R1star, R2star, R3star, 'o', 'Color','k', 'MarkerSize',mksz/3);
    end
end
set(gca, 'XScale','log');
set(gca, 'YScale','log');
set(gca, 'ZScale','log');
xlim([.2,10^2]);ylim([.2,10^2]);zlim([.2,10^2]);
xticks([1,10,100]);yticks([1,10,100]);
xlabel('R_1 activity (a.u.)');
ylabel('R_2 activity (a.u.)');
zlabel('R_3 activity (a.u.)');
view([0,90]);
legend([lgd1, lgd2, lgd3, lgdvec(1)],{'dR_1/dt = 0','dR_2/dt = 0','dR_2/dt = 0','unstable points'}, 'FontName','Times New Roman', ...
    'FontAngle','Italic','FontSize',fontsize-6, 'Location','northeastoutside','Box','off');
aspect3 = [2.8 2.54]; % 1:1 for phase plane
aspect5 = [4.3 2.54]; % 1:1 for phase plane with legend outside
savefigs(h, filename, outdir, fontsize, aspect5);
grid on;
view([-52,8]);
savefigs(h, [filename '3D'], outdir, fontsize, aspect5*2);
savefig(h,fullfile(outdir,[filename '3D.fig']));
