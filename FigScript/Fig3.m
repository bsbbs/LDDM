% Fig 3. property #1, divisive normalization
%% Inputs
% value inputs in the ppaer Louie et al., 2011, Figure
% Vin = [65, 195, 260, 390
%     130, 130, 130, 130]';
% Vout = [260, 260, 260, 260, 260
%     130, 163, 195, 228, 260]'; as in paper Louie et al., 2011
Vin = [linspace(0,2,5)*scale0;
    ones(1,5)*scale0]';
Vout = [ones(1,5)*scale0;
    linspace(0,2,5)*scale0]';
h = figure; hold on;
filename = 'Fig3a_InputMatrix';
mygray5 = flip(gray(5 + 2));
for vi = 1:length(Vin)
    Vinput = Vin(vi,:);
    plot(Vinput(2),Vinput(1),'k.','MarkerSize',mksz/2,'Color',mygray5(vi+1,:));
end
mygray5 = flip(gray(5 + 2));
for vo = 1:length(Vout)
    Vinput = Vout(vo,:);
    plot(Vinput(2),Vinput(1),'k.','MarkerSize',mksz/2,'Color',mygray5(vo+1,:));
end
ylabel('V_1 (a.u.)');xlabel('V_2 (a.u.)');
xlim([0,500]);ylim([0,500]);
xticks([0,500]);yticks([0,500]);
set(gca,'FontSize',fontsize-7);
set(gca,'TickDir','out');
set(gca,'LineWidth',1);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 3.1 2.8]/2.4;
saveas(h,fullfile(plotdir,[filename,'.eps']),'epsc2');
%% panel a, dynamic of neural firing rates
a = 15*eye(2);
b = zeros(2);
w = ones(2);
initialvals = [2,2;4,4;0,0]*0;
predur = .1;
presentt = 0;
dur = 1;
stimdur = dur;
sgm = 0;
triggert = Inf;
thresh = Inf;
stoprule = 0;
h = figure;
filename = 'Fig3a_Dynamic';
mygray5 = flip(gray(5 + 2));
subplot(2,1,1); hold on;
R1vecVin = [];
for vi = 1:length(Vin)
    Vprior = [0, 0];
    Vinput = Vin(vi,:) + B0;
    [choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
        sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray5(vi+1,:), 'LineWidth',lwd);
    R1vecVin(vi) = mean(R((predur+presentt)/dt:((predur+presentt)/dt+1000),1));
end
ylim([-1,max(ylim)*1.3]);
yticks([0, 30, 60]);
xticks([(predur+presentt)/dt,predur/dt + 1000]);
xticklabels({});
savefigs(h, filename, plotdir, fontsize, [2, 4]);

subplot(2,1,2); hold on;
R1vecVout = [];
for vo = 1:length(Vout)
    Vprior = [0, 0];
    Vinput = Vout(vo,:) + B0;
    [choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
        sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray5(vo+1,:), 'LineWidth',lwd);
    R1vecVout(vo) = mean(R((predur+presentt)/dt:((predur+presentt)/dt+1000),1));
end
ylim([-1,max(ylim)*1.3]);
yticks([0, 30, 60]);
xticks([(predur+presentt)/dt,predur/dt + 1000]);
xticklabels({'0','1.0'});
ylabel('R_1 Activity (a.u.)');
xlabel('Time (a.u.)');
savefigs(h, filename, plotdir, fontsize, [2,4]);
%% panel b_left, nullclines for R1 and R2 under equal inputs
rng('default'); rng(5);
cp = 0;
V = [1+cp, 1-cp]*scale0 + B0;
a = a0;
b = 0;
w = 1;
v = 1;
% - time couse R1 & R2
sgm = 4; dur = 1;
presentt = dt; triggert = Inf; thresh = Inf;
initialvals = zeros(3,2); stimdur = dur; stoprule = 0;
% - Nullclines R1*-R2* space
h = figure; hold on;
filename = 'Fig3bL';
R2 = linspace(.1,40,200);
R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
plot([min(R2), V(2)/(1-a)],[V(1)/(1-a), V(1)/(1-a)],'--k','LineWidth',1);
Line1 = [R1' R2'];
R1 = linspace(.1,40,200);
R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
plot([V(2)/(1-a), V(2)/(1-a)],[min(R1), V(1)/(1-a)],'--k','LineWidth',1);
Line2 = [R1' R2'];
if 1
    plot([0.1,1000],[0.1,1000],'--k');
end
b = b*w;
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
        plot(R2star, R1star, '.', 'Color','#073b4c', 'MarkerSize',mksz);
        plot(R2star, R1star, 'o', 'Color','white', 'MarkerSize',mksz/3);
    else
        plot(R2star, R1star, '.', 'Color', '#06d6a0', 'MarkerSize',mksz);
    end
end
set(gca, 'XScale','log');
set(gca, 'YScale','log');
xlim([.1,5*10^3]);ylim([.1,5*10^3]);
xticks([1,10,100,1000]);yticks([1,10,100,1000]);
xlabel('R_2 activity (a.u.)');
ylabel('R_1 activity (a.u.)');
savefigs(h, filename, plotdir, fontsize, [2.8 2.54]);

%% panel b_middle, nullcines for R1 and R2 under moderately unequal inputs
rng('default'); rng(5);
cp = .512;
V = [1+cp 1-cp]*scale0 + B0; % 51.2 %
a = a0;
b = 0;
w = 1;
v = 1;
% - time couse R1 & R2
sgm = 4; dur = 1;
presentt = dt; triggert = Inf; thresh = Inf;
initialvals = zeros(3,2); stimdur = dur; stoprule = 0;
% - Nullclines R1*-R2* space
h = figure; hold on;
filename = 'Fig3bM';
R2 = linspace(.01,70,200);
R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
plot([min(R2), V(2)/(1-a)],[V(1)/(1-a), V(1)/(1-a)],'--k','LineWidth',1); %upbound
Line1 = [R1' R2'];
R1 = linspace(.01,60,200);
R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
plot([V(2)/(1-a), V(2)/(1-a)],[min(R1), V(1)/(1-a)],'--k','LineWidth',1); %upbound
Line2 = [R1' R2'];
if b == 0 || V(1) ~= V(2)
    plot([0.1,1000],[0.1,1000],'--k');
end
b = b*w;
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
        plot(R2star, R1star, '.', 'Color','#073b4c', 'MarkerSize',mksz);
        plot(R2star, R1star, 'o', 'Color','white', 'MarkerSize',mksz/3);
    else
        plot(R2star, R1star, '.', 'Color', '#06d6a0', 'MarkerSize',mksz);
    end
end
set(gca, 'XScale','log');
set(gca, 'YScale','log');
xlim([.1,5*10^3]);ylim([.1,5*10^3]);
xticks([1,10,100,1000]);yticks([1,10,100,1000]);
xlabel('R_2 activity (a.u.)');
ylabel('R_1 activity (a.u.)');
% legend([lgd1, lgd2],{'dR_1/dt = 0','dR_2/dt = 0'},...
%     'Location','NorthEast','FontSize',fontsize-5, 'FontName','Times New Roman', ...
%     'FontAngle','italic','NumColumns',1,'Box','off');
savefigs(h, filename, plotdir, fontsize, [2.8 2.54]);

%% panel b_right, nullcines for R1 and R2 under extremely unequal inputs
rng('default'); rng(5);
cp = 1;
V = [1+cp 1-cp]*scale0 + B0; % 100 %
a = a0;
b = 0;
w = 1;
v = 1;
% - time couse R1 & R2
sgm = 4; dur = 1;
presentt = dt; triggert = Inf; thresh = Inf;
initialvals = zeros(3,2); stimdur = dur; stoprule = 0;
% - Nullclines R1*-R2* space
h = figure; hold on;
filename = 'Fig3bR';
R2 = linspace(.01,32,200);
R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
plot([min(R2), V(2)/(1-a)],[V(1)/(1-a), V(1)/(1-a)],'--k','LineWidth',1); %upbound
Line1 = [R1' R2'];
R1 = linspace(.01,35,200);
R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
plot([V(2)/(1-a), V(2)/(1-a)],[min(R1), V(1)/(1-a)],'--k','LineWidth',1); %upbound
Line2 = [R1' R2'];
if b == 0 || V(1) ~= V(2)
    plot([0.1,30],[0.1,30],'--k');
end
b = b*w;
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
        plot(R2star, R1star, '.', 'Color','#073b4c', 'MarkerSize',mksz);
        plot(R2star, R1star, 'o', 'Color','white', 'MarkerSize',mksz/3);
    else
        plot(R2star, R1star, '.', 'Color', '#06d6a0', 'MarkerSize',mksz);
    end
end
set(gca, 'XScale','log');
set(gca, 'YScale','log');
xlim([.1,5*10^3]);ylim([.1,5*10^3]);
xticks([1,10,100,1000]);yticks([1,10,100,1000]);
xlabel('R_2 activity (a.u.)');
ylabel('R_1 activity (a.u.)');
legend([lgd1, lgd2],{'dR_1/dt = 0','dR_2/dt = 0'},...
    'Location','NorthEast','FontSize',fontsize-5, 'FontName','Times New Roman', ...
    'FontAngle','italic','NumColumns',1,'Box','off');
savefigs(h, filename, plotdir, fontsize, [2.8 2.54]);
