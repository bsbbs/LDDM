% Fig 8. Persistent activity
%% Inputs
Vin = [linspace(0,2,5)*scale0;
    ones(1,5)*scale0]';
Vout = [ones(1,5)*scale0;
    linspace(0,2,5)*scale0]';
%% panel a, dynamic of neural firing rates
a = a0*eye(2);
b = zeros(2);
w = ones(2);
initialvals = [2,2;4,4;0,0]*0;
predur = .1;
presentt = 0;
dur = 1.7;
stimdur = .8;
sgm = 0;
triggert = Inf;
thresh = Inf;
stoprule = 0;
h = figure;
filename = 'Fig7a';
mygray5 = flip(gray(5 + 2));
subplot(2,1,1); hold on;
R1vecVin = [];
for vi = 1:length(Vin)
    Vprior = [0, 0];
    Vinput = Vin(vi,:) + B0;
    [choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
        sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray5(vi+1,:), 'LineWidth',1.5);
    R1vecVin(vi) = mean(R((predur+presentt)/dt:((predur+presentt)/dt+1000),1));
end
ylim([-1,70]);
yticks([0, 30, 60]);
xticks([(predur+presentt)/dt, (predur + presentt + stimdur)/dt]);
xlim([1, (predur + dur)/dt]);
xticklabels({});
savefigs(h, filename, plotdir, fontsize, [3, 4]);

subplot(2,1,2); hold on;
R1vecVout = [];
for vo = 1:length(Vout)
    Vprior = [0, 0];
    Vinput = Vout(vo,:) + B0;
    [choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
        sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    lgd1 = plot(R(:,1), 'k-', 'Color', mygray5(vo+1,:), 'LineWidth',1.5);
    R1vecVout(vo) = mean(R((predur+presentt)/dt:((predur+presentt)/dt+1000),1));
end
ylim([-1,70]);
yticks([0, 30, 60]);
xticks([(predur+presentt)/dt, (predur + presentt + stimdur)/dt]);
xlim([1, (predur + dur)/dt]);
xticklabels({'on','off'});
ylabel('R_1 Activity (a.u.)');
xlabel('Time (a.u.)');
savefigs(h, filename, plotdir, fontsize - 2, [3,4]);
%% panel b, nullclines for R1 and R2 under any inputs
i = 2;
Taup = [1,1,1];
alpha = eye(2)*a0;
if i == 1
    w = [1,.5;.5,1];
elseif i == 2
    w = [1,1;1,1];
elseif i == 3
    w = [1,2;2,1];
end
beta = 0; % 0, .4, 1.2
Dt = 1;
if w(1,1) == w(1,2)
    if beta > 0 && beta < 1
        axismax = 25;
    elseif beta == 0
        axismax = 18;
    elseif beta >= 1
        axismax = 50;
    end
else
    if w(1,1) > w(1,2)
        axismax = (alpha(1,1) - 1)/w(1,2)*1.1;
    elseif w(1,1) < w(1,2)
        axismax = (alpha(1,1) - 1)/w(1,1)*1.1;
    end
end
[R1,R2] = meshgrid(linspace(0,axismax,15),linspace(0,axismax,15));
R1(1,1) = R1(1,2)/3;
R2(1,1) = R2(2,1)/3;
G1 = R1*w(1,1) + R2*w(1,2) - beta*R1;
G2 = R1*w(2,1) + R2*w(2,2) - beta*R2;
G1(G1<0) = 0;
G2(G2<0) = 0;
dR1 = (-R1 + (R1.*alpha(1,1))./(1+G1)) * Dt/Tau(1);
dR2 = (-R2 + (R2.*alpha(2,2))./(1+G2)) * Dt/Tau(2);
norm = sqrt(dR1.^2 + dR2.^2);
rate = (norm.^.6)./(norm);
dR1 = dR1.*rate;
dR2 = dR2.*rate;
h = figure;
filename = 'Fig7b';
hold on;
if w(2,2) - beta > 0 && w(1,2) == w(1,1)
    R2Line = [0,(alpha(2,2) - 1)/(w(2,2) - beta)];
elseif w(2,2) - beta <= 0 && w(1,2) == w(1,1)
    R2Line = [0,axismax];
elseif beta == 0 && w(1,2) ~= w(1,1)
    R2Line = [0,(alpha(1,1)-1)/w(1,2)];
end
R1Line = -w(1,2)/(w(1,1)-beta).*R2Line+(alpha(1,1) - 1)/(w(1,1)-beta);
lg1 = plot(R2Line,R1Line,'-','Color',colorpalette{4},'LineWidth',lwd);% dR1/dt = 0
if w(2,2) - beta > 0 && w(1,2) == w(1,1)
    R1Line = [0,(alpha(1,1) - 1)/(w(1,1) - beta)];
elseif w(2,2) - beta <= 0 && w(1,2) == w(1,1)
    R1Line = [0,axismax];
elseif beta == 0 && w(1,2) ~= w(1,1)
    R1Line = [0,(alpha(2,2)-1)/w(2,1)];
end
R2Line = -w(2,1)/(w(2,2)-beta).*R2Line+(alpha(2,2) - 1)/(w(2,2)-beta);
lg2 = plot(R2Line,R1Line,'--','Color',colorpalette{1},'LineWidth',lwd);% dR2/dt = 0
if w(1,1) == w(1,2)
    if beta > 0 && beta < 1
        plot([(alpha(1,1)-1)/(w(1,1) - beta),(alpha(1,1)-1)/(w(1,1) - beta)],[0, (alpha(2,2)-1)/(w(2,2) - beta)],'k--','LineWidth',1); % constraints
        plot([0, (alpha(1,1)-1)/(w(1,1) - beta)],[(alpha(2,2)-1)/(w(2,2) - beta),(alpha(2,2)-1)/(w(2,2) - beta)],'k--','LineWidth',1); % constraints
        xlim([-1,axismax]);ylim([-1,axismax]);
        xticks([0, (alpha(1,1)-1)/(w(1,1)+w(1,2)-beta), (alpha(1,1)-1)/(w(1,1) - beta)]);
        yticks([0, (alpha(2,2)-1)/(w(2,2)+w(2,1)-beta), (alpha(2,2)-1)/(w(2,2) - beta)]);
        xticklabels({'0','$\frac{\alpha-1}{2\omega-\beta}$','$\frac{\alpha-1}{\omega-\beta}$'});
        yticklabels({'0','$\frac{\alpha-1}{2\omega-\beta}$','$\frac{\alpha-1}{\omega-\beta}$'});
    elseif beta == 0
        xlim([-1,axismax]);ylim([-1,axismax]);
        xticks([0, (alpha(1,1)-1)/(w(1,1) - beta)]);
        yticks([0, (alpha(2,2)-1)/(w(2,2) - beta)]);
        xticklabels({'0','$\frac{\alpha-1}{\omega}$'});
        yticklabels({'0','$\frac{\alpha-1}{\omega}$'});
    elseif beta >= 1
        xlim([-1,axismax]);ylim([-1,axismax]);
        xticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2) - beta)]);
        yticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1) - beta)]);
        xticklabels({'0','$\frac{\alpha-1}{2\omega-\beta}$'});
        yticklabels({'0','$\frac{\alpha-1}{2\omega-\beta}$'});
    end
elseif w(1,1) ~= w(1,2) && beta == 0
    if w(1,1) > w(1,2) % convergent
        xlim([-1,axismax]);ylim([-1,axismax]);
        xticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2)),(alpha(1,1)-1)/w(1,2)]);
        yticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1)),(alpha(2,2)-1)/w(2,1)]);
        xticklabels({'0','$\frac{\alpha-1}{w+v}$','$\frac{\alpha-1}{v}$'});
        yticklabels({'0','$\frac{\alpha-1}{w+v}$','$\frac{\alpha-1}{v}$'});
    elseif w(1,1) < w(1,2) % divergent
        xlim([-1,axismax]);ylim([-1,axismax]);
        xticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2)),(alpha(1,1)-1)/w(1,1)]);
        yticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1)),(alpha(2,2)-1)/w(2,2)]);
        xticklabels({'0','$\frac{\alpha-1}{w+v}$','$\frac{\alpha-1}{w}$'});
        yticklabels({'0','$\frac{\alpha-1}{w+v}$','$\frac{\alpha-1}{w}$'});
    end
end
set(gca,'TickLabelInterpreter', 'latex');
xlabel('R_2 activity (a.u.)');
ylabel('R_1 activity (a.u.)');
quiver(R2,R1,dR2,dR1,1.3,'r','LineWidth',1);
legend({'\color[rgb]{.0667,.5412,.6980}dR_1/dt = 0','\color[rgb]{.9373,.2784,.4353}dR_2/dt = 0','\color{red}Change rate'}, 'FontName','Times New Roman', ...
    'FontAngle','Italic','FontSize',fontsize-5, 'Location','northeastoutside','Box','off');
savefigs(h, filename, plotdir, fontsize, [4.1, 2.54]);
%% panel c, coded ratio, comparing with the Wang type model
h = figure; hold on;
filename = 'Fig7c';
% Wang type model
paramspecify_WongWang;
dt = .001; sgm = 0; dur = 4; presentt = dt;stimdur = 2;stoprule = 0;
name = sprintf('CodedRatioWong06_Sim%i_dur%1.1f_sgm%1.2f',length(V1Iter), dur, sigma);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iter)
        fprintf('V1 %3.1f',V1Iter(ii));
        fprintf('.');
        Vinput = [V1Iter(ii), V2Iter(ii)];
        rescale = 256/mean(Vinput);
        [nu_wind, s_wind, rt, choice, H, S] = wong06(Vinput*rescale,miu0,sgm,I0,JN,...
            gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
        R1rp(ii) = nu_wind(round((presentt+stimdur)/dt) - 1,1);
        R2rp(ii) = nu_wind(round((presentt+stimdur)/dt) - 1,2);
        R1wm(ii) = nu_wind(round((dur)/dt) - 1,1);
        R2wm(ii) = nu_wind(round((dur)/dt) - 1,2);
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm','V1Iter','V2Iter');
else
    load(output);
end
lgd1 = plot(V1Iter(:)./(V1Iter(:)+V2Iter(:)), R1wm(:)./(R1wm(:)+R2wm(:)),'.','Color',colorpalette{1},'MarkerSize',mksz/2);
plot(V1Iter(:)./(V1Iter(:)+V2Iter(:)), R1wm(:)./(R1wm(:)+R2wm(:)),'-','Color',colorpalette{1});
% Disinhibitory model
sigma = 0;
alpha = eye(2)*15;
w = [1,1; 1,1];
dur = 4;
stimdur = 2;
name = sprintf('CodedRatioSim%i_dur%1.1f_W%1.2f%1.2f_alpha%1.1f_sgm%1.2f',length(V1Iter), dur, w(1,1), w(1,2), alpha(1), sigma);
output = fullfile(datadir,[name '.mat']);
if exist(output,'file')
    load(output);
else
    error('Please check simulation in Fig3 - panel d');
end
lgd2 = plot(V1Iter(:)./(V1Iter(:)+V2Iter(:)), R1wm(:)./(R1wm(:)+R2wm(:)),'.','Color',colorpalette{5},'MarkerSize',mksz/2);
plot(V1Iter(:)./(V1Iter(:)+V2Iter(:)), R1wm(:)./(R1wm(:)+R2wm(:)),'-','Color',colorpalette{5});
xlabel('Input ratio (V_1 vs. V_2)');ylabel('Coded ratio (R_1 vs. R_2)');
xticks([0,.25,.5,.75,1]); yticks([0,.25,.5,.75,1]);
legend([lgd2, lgd1],{'\color[rgb]{.0275,.2314,.298}The hybird model',...
    '\color[rgb]{.9373,.2784,.4353}Wong&Wang, 2006'},'Location','SouthEast','FontSize',fontsize-5,'Box','off');
savefigs(h, filename, outdir, fontsize, aspect3);