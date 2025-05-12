%% Singleton - representation
%% Common parameters
N = 1;
w = w0*ones(N);
a = a0*eye(N);
b = b0*eye(N);
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
mygray = [.62, .62, .62; 0, 0, 0];
mygreen = [.62, 1, .62; 0, 1, 0];
myred = [1, .62, .62; 1, 0, 0];
myalpha = [0.7, 0];
cp = 1;
stoprule = 1;


%% Rate-based Plasticity 
Vprior = 0;
Vinput = 100;
BR = 30;
BG = 30;
dur = 1500; %s
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
R0 = BR*ones(N,1);
D0 = 0*ones(N,1);
G0 = BG*ones(N,1);
initialvals = [R0'; G0'; D0'];
winput0 = rand(1,N)*0;
wrg0 = rand(N,N)*0;
wgr0 = rand(1,N)*0;
sgmR = 0; %12;
sgmG = 0; %12;
sgmInput = 0; %66;
[R, G, D, Vcourse, winputc, wrgc, wgrc] = LDDM_RndInputPlastic(Vprior, Vinput, BR, BG, winput0, wrg0, wgr0, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
winputc(end) % 3.1250
wrgc(end) % .7352; .6357
wgrc(end) % .2229; .1927
%%
h = figure; 
filename = sprintf('MeanFR_N%i', N);
rate = 1000;
x1 = 1:rate:(6*60/dt);
x2 = (24*60/dt):rate:(25*60/dt);
x = [x1, x2];
hold on;
ldg = [];
ldg(1) = plot(x1*dt,R(x1), 'g-', 'LineWidth',lwd);
ldg(2) = plot(x1*dt,G(x1), 'r-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,R(x2), 'g-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,G(x2), 'r-', 'LineWidth',lwd);
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ylim([25,55]);
legend(ldg, {'R','G'}, "Location","east");
xlabel('Time (mins)');
ylabel('Firing rates (Hz)');
mysavefig(h, filename, plotdir, fontsize, [2.4,1.86], .75);


h = figure; 
filename = sprintf('WDynamic_N%i', N);
hold on;
ldg = [];
yyaxis left;
ldg(1) = plot(x1*dt, winputc(x1), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
ldg(2) = plot(x1*dt, wrgc(x1), 'g-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,winputc(x2), '-', 'Color', [.38,.38,.38], 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wrgc(x2), 'g-', 'LineWidth',lwd);
ylabel('Excitatory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ax = gca;
ax.YColor = [0,0,0];
yyaxis right;
ldg(3) = plot(x1*dt, wgrc(x1), 'r-', 'LineWidth',lwd);
plot(50+(x2-min(x2)+max(x1))*dt,wgrc(x2), 'r-', 'LineWidth',lwd);
ylabel('Inhibitory weight');
xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
xticklabels({'0','2','4','6','24','25'});
xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
ax.YColor = [0,0,0];
legend(ldg, {'w_{Input}','w_{R to G}', 'w_{G to R}'}, "Location","east");
xlabel('Time (mins)');
mysavefig(h, filename, plotdir, fontsize, [2.59,1.86], .75);
%% Example trace
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
task = sprintf('Trace_N%i',N);
h = figure;
hold on;
filename = sprintf('%s_BeforetoAfter',task);
sgmR = 0;
sgmG = 0;
sgmInput = 0;
winputs = [0, winputc(end)]; %3.1250;
wrgs = [0, wrgc(end)]; % .6357
wgrs = [0, wgrc(end)]; % .1927
for i = [1,2]
    winput = winputs(i);
    wrg = wrgs(i);
    wgr = wgrs(i);
    R0 = 0*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = 0*ones(N,1);
    initialvals = [R0'; G0'; D0']; %zeros(3,N);
    Vprior = zeros(size(cp));
    Vinput = 100*ones(1,N).*cp;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    subplot(2,2,i);hold on;
    ldg = [];
    x = (1:numel(R(:,1)))*dt;
    ldg(1) = scatter(x(1:5:end), R(1:5:end,1), 10, 'MarkerFaceColor', 'none', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor','green');
    ldg(2) = scatter(x(1:5:end), G(1:5:end,1), 4, 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor','none');
    xlim([-50*dt, dur]);
    ylim([-5, 65]);
    ylabel('Firing rates (Hz)');
    xlabel('Time (s)');
    if i == 2
        legend(ldg, {'R','G'}, 'Location', 'southeast');
    end
    mysavefig(h, filename, plotdir, fontsize, [4.8,4.4]);
    
    subplot(2,2,2+i);hold on;
    scatter(G(1:5:end,1), R(1:5:end,1), 8, 'MarkerFaceColor', 'black', 'MarkerFaceAlpha', .3, 'MarkerEdgeColor','none');
    % scatter(G(:,1), R(:,1), '-', 'Color', mygray(i,:), 'LineWidth', lwd);
    %xlim([-max([R; G])*.08, max(G)*1.4]);
    %ylim([-max([R; G])*.08, max(R)*1.4]);
    xlim([-5, 65]);
    ylim([-5, 65]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    mysavefig(h, filename, plotdir, fontsize, [4.8,4.4]);
end
%% Vector field
task = sprintf('Vctrfld_N%i',N);
filename = sprintf('%s_BeforetoAfter',task);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
sgmR = 0;
sgmG = 0;
sgmInput = 0;
BR = 30;
BG = 30;
R0 = 3:5:70;
G0 = 3:5:100;
[Rm, Gm] = meshgrid(R0, G0);
h = figure;
for i = [1,2]
    subplot(1,2,i);hold on;
    winput = winputs(i);
    wrg = wrgs(i);
    wgr = wgrs(i);
    R0 = 0*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = 0*ones(N,1);
    initialvals = [R0'; G0'; D0']; %zeros(3,N);
    Vprior = zeros(size(cp));
    Vinput = 100*ones(1,N).*cp;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    dRm = -Rm + (winput.*Vinput + a.*Rm + BR)./(1 + wgr*Gm);
    dGm = -Gm + wrg.*Rm + BG;
    rate = 0.12;
    dRm = dRm.*rate;
    dGm = dGm.*rate;
    quiver(Gm,Rm,dGm,dRm,'off','Color', mygray(1,:),'LineWidth',1);
    scatter(G(1:5:end,1), R(1:5:end,1), 8, 'MarkerFaceColor', colorpalette{4}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor','none');
    xlim([0, 65]);
    ylim([0, 65]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    mysavefig(h, filename, plotdir, fontsize, [4.8,2.1], 0);
end

%% Vector field Surface
task = sprintf('VctrfldSurf_N%i',N);
filename = sprintf('%s_BeforetoAfter',task);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
sgmR = 0;
sgmG = 0;
sgmInput = 0;
BR = 30;
BG = 30;

h = figure;
for i = [1,2]
    subplot(1,2,i);hold on;
    winput = winputs(i);
    wrg = wrgs(i);
    wgr = wgrs(i);
    if i == 1
        R0 = 20:.5:44;
        G0 = 20:.5:40;
    elseif i == 2
        R0 = 20:.5:44;
        G0 = 40:.5:60;
    end
    [Rm, Gm] = meshgrid(R0, G0);
    dRm = -Rm + (winput.*Vinput + a.*Rm + BR)./(1 + wgr*Gm);
    dGm = -Gm + wrg.*Rm + BG;
    rate = 0.12;
    dRm = dRm.*rate;
    dGm = dGm.*rate;
    F = sqrt(dRm.^2 + dGm.^2);
    surf(Gm,Rm,F,'EdgeColor','none');
    grid on;
    xlim([0, 65]);
    ylim([0, 65]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    % view([30,90]);
    mysavefig(h, filename, plotdir, fontsize, [4.8,2.1], 0);
end
%% Distribution of representation
task = sprintf('Rprsnt_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 60; % second
h = figure;
filename = sprintf('%s_BeforetoAfter',task);
names = ["Early noise", "Late noise", "Interneuronal noise"];
for testi = 1:3
    switch testi
        case 1
            sgmInput = 66;
            sgmR = 0;
            sgmG = 0;
        case 2
            sgmInput = 0;
            sgmR = 18;
            sgmG = 0;
        case 3
            sgmInput = 0;
            sgmR = 0;
            sgmG = 18;
    end
    subplot(1,3,testi);
    hold on;
    % early noise
    rng(2025);
    mylgd = [];
    for i = [1,2]
        winput = winputs(i);
        wrg = wrgs(i);
        wgr = wgrs(i);
        R0 = eqlb*ones(N,1);
        D0 = 0*ones(N,1);
        G0 = wrg*R0+BG;
        initialvals = [R0'; G0'; D0'];
        Vprior = 100*ones(size(cp));
        Vinput = 100*ones(1,N).*cp;
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        hst1 = histogram(R(5000:end,1),...
            'EdgeColor', 'none', 'FaceColor', mygray(i,:), 'FaceAlpha', .3, 'Normalization', 'pdf');
        pd1 = fitdist(R(:,1),'kernel','Kernel','normal');
        x = hst1.BinEdges;
        y1 = pdf(pd1,x);
        mylgd(i) = plot(x,y1, '-', 'Color', colorpalette{i}, 'LineWidth', lwd);
    end
    ylabel('Density');
    xlabel('Activity (Hz)');
    lgd = legend(mylgd, {'Before','After'},...
        'Box','off','Location','east',...
        'FontName','Arial','FontSize',fontsize-4);
    title(lgd, 'LTP');
    title(names(testi));
    mysavefig(h, filename, plotdir, fontsize, [9,2.2]);
end

%%
% Late noise
sgmR = 18;
sgmG = 0;
sgmInput = 0;
subplot(1,3,2);
hold on;
rng(2025);
mylgd = [];
for i = 1:length(potentiation)
    eltp = potentiation(i);
    iltp = potentiation(i);
    balance_scale = iltp*eltp*mean(w*ones(N,1))*eqlb^2 + (iltp*BG+1-mean(a*ones(N,1)))*eqlb - BR;
    R0 = eqlb*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = eltp*w*R0+BG;
    initialvals = [R0'; G0'; D0'];
    Vprior = ones(size(cp));
    Vinput = ones(1,N).*cp;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, balance_scale, w, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);

    hst1 = histogram(R(5000:end,1),...
        'EdgeColor', 'none', 'FaceColor', colorpalette{i}, 'FaceAlpha', .3, 'Normalization', 'pdf');
    pd1 = fitdist(R(:,1),'kernel','Kernel','normal');
    x = hst1.BinEdges;
    y1 = pdf(pd1,x);
    plot(x, y1, '-', 'Color', colorpalette{i}, 'LineWidth', lwd);
end
ylabel('Density');
xlabel('Activity (Hz)');
title('Late Noise');
mysavefig(h, filename, plotdir, fontsize, [9,2.2]);

% Interneuronal noise
sgmR = 0;
sgmG = 2;
sgmInput = 0;
subplot(1,3,3);
hold on;
rng(2025);
mylgd = [];
for i = 1:length(potentiation)
    eltp = potentiation(i);
    iltp = potentiation(i);
    balance_scale = iltp*eltp*mean(w*ones(N,1))*eqlb^2 + (iltp*BG+1-mean(a*ones(N,1)))*eqlb - BR;
    R0 = eqlb*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = eltp*w*R0+BG;
    initialvals = [R0'; G0'; D0'];
    Vprior = ones(size(cp));
    Vinput = ones(1,N).*cp;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, balance_scale, w, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);

    hst1 = histogram(R(5000:end,1),...
        'EdgeColor', 'none', 'FaceColor', colorpalette{i}, 'FaceAlpha', .3, 'Normalization', 'pdf');
    pd1 = fitdist(R(:,1),'kernel','Kernel','normal');
    x = hst1.BinEdges;
    y1 = pdf(pd1,x);
    plot(x, y1, '-', 'Color', colorpalette{i}, 'LineWidth', lwd);
end
% xlim([47,56]);
ylabel('Density');
xlabel('Activity (Hz)');
title('Interneuronal Noise');
mysavefig(h, filename, plotdir, fontsize, [9,2.2]);
%% Surrogated Fokker-Plank probability distribution
task = sprintf('PrbDstrbS_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 120; % second
stoprule = 1;
sgmG = 18;
myxrng = 10;
for j = 2
    switch j
        case 1
            sgmInput = 66;
            sgmR = 0;
            myylim = [0,60];
            myyrng = 30;
        case 2
            sgmInput = 0;
            sgmR = 26;
            myylim = [20,40];
            myyrng = 18;
    end
    
    for i = [1,2]
        h = figure;
        set(h, 'Units', 'inches','Position', [1,1,2.4,2.2]);
        filename = sprintf('%s%i%i',task, i, j);
        hold on;
        winput = winputs(i);
        wrg = wrgs(i);
        wgr = wgrs(i);
        R0 = eqlb*ones(N,1);
        D0 = 0*ones(N,1);
        G0 = wrg*R0+BG;
        initialvals = [R0'; G0'; D0'];
        Vprior = 100*ones(size(cp));
        Vinput = 100*ones(1,N).*cp;
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        data(:,1) = G(5000:end,1);
        data(:,2) = R(5000:end,1);
        mu = mean(data);
        Sigma = cov(data);
        % Grid for evaluation
        G_vals = linspace(mu(1)-3*sqrt(Sigma(1,1)), mu(1)+3*sqrt(Sigma(1,1)), 200);
        R_vals = linspace(mu(2)-3*sqrt(Sigma(2,2)), mu(2)+3*sqrt(Sigma(2,2)), 200);
        [G_grid, R_grid] = meshgrid(G_vals, R_vals);
        % Evaluate the 2D Gaussian PDF over the grid
        F = mvnpdf([G_grid(:) R_grid(:)], mu, Sigma);
        F = reshape(F, size(G_grid));

        % Plot the Gaussian density as shaded contour
        % contourf(G_grid, R_grid, max(F(:))-F, 20, 'LineColor', 'none'); % shadow effect
        surf(G_grid, R_grid, max(F(:))-F,'EdgeColor','none');
        colormap('gray');
        % Add contour lines for clarity
        % contour(G_grid, R_grid, max(F(:))-F, 8, 'LineColor', 'k', 'LineWidth', 1);
        % Scatter plot of the data points
        % scatter(G(5000:end,1), R(5000:end,1), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colorpalette{i}, 'MarkerFaceAlpha', .3);
        ylim([mu(2)-myyrng, mu(2) + myyrng]);
        xlim([mu(1)-myxrng, mu(1)+myxrng]);
        xlabel('G firing rate (Hz)');
        ylabel('R firing rate (Hz)');
        mysavefig(h, filename, plotdir, fontsize, [2.4,2.2]);
        exportgraphics(h, fullfile(plotdir, [filename, '.pdf']), 'ContentType','vector');
    end
end
%% Signal-to-Noise ratio over training 
task = sprintf('CovtCourse_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 120; % second
stoprule = 1;
sgmG = 18;
myxrng = 10;
for j = 1:2
    switch j
        case 1
            sgmInput = 66;
            sgmR = 0;
            myylim = [0,60];
            myyrng = 30;
        case 2
            sgmInput = 0;
            sgmR = 26;
            myylim = [20,40];
            myyrng = 18;
    end
    rate = 10000;
    x1 = 1:rate:(6*60/dt);
    x2 = (24*60/dt):rate:(25*60/dt);
    tvec = [x1, x2];
    % tvec = 1:10000:length(winputc);
    SNR = nan(length(tvec),2);
    SD = nan(length(tvec),3);
    r = nan(length(tvec), 1);
    for ti = 1:length(tvec)
        winput = winputc(tvec(ti));
        wrg = wrgc(tvec(ti));
        wgr = wgrc(tvec(ti));
        R0 = eqlb*ones(N,1);
        D0 = 0*ones(N,1);
        G0 = wrg*R0+BG;
        initialvals = [R0'; G0'; D0'];
        Vprior = 100*ones(size(cp));
        Vinput = 100*ones(1,N).*cp;
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
            sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        data(:,1) = G(5000:end,1);
        data(:,2) = R(5000:end,1);
        mu = mean(data);
        Sigma = cov(data);
        SD(ti,:) = [Sigma(1,1), Sigma(2,2), Sigma(1,2)];
        SNR(ti,:) = mu./[Sigma(1,1), Sigma(2,2)];
        r(ti) = Sigma(1,2)/sqrt(Sigma(1,1)*Sigma(2,2));
    end

    h = figure;
    set(h, 'Units', 'inches','Position', [1,1,2.4,2.2]);
    filename = sprintf('%s%i',task, j);
    hold on;
    ldg = [];
    yyaxis left;
    ldg(1) = plot(x1*dt,SD(1:length(x1),2), 'g-', 'LineWidth',lwd);
    ldg(2) = plot(x1*dt,SD(1:length(x1),1), 'r-', 'LineWidth',lwd);
    plot(50+(x2-min(x2)+max(x1))*dt,SD(length(x1)+1:end,2), 'g-', 'LineWidth',lwd);
    plot(50+(x2-min(x2)+max(x1))*dt,SD(length(x1)+1:end,1), 'r-', 'LineWidth',lwd);
    ylabel('S.D. of activities');
    ax = gca;
    ax.YColor = [0,0,0];
    yyaxis right;
    ldg(3) = plot(x1*dt,r(1:length(x1)), '--', 'Color', [.38,.38,.38], 'LineWidth',lwd);
    plot(50+(x2-min(x2)+max(x1))*dt,r(length(x1)+1:end), '--', 'Color', [.38,.38,.38], 'LineWidth',lwd);
    ylabel('Covariance');
    xticks([[0,2,4,6]*60, 24*60+(-min(x2)+max(x1))*dt+50,25*60+(-min(x2)+max(x1))*dt+50]);
    xticklabels({'0','2','4','6','24','25'});
    xlim([-20,1500+(-min(x2)+max(x1))*dt+50]);
    ax.YColor = [0,0,0];
    xlabel('Time (mins)');
    legend(ldg, {'R', 'G', 'Cov(R,G)'}, "Location","east");
    mysavefig(h, filename, plotdir, fontsize, [2.5,1.57]);
    % exportgraphics(h, fullfile(plotdir, [filename, '.pdf']), 'ContentType','vector');
end
