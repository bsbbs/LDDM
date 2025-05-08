%% Singleton - representation
%% Common parameters
N = 1;
w = w0*ones(N);
a = a0*eye(N);
b = b0*eye(N);
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
cp = 1;
%% Temporal dynamics
task = sprintf('Rprsnt_Dynamic_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 60; % second
stoprule = 1;
balance_scale0 = w*ones(N,1)*eqlb^2 + (BG+1-a*ones(N,1))*eqlb - BR;
sgmR = 0;
sgmG = 0;
sgmInput = sgmInput_rprsnt*mean(balance_scale0);
potentiation = linspace(0,1,5);

h = figure;
hold on;
filename = sprintf('%s_from%1.1fto%1.1f',task, min(potentiation), max(potentiation));
% early noise
rng(2025);
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
    x = (1:numel(R(:,1)))*dt;
    plot(x, R(:,1), '-', 'LineWidth', lwd, 'Color', myred(i+1,:));
    xlim([-50*dt, dur]);
end
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [2.9,2.2]);
%% Distribution of representation
task = sprintf('Rprsnt_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 60; % second
stoprule = 1;
balance_scale0 = mean(w*ones(N,1))*eqlb^2 + (BG+1-mean(a*ones(N,1)))*eqlb - BR;
sgmR = 0;
sgmG = 0;
sgmInput = sgmInput_rprsnt*mean(balance_scale0);
potentiation = linspace(0,1,5);

h = figure;
subplot(1,3,1);
hold on;
filename = sprintf('%s_from%1.1fto%1.1f',task, min(potentiation), max(potentiation));
% early noise
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
    mylgd(i) = plot(x,y1, '-', 'Color', colorpalette{i}, 'LineWidth', lwd);
end
ylabel('Density');
xlabel('Activity (Hz)');
lgd = legend(mylgd, {'0','.25','.50','.75','1'},...
    'Box','off','Location','east',...
    'FontName','Arial','FontSize',fontsize-4);
title(lgd, 'LTP');
title('Early Noise');
mysavefig(h, filename, plotdir, fontsize, [9,2.2]);

% Late noise
sgmR = 2;
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
%% Rate-based Plasticity 


%% Example trace
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 1; % second
task = sprintf('Trace_N%i',N);
h = figure;
hold on;
filename = sprintf('%s_from%1.1fto%1.1f',task, min(potentiation), max(potentiation));
sgmR = 0;
sgmG = 0;
sgmInput = 0;
p = 0;
for i = [1,2,length(potentiation)]
    eltp = potentiation(i);
    iltp = potentiation(i);
    balance_scale = iltp*eltp*mean(w*ones(N,1))*eqlb^2 + (iltp*BG+1-mean(a*ones(N,1)))*eqlb - BR;
    R0 = 5*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = w*R0;
    initialvals = [R0'; G0'; D0']; %zeros(3,N);
    Vprior = zeros(size(cp));
    Vinput = ones(1,N).*cp;
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, balance_scale, w, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    p = p + 1;
    subplot(2,3,p);hold on;
    x = (1:numel(R(:,1)))*dt;
    plot(x, G(:,1), '--', 'LineWidth', 2, 'Color', colorpalette{i});
    plot(x, R(:,1), '-', 'LineWidth', 2, 'Color', colorpalette{i});
    xlim([-50*dt, dur]);
    ylim([0, max([R; G])*1.2]);
    ylabel('Firing rates (Hz)');
    xlabel('Time (s)');
    mysavefig(h, filename, plotdir, fontsize, [8.2,2.2*2]);
    plot(G(:,1), R(:,1), '-', 'Color', colorpalette{i}, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('R firing rate (Hz)');
    subplot(2,3,3+p);hold on;
    plot(G(:,1), R(:,1), '-', 'Color', colorpalette{i}, 'LineWidth', 2);
    xlim([0, max(G)*1.4]);
    ylim([0, max(R)*1.4]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    mysavefig(h, filename, plotdir, fontsize, [8.2,2.2*2]);
end
%% Vector field
task = sprintf('Vctrfld_N%i',N);
filename = sprintf('%s_from%1.1fto%1.1f',task, min(potentiation), max(potentiation));
BR = 30;
BG = 30;
R0 = 3:7:70;
G0 = 3:7:100;
[Rm, Gm] = meshgrid(R0, G0);
h = figure;
p = 0;
for i = [1,2,length(potentiation)]
    p = p + 1;
    subplot(1,3,p);hold on;
    eltp = potentiation(i);
    iltp = potentiation(i);
    balance_scale = iltp*eltp*mean(w*ones(N,1))*eqlb^2 + (iltp*BG+1-mean(a*ones(N,1)))*eqlb - BR;
    Vinput = ones(1,N).*cp;
    dRm = -Rm + (Vinput*balance_scale + a.*Rm + BR)./(1 + Gm*iltp);
    dGm = -Gm + eltp*w.*Rm + BG;
    rate = 0.15;
    dRm = dRm.*rate;
    dGm = dGm.*rate;
    quiver(Gm,Rm,dGm,dRm,'off','Color',colorpalette{i},'LineWidth',1);
    xlim([0, max(G0)]);
    ylim([0, max(R0)]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    mysavefig(h, filename, plotdir, fontsize, [8.2,2.2], 0);
end
%% Simulated Fokker-Plank probability distribution
task = sprintf('PrbDstrb_N%i',N);
predur = 0;
presentt = 0;
stimdur = Inf;
triggert = Inf;
dur = 120; % second
stoprule = 1;
h = figure;
filename = sprintf('%s_from%1.1fto%1.1f',task, min(potentiation), max(potentiation));
balance_scale0 = mean(w*ones(N,1))*eqlb^2 + (BG+1-mean(a*ones(N,1)))*eqlb - BR;
sgmInput = sgmInput_rprsnt*mean(balance_scale0);
sgmR = 0;
sgmG = 12;
p = 0;
for i = [1,2,length(potentiation)]
    p = p + 1;
    subplot(1,3,p);hold on;
    eltp = potentiation(i);
    iltp = potentiation(i);
    balance_scale = iltp*eltp*mean(w*ones(N,1))*eqlb^2 + (iltp*BG+1-mean(a*ones(N,1)))*eqlb - BR;
    R0 = eqlb*ones(N,1);
    D0 = 0*ones(N,1);
    G0 = eltp*w*R0+BG;
    initialvals = [R0'; G0'; D0'];
    Vprior = ones(size(cp));
    Vinput = (ones(1,N).*cp);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, balance_scale, w, a, b,...
        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    data(:,1) = G(5000:end,1);
    data(:,2) = R(5000:end,1);
    mu = mean(data);
    Sigma = cov(data);
    % Grid for evaluation
    G_vals = linspace(mu(1)-4*sqrt(Sigma(1,1)), mu(1)+4*sqrt(Sigma(1,1)), 200);
    R_vals = linspace(mu(2)-4*sqrt(Sigma(2,2)), mu(2)+4*sqrt(Sigma(2,2)), 200);
    [G_grid, R_grid] = meshgrid(G_vals, R_vals);
    % Evaluate the 2D Gaussian PDF over the grid
    F = mvnpdf([G_grid(:) R_grid(:)], mu, Sigma);
    F = reshape(F, size(G_grid));

    % Plot the Gaussian density as shaded contour
    % contourf(G_grid, R_grid, max(F(:))-F, 20, 'LineColor', 'none'); % shadow effect
    surf(G_grid, R_grid, max(F(:))-F,'EdgeColor','none');
    colormap('gray');
    % Add contour lines for clarity
    % contour(R_grid, G_grid, F, 5, 'LineColor', 'k', 'LineWidth', 1);
    % Scatter plot of the data points
    % scatter(G(5000:end,1), R(5000:end,1), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colorpalette{i}, 'MarkerFaceAlpha', .3);
    ylim([5, 60]);
    % xlim([5, 60]);
    xlabel('G firing rate (Hz)');
    ylabel('R firing rate (Hz)');
    mysavefig(h, filename, plotdir, fontsize, [8.2,2.2]);

    %     hst1 = histogram(R(5000:end,1),...
    %         'EdgeColor', 'none', 'FaceColor', colorpalette{i}, 'FaceAlpha', .3, 'Normalization', 'pdf');
    %     pd1 = fitdist(R(:,1),'kernel','Kernel','normal');
    %     x = hst1.BinEdges;
    %     y1 = pdf(pd1,x);
    %     plot(x,y1, '-', 'Color', colorpalette{i}, 'LineWidth', lwd);
end