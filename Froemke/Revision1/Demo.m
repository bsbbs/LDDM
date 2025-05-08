% Demo

%% define paths
Homedir = 'C:\Users\Bo\Documents\GitHub';
% Homedir = '~';
% addpath(fullfile(Homedir,'Documents','LDDM','CoreFunctions'));
addpath(fullfile(Homedir,'LDDM','Froemke','Revision1'));
addpath(fullfile(Homedir,'LDDM','Froemke','utils'));
Drpbx = "C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project0";
out_dir = fullfile(Drpbx, 'Revision1');
if ~exist("out_dir",'dir')
    mkdir(out_dir);
end
plotdir = fullfile(out_dir,'Graphics');
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Simdir = fullfile(out_dir,'SimRslts');
if ~exist(Simdir,'dir')
    mkdir(Simdir);
end
%% parameters for visulization
fontsize = 14;
mksz = 25;
lwd = 1.5;
cp = [.032, .064,.128, .256];
mygray = flip(gray(numel(cp) + 2));
myred = mygray; myred(:,1) = 1;
myblue = mygray; myblue(:,3) = 1;
%% parameters for simulation
dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauD = .1; % second
Tau = [tauR, tauG, tauD];
w0 = 1;
a0 = 0;
b0 = 1.5; % 1.1;
c_rprsnt = .064;
sgmInput_rprsnt = 1/3;
c_choice = .032;
sgmInput_choice = .75;
predur = 0;
presentt = dt;
thresh = 70; % Hz
stoprule = 1;
eqlb = 32;
BG = 30;
BR = 30;

%% Example dynamics
%% Representation dynamic
task = 'VR_RT';
a = a0*eye(2);
b = b0*eye(2);
stimdur = Inf;
triggert = Inf;
dur = 1.8; % second
% Dynamics
sgm = .01;
sgmInput = sgmInput_rprsnt/5;
potentiation = linspace(0.1,1,5); %[1:.5:2.5];
h = figure;
subplot(2,1,1);
hold on;
filename = sprintf('timeCourse_%s_from%1.1fto%1.1f_sgm%2.2f_sgmInput%.2f',task, min(potentiation), max(potentiation), sgm, sgmInput);
rng(5);
scale0 = w*[1;1]*eqlb^2 + (BG+1-a*[1;1])*eqlb - BR;
for i = 1:length(potentiation)
    eltp = potentiation(i);
    iltp = potentiation(i);

    w = w0*ones(2);
    a = a0*eye(2);
    b = b0*eye(2);

    scale = iltp*eltp*w*[1;1]*eqlb^2 + (iltp*BG+1-a*[1;1])*eqlb - BR; %2*w0*eqlb.^2 + (1-a0).*eqlb;
    R0 = eqlb*[1;1];
    D0 = 0*[1;1];
    G0 = eltp*w*R0+BG;
    initialvals = [R0'; G0'; D0'];

    STDP_GR = iltp;
    Vinput = [1 + c_rprsnt, 1 - c_rprsnt].*scale';
    Vprior = ones(size(Vinput)).*scale';

    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, w, a, b,...
        sgm, sgmInput*mean(scale0), Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    x = (1:numel(R(:,1)))*dt;
    plot(x, R(:,1), '-', 'LineWidth', lwd, 'Color', myred(i+1,:));
    plot(x, R(:,2), '-', 'LineWidth', lwd, 'Color', myblue(i+1,:));
    xlim([-50*dt, dur]);
    % plot(x, G(:,1), 'r--', 'LineWidth', lwd);
    % plot(x, G(:,2), 'b--', 'LineWidth', lwd);
    % plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
    % plot(x, D(:,2), 'b-.', 'LineWidth', lwd);
end
%yticks([10, 15, 20, 25]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
legend({'R_1','R_2'},'Box','off','Location','east',...
    'FontName','Arial', 'FontAngle','italic',...
    'FontSize',fontsize-4);
title('Representation');

mysavefig(h, filename, plotdir, fontsize, [2.9,4]);

% Reaction-time task
stimdur = Inf;
triggert = presentt;
dur = 1.8; % second
rng(8);
subplot(2,1,2); hold on;
for i = 1:length(potentiation)
    eltp = potentiation(i);
    iltp = potentiation(i);

    w = w0*ones(2);
    a = a0*eye(2);
    b = b0*eye(2);

    scale = iltp*eltp*w*[1;1]*eqlb^2 + (iltp*BG+1-a*[1;1])*eqlb - BR; %2*w0*eqlb.^2 + (1-a0).*eqlb;
    R0 = eqlb*[1;1];
    D0 = 0*[1;1];
    G0 = w*eltp*R0+BG;
    Vprior = scale';
    initialvals = [R0'; G0'; D0'];

    Vprior = ones(size(Vinput)).*scale';
    rtvec = [];
    Vinput = [1 + c_rprsnt, 1 - c_rprsnt].*scale';
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv1(Vprior, Vinput, BR, BG, eltp, iltp, w, a, b,...
        sgm, sgmInput*mean(scale0), Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    rtvec(vi) = rt;
    x = (1:numel(R(:,1)))*dt;
    plot(x, R(:,1), '-', 'LineWidth', lwd, 'Color', myred(i+1,:));
    plot(x, R(:,2), '-', 'LineWidth', lwd, 'Color', myblue(i+1,:));
    % plot(x, G(:,1), 'r--', 'LineWidth', lwd);
    % plot(x, G(:,2), 'b--', 'LineWidth', lwd);
    % plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
    % plot(x, D(:,2), 'b-.', 'LineWidth', lwd);

    % plot(rtvec, thresh*ones(size(rtvec)),'k-','LineWidth',.5);
    % text(mean(rtvec), thresh*1.05,'Threshold','FontName','Arial');
end
xlim([-50*dt, dur]);
%ylim([10,thresh]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
mysavefig(h, filename, plotdir, fontsize, [2.9,4]);
