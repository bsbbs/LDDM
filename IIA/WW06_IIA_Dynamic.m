%% panel a, dynamic of neural firing rates
outdir = './rslts/WW06_IIA_Dynamic';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;
aspect = [7, 6];
%% initalize parameters as in the paper of Wong & Wang, 2006
miu0 = 30; % Hz
sgm = .02; % nA
I0 = .3255; %nA
JN = [.2609, -.0497, -.0497;
    -.0497, .2609, -.0497;
    -.0497, -.0497, .2609]; %nA
gamma = .641;
tauS = .1; % second
tauAMPA = .002; % second
% simulation parameters
dur = 1.1; % second
dt = .001; % second
presentt = dt;
stimdur = dur;
thresh = 20; % Hz
initialvals = [2, 2, 2; .1, .1, .1]; % for H and S variables
stoprule = 1;
sims = 10240*2; % number of iterations

%c = [3.2 12.8, 25.6, 38.4 51.2]'/100;
c3 = [0:.2:1]';
c1 = 1+.256;
c2 = 1;
cp = [ones(size(c3))*c1, c2*ones(size(c3)), c3];
mygray = flip(gray(length(c3) + 1));
mycol2 = colormap(winter(length(c3)));
task = 'RT';

%% in space of R1-R2
sgm = 0;
simname = sprintf('WW06_DynamicR1R2_%s_c1%1.2f_c2%1.2f_%ic3_init%1.2f_sgm%2.1f',...
    task,c1,c2,length(c3),initialvals(1,1),sgm);
lim = 0;
h = figure; hold on;
filename = sprintf('%s',simname);
lgd3 = [];
for vi = 1:numel(c3)
    cpinput = cp(vi,:);
    [~, ~, R, ~, ~, ~] = wong06(cpinput,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
    lgd3(vi) = plot(R(:,1), R(:,2), '-', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
    plot([0,R(round((presentt + .6)/dt),1)],[0,R(round((presentt + .6)/dt),2)],'.-', 'Color', mycol2(vi,:), 'LineWidth',lwd/2,'MarkerSize',mksz*2);
    lim = max([lim, max(R(:,1:2))]);
end
if 1
    plot([0, 2*thresh],[thresh,thresh], 'k--');
    plot([thresh,thresh],[0, 2*thresh], 'k--');
end
plot([0,thresh],[0,thresh], 'k--');
ylim([0,lim*1.2]);
ylabel('R2 Activity (a.u.)');
xlim([0,lim*1.2]);
xlabel('R1 Activity (a.u.)');

lgd = legend(lgd3,cellstr(num2str(c3)),...
    'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
    'FontAngle','italic','NumColumns',1,'Box','off');
title(lgd, 'V_3');

savefigs(h, filename, outdir, fontsize, aspect);


%% in the space of time - R
sgm = 0;
simname = sprintf('WW06_Dynamic_%s_c1%1.2f_c2%1.2f_%ic3_init%1.2f_sgm%2.1f',...
    task,c1,c2,length(c3),initialvals(1,1),sgm);

h = figure; hold on;
filename = sprintf('%s',simname);
clear lgd3;
for vi = 1:numel(c3)
    cpinput = cp(vi,:);
    [~, ~, R, ~, ~, ~] = wong06(cpinput,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
    lgd3(vi) = plot(R(:,3), 'k-.', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
    lgd2(vi) = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
    lgd1(vi) = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
end

plot([1.2, dur]/dt,[thresh,thresh], 'k-');
% text(1200,thresh*1.1,'threshold');
yticks([0]);
yticklabels({'0'});
ylabel('Activity (a.u.)');
xticks([]);
xlabel('Time (s)');

lgd = legend(lgd3,cellstr(num2str(c3)),...
    'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
    'FontAngle','italic','NumColumns',1,'Box','off');
title(lgd, 'V_3');

savefigs(h, filename, outdir, fontsize, aspect);


