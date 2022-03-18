% Fig 2
%% panel b, Example dynamics of the LDDM
a = 15*eye(2,2);
b = 1.1*eye(2,2);
w = ones(2,2);
cp = .512;
initialvals = [10,10;20,20;0,0]*0;
predur = 0;
presentt = 0;
dur = 1.9;
stimdur = dur;
sgm = 0;
triggert = .9;
thresh = 70;
stoprule = 1;


h = figure; hold on;
filename = 'Fig2b';
Vprior = [0, 0] + B0;
Vinput = [1 + cp, 1 - cp]*scale0 + B0;
[choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
lgd4 = plot(G(:,2), '--', 'Color',"#659FBE", 'LineWidth',lwd);
lgd3 = plot(G(:,1), '-', 'Color',"#659FBE", 'LineWidth',lwd/2);
lgd6 = plot(I(:,2),'g--', 'Color','#06d6a0', 'LineWidth',lwd);
lgd5 = plot(I(:,1),'g-', 'Color','#06d6a0', 'LineWidth',lwd/2); % old color 669900
lgd2 = plot(R(:,2), 'k--', 'LineWidth',lwd);
lgd1 = plot(R(:,1), 'k-', 'LineWidth',lwd);
plot([predur + triggert + rt/2, predur + triggert + rt]/dt,  [thresh, thresh],'k-','LineWidth',.5);
ylim([-5,max(ylim)]);
yticks([0, 70]);
xlim([-50,(predur + dur)/dt]);
xticks([presentt/dt, triggert/dt]);
xticklabels({'stim on','action'});
yticklabels({'0','70'});
ylabel('Activity (a.u.)');
xlabel('Time (a.u.)');
legend([lgd1,lgd3,lgd5, lgd2,lgd4,lgd6],{'R_1','G_1','D_1', 'R_2','G_2','D_2'},...
    'Location','North','FontSize',fontsize-5, 'FontName','Times New Roman', ...
    'FontAngle','italic','NumColumns',2,'Box','off');
savefigs(h, filename, plotdir, fontsize, [3, 2.2]);