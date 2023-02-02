function h = PlotmeanRTACC(RoitmanDataDir, rtmat, choicemat, plot_dir, name)
lwd = 1;
mksz = 3;
fontsize = 11;

addpath(RoitmanDataDir);
ColumnNames608
load T1RT.mat;
x(:,R_RT) = x(:,R_RT)/1000;
cohlist = unique(x(:,R_COH));
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
end

Cohr = [0 32 64 128 256 512]/1000;
cplist = Cohr*100;
cplist(1) = 1.1;

acc = [];
meanrtc = [];
meanrtw = [];
for ii = 1:6
    acc(ii) = sum(choicemat(:,ii)==1)/(sum(choicemat(:,ii)==1) + sum(choicemat(:,ii)==2));
    meanrtc(ii) = mean(rtmat(choicemat(:,ii)==1,ii));
    meanrtw(ii) = mean(rtmat(choicemat(:,ii)==2,ii));
end

h = figure;
filename = sprintf('RT&ACC_%s',name);
subplot(2,1,1);
hold on;
plot(cplist, accr*100, 'xk', 'MarkerSize', mksz+1);
plot(cplist, acc*100,'-k','LineWidth',lwd);
ylim([.45,1]*100);
yticks([50,100]);
xlim([1,100]);
xticks([1,10,100]);
xticklabels({'0','10','100'});
ylabel('Correct (%)');
xlabel('Input coherence (%)');
set(gca, 'XScale', 'log');
legend({'data','model'},'NumColumns',1,'Location','SouthEast','FontSize',fontsize-2);
legend('boxoff');
mysavefig(h,filename,plot_dir,fontsize,[2,3.0]);

subplot(2,1,2);
hold on;
lg1 = plot(cplist, meanrtcr, '.k', 'MarkerSize', mksz*3);
lg2 = plot(cplist, meanrtc, '-k','LineWidth',lwd);
lg3 = plot(cplist, meanrtwr, 'ok', 'MarkerSize', mksz);
lg4 = plot(cplist, meanrtw, '--k','LineWidth',lwd);
xlim([1,100]);
xticks([1,10,100]);
xticklabels({'0','10','100'});
yticks([.4,1]);
ylim([.4, 1]);
ylabel('RT (secs)');
xlabel('Input coherence (%)');
set(gca, 'XScale', 'log');
mysavefig(h,filename,plot_dir,fontsize,[2,3.0]);
end