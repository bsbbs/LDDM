function [h, h2, mRT] = PlotFitDynmcs(RoitmanDataDir, sm_mr1c, sm_mr2c, sm_mr1cD, sm_mr2cD, plot_dir, name, realdata)
load(fullfile(RoitmanDataDir,'DynmcsData.mat'));
aspect = [3, 2.5];
fontsize = 11;
lwd = 1;
mygray = flip(gray(8));
% colvec = flip({[218,166,109]/256,[155 110
% 139]/256,'#32716d','#af554d','#708d57','#3b5d64'}); % similar to
% Roitman paper
filename = sprintf('FittedDynmcs_%s',name);
yup = max([max(sm_mr1c), max(sm_mr2c), max(sm_mr1cD(1:find(sac_ax > 0, 1,"first"),:)), max(sm_mr2cD(1:find(sac_ax > 0, 1,"first"),:))]);
ylow = min([min(sm_mr1c), min(sm_mr2c), min(sm_mr1cD(1:find(sac_ax > 0, 1,"first"),:)), min(sm_mr2cD(1:find(sac_ax > 0, 1,"first"),:))]);
yrange = yup - ylow;
yuplim = yup + yrange/20;
ylowlim = ylow - yrange/20;
dot_tick = find(~isnan(sm_mr1c(:,6)), 1, "last");
dot_tick_x = dot_ax(dot_tick);
mRT = (dot_tick_x + 90)/1000;
h = figure;
subplot(1,2,1);hold on;
plot([dot_tick_x,dot_tick_x]/1000,[ylow,yup],'-k');
for ci = 1:6
    lg(ci) = plot(dot_ax/1000, sm_mr1c(:,ci),'Color',mygray(ci+2,:),'LineWidth',lwd);
    plot(dot_ax/1000, sm_mr2c(:,ci),'--','Color',mygray(ci+2,:),'LineWidth',lwd);
end
ylim([ylowlim, yuplim]);
ylabel('Firing rate (sp/s)');
xlabel('Time (secs)');
xlim([-.09, .8]);
xticks([0:.2:.8]-.09);
xticklabels({'0','.2','.4','.6','.8'});
mysavefig(h,filename,plot_dir,fontsize,aspect);

subplot(1,2,2);hold on;
plot([0.0,0.0],[ylow,yup],'-k');
for ci = 1:6
    lg(ci) = plot(sac_ax/1000, sm_mr1cD(:,ci),'Color',mygray(ci+2,:),'LineWidth',lwd);
    plot(sac_ax/1000, sm_mr2cD(:,ci),'--','Color',mygray(ci+2,:),'LineWidth',lwd);
end
xlim([-.8, .03]);
xticks([-.8:.2:0]+.03);
xticklabels({'-.8','-.6','-.4','-.2','0'});
ylim([ylowlim, yuplim]);
yticks([]);
set(gca,'ycolor',[1 1 1]);
legend(flip(lg),flip({'0','3.2','6.4','12.8','25.6','51.2'}),'Location','best','FontSize',fontsize-2);
mysavefig(h,filename,plot_dir,fontsize,aspect);
saveas(h,fullfile(plot_dir,[filename, '.fig']),'fig');
%% plot firing rates at position a,b,c,d
% colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
% rescale to the threshold
sac_tick = find(sac_ax == 0);
y = sm_mr1cD(sac_tick,:);
rate = mean(c)/mean(y);
%
Cohr = [0 32 64 128 256 512]/1000; % percent of coherence
x = Cohr*100;
h2 = figure;
filename = sprintf('abcd_%s',name);
subplot(2,1,1);hold on;
y = sm_mr1c(dot_tick,:)*rate;
if realdata
    plot(x, a,'bx', 'MarkerSize',8);
    sa = y;
end
scatter(x, y, 14, mygray(3:8,:), 'filled');
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-', 'LineWidth', lwd);
y = sm_mr2c(dot_tick,:)*rate;
if realdata
    plot(x, b, 'rx', 'MarkerSize',8);
    sb = y;
    RMSE = sqrt(sum(([a,b] - [sa, sb]).^2)/numel([a,b]));
    text(x(end-1), mean([a,b]), sprintf('RMSE = %1.2f', RMSE));
end
scatter(x, y, 14, mygray(3:8,:));
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k--', 'LineWidth',lwd);
xlim([-4,55.2]);
xticks([0:10:50]);
xticklabels({});
if realdata
    ylabel('Rescaled firing rates (sp/s)');
else
    ylabel('Firing rates (sp/s)');
end
mysavefig(h2, filename, plot_dir, fontsize, [2 3]);


subplot(2,1,2);hold on;
sac_tick = find(sac_ax == 0);
y = sm_mr1cD(sac_tick,:)*rate;
if realdata
    plot(x, c,'bx','MarkerSize',8);
    sc = y;
end
scatter(x, y, 14, mygray(3:8,:), 'filled');
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k-','LineWidth',lwd);
y = sm_mr2cD(sac_tick,:)*rate;
if realdata
    plot(x, d,'rx','MarkerSize',8);
    sd = y;
    RMSE = sqrt(sum(([c,d] - [sc, sd]).^2)/numel([c,d]));
    text(x(end-1), mean([c,d]), sprintf('RMSE = %1.2f', RMSE));
end
scatter(x, y, 14, mygray(3:8,:));
p = polyfit(x,y,1);
mdl = fitlm(x,y,'linear')
plot(x,p(1)*x+p(2),'k--','LineWidth',lwd);
xlim([-4,55.2]);
xticks([0:10:50]);
xlabel('Input strength (% coh)');
if realdata
    ylabel('Rescaled firing rates (sp/s)');
else
    ylabel('Firing rates (sp/s)');
end
mysavefig(h2, filename, plot_dir, fontsize, [2 3]);
end