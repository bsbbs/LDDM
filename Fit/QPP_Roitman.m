function [h] = QPP_Roitman(dataBhvr, rtmat, choicemat, plot_dir, name)
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
lwd = 1.0;
mksz = 3;
fontsize = 11;
x = dataBhvr.proportionmat;
y = dataBhvr.q;
qntls = dataBhvr.qntls;
h = figure; hold on;
acc = [];
for vi = 1:length(x)
    xc = x(vi)*ones(size(y(:,1,vi)));
    xw = 1 - x(vi)*ones(size(y(:,2,vi)));
    lgc = plot(xc,y(:,1,vi),'x','Color',colorpalette{4}, 'MarkerSize',mksz+1,'LineWidth',lwd);
    lge = plot(xw,y(:,2,vi),'x','Color', colorpalette{1}, 'MarkerSize',mksz+1,'LineWidth',lwd);
    % fitted value
    En(vi) = numel(rtmat(:,vi));
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    acc(vi) = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
    q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
end
for qi = 1:size(q,1)
    xq = [flip(1-acc), acc]';
    plot(xq,[squeeze(flip(q(qi,2,:)));squeeze(q(qi,1,:))],'k-o','MarkerSize',mksz,'LineWidth',lwd/2);
end
legend([lge,lgc],{'error','correct'},"NumColumns",2,'Location','northeast','FontSize',fontsize-2);
legend('box','off');
xlim([-.05 1.05]);
ylim([.2, 1.4]);
yticks([.2:.4:1.4]);
xlabel('Proportion');
ylabel('RT (s)');
filename = sprintf('QPPlot_%s',name);
mysavefig(h, filename, plot_dir, fontsize, [2.5 2.5]);
end