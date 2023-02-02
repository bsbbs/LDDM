function h = RT_Dstrbtn(dataBhvr, rtmat, choicemat, plot_dir, name)
lwd = 1;
fontsize = 11;
colorpalette = {'#ef476f','#ffd166','#06d6a0','#118ab2','#073b4c'};
aspect8 = [2, 6.4]; % for the long format RT distribution fitting panels
rate = length(rtmat)/1024;
maxrt = max(max(rtmat));
minrt = min(min(rtmat));
bank1 = [];
bank2 = [];
for ii = 1:6
    gap = (dataBhvr.rtrange(ii,2) - dataBhvr.rtrange(ii,1))/dataBhvr.numbins;
    BinEdge = [minrt:gap:(maxrt+gap)];
    hg = histogram(rtmat(choicemat(:,ii)==1,ii),'BinEdges',BinEdge);
    bank1{ii} = hg.Values/rate;
    hg = histogram(rtmat(choicemat(:,ii)==2,ii),'BinEdges',BinEdge);
    bank2{ii}= hg.Values/rate;
    BinMiddle{ii} = hg.BinEdges(1:end-1) + hg.BinWidth/2;
end

h = figure;
filename = sprintf('RTDistrb_%s.eps',name);
for ii = 1:6
    subplot(6,1,ii);hold on;
    bar(dataBhvr.bincenter(ii,1:30),dataBhvr.histmat(ii,1:30)*1024,'FaceColor',colorpalette{3},'EdgeAlpha',0);
    bar(dataBhvr.bincenter(ii,1:30),-dataBhvr.histmat(ii,31:60)*1024,'FaceColor',colorpalette{2},'EdgeAlpha',0,'EdgeColor','none');
    plot(BinMiddle{ii},bank1{ii},'Color',colorpalette{4},'LineWidth',lwd);
    plot(BinMiddle{ii},-bank2{ii},'Color',colorpalette{1},'LineWidth',lwd);
    if ii == 7
        legend({'','','Correct','Error'},'NumColumns',2,'Location','North');
        legend('boxoff');
    end
    ylim([-60,100]);
    yticks([-50:50:100]);
    yticklabels({'50','0','50','100'});
    xlim([100 1762]/1000);
    xticks([.5,1.0,1.5]);
    if ii == 6
        xticklabels({'.5','1.0','1.5'});
        xlabel('Reaction time (s)');
    else
        xticklabels({});
    end
    if ii == 1
        ylabel(' ');
    end
    set(gca, 'box','off');
    mysavefig(h, filename, plot_dir, fontsize, aspect8);
end
end