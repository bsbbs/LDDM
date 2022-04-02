function savefigs(h, filename, outdir, fontsize, aspect)
set(gca,'FontSize',fontsize);
set(gca,'FontName','Arial')
set(gca,'TickDir','in');
set(gca,'LineWidth',1); 
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', fontsize-2)
set(xl, 'FontSize', fontsize);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', fontsize-2)
set(yl, 'FontSize', fontsize);
axis tight;
% 
% xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
% yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 aspect];
saveas(h,fullfile(outdir,sprintf('%s.eps',filename)),'epsc2');