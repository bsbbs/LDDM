function H = Vslz_ParamSpc(x,y,namex,namey,filename,ACC,meanRT,meandR,thresh,modelname,vi)
dmeanRT = squeeze(meanRT(vi,:,:,1) - meanRT(vi,:,:,2));
dACC = squeeze(ACC(vi,:,:,1) - ACC(vi,:,:,2));
[xmat, ymat] = meshgrid(x,y);
switch modelname
    case 'LcD_GABA'
        % xrng = [0,55];
        % yrng = [0,4];
        xrng = [min(x),max(x)];
        yrng = [min(y),4.7];
        tmp =  -max(abs(dACC(:,y <= yrng(2))'));
        dACC(:,end) = tmp; % set a max negative value outside of the visible space
        dACC(:,end-1) = -tmp; % set a max positive value outside of the visible space
        dACC(abs(dACC) > max(abs(tmp))) = max(abs(tmp)); % flatten of all spikes of extreme values
        dmeanRT(:,end) = -max(abs(dmeanRT(:,y <= yrng(2))')); % set a max negative value outside of the visible space
        dmeanRT(:,end-1) = max(abs(dmeanRT(:,y <= yrng(2))')); % set a max positive value outside of the visible space
    case 'LcDFD_GABA'
        xrng = [min(x),max(x)];
        yrng = [min(y),max(y)];
    case 'XJ_GABA'
        xrng = [min(x),max(x)];
        yrng = [min(y),max(y)];
    case 'AymW_GABA'
        xrng = [min(x),max(x)];
        yrng = [min(y),max(y)];
    case 'AymWFD_GABA'
        xrng = [min(x),max(x)];
        yrng = [min(y),max(y)];
end

h=figure;   colormap(jet);
for GABAi = 1:2
    subplot(1,2,GABAi);
    % imagesc(x,y,squeeze(ACC(vi,:,:,GABAi)));
    % set(gca,'YDir','normal');
    s = surf(xmat,ymat,squeeze(ACC(vi,:,:,GABAi))'*100);
    s.EdgeColor = 'none';
    view([0,90]);
    xlim([min(x),max(x)]);ylim(yrng);
    xlabel(namex);ylabel(namey);
    c = colorbar;
    ylabel(c, 'Accuracy [%]');
    set(gca,'FontSize',18);
end
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 12 4];
saveas(h,fullfile('./graphics', modelname, sprintf('ACC_V%i_%s.eps',vi,filename)),'epsc2');
H = h;

h=figure;   colormap(jet);
for GABAi = 1:2
    subplot(1,2,GABAi);
    %imagesc(x,y,squeeze(meanRT(vi,:,:,GABAi)));
    %set(gca,'YDir','normal');
    s = surf(xmat,ymat,squeeze(meanRT(vi,:,:,GABAi))');
    s.EdgeColor = 'none';
    view([0,90]);
    xlim([min(x),max(x)]);ylim(yrng);
    xlabel(namex);ylabel(namey);
    c = colorbar;
    set(gca, 'ColorScale','log');
    ylabel(c, 'RT [s]');
    set(gca,'FontSize',18);
end
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 12 4];
saveas(h,fullfile('./graphics',modelname, sprintf('RT_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];

%% Difference of RT and ACC, whole area(unmasked)
h=figure;   %colormap(jet);
% subplot(1,2,1);
% imagesc(x,y,squeeze(meanRT(vi,:,:,2) - meanRT(vi,:,:,1)));
% set(gca,'YDir','normal');
% xmask = x>=xrng(1) & x<= xrng(2);
% ymask = y>=yrng(1) & y<= yrng(2);
% maskedx = x(xmask);
% maskedy = y(ymask);
% [maskedxmat, maskedymat] = meshgrid(maskedx,maskedy);
s = surf(xmat,ymat,dmeanRT');
s.EdgeColor = 'none';
view([0,90]);
xlim(xrng);ylim(yrng);
% xlim([min(x),max(x)]);ylim([min(y),max(y)])
xlabel(namex);ylabel(namey);
colormap(bluewhitered);
c = colorbar;
ylabel(c, '\DeltaRT [s]');
%freezeColors;
set(gca,'FontSize',18);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile('./graphics',modelname, sprintf('DeltaRT_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];

h=figure; % subplot(1,2,2);
% imagesc(x,y,squeeze(ACC(vi,:,:,2) - ACC(vi,:,:,1)));
% set(gca,'YDir','normal');
s = surf(xmat,ymat,dACC'*100);
s.EdgeColor = 'none';
view([0,90]);
xlim(xrng);ylim(yrng);
% xlim([min(x),max(x)]);ylim([min(y),max(y)])
xlabel(namex);ylabel(namey);
colormap(bluewhitered);
c = colorbar;
ylabel(c, '\DeltaAccuracy [%]');
set(gca,'FontSize',18);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile('./graphics',modelname, sprintf('DeltaACC_V%i_%s.eps',vi,filename)),'epsc2');
% saveas(h,fullfile(outdir, sprintf('Delta_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];

%% difference of R1-R2 at the time point of chioces were made
h=figure;   colormap(jet);
for Gabai = 1:2
    subplot(1,2,Gabai);
    s = surf(xmat,ymat,squeeze(meandR(vi,:,:,Gabai))');
    s.EdgeColor = 'none';
    view([0,90]);
    xlim([min(x),max(x)]);ylim([min(y),max(y)])
    xlabel(namex);ylabel(namey);
    c = colorbar;
    ylabel(c, '\DeltaR [Hz]');
    set(gca,'FontSize',18);
end
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 12 4];
saveas(h,fullfile('./graphics',modelname, sprintf('dR_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];

%% Difference of RT and ACC, masked by the rule of R1-R2 bifurcation
switch modelname
    case '.\graphics\LcD_GABA'
        xrng = [0,55];
        yrng = [0,4];
    case '.\graphics\XJ_GABA'
        xrng = [0,.7];
        yrng = [min(y),max(y)];
    case '.\graphics\AymW_GABA'
        xrng = [min(x),max(x)];
        yrng = [min(y),max(y)];
end
mask = nan(size(xmat));
dR1 = squeeze(meandR(vi,:,:,1))' > thresh/3;
dR2 = squeeze(meandR(vi,:,:,2))' > thresh/3;
mask(dR1 & dR2) = 1;
h=figure;   %colormap(jet);
% subplot(1,2,1);
% imagesc(x,y,squeeze(meanRT(vi,:,:,2) - meanRT(vi,:,:,1)));
% set(gca,'YDir','normal');
s = surf(xmat,ymat,dmeanRT'.*mask);
s.EdgeColor = 'none';
view([0,90]);
xlim(xrng);ylim(yrng);
xlabel(namex);ylabel(namey);
colormap(bluewhitered);
% freezeColors;
c = colorbar;
ylabel(c, '\DeltaRT [s]');
set(gca,'FontSize',18);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile('./graphics',modelname, sprintf('MaskedDeltaRT_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];

h=figure; %subplot(1,2,2);
% imagesc(x,y,squeeze(ACC(vi,:,:,2) - ACC(vi,:,:,1)));
% set(gca,'YDir','normal');
s = surf(xmat,ymat,dACC'.*mask*100);
s.EdgeColor = 'none';
view([0,90]);
xlim(xrng);ylim(yrng)
xlabel(namex);ylabel(namey);
colormap(bluewhitered);
c = colorbar;
ylabel(c, '\DeltaAccuracy [%]');
set(gca,'FontSize',18);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 5.3 4];
saveas(h,fullfile('./graphics',modelname, sprintf('MaskedDeltaACC_V%i_%s.eps',vi,filename)),'epsc2');
% saveas(h,fullfile(outdir, sprintf('MaskedDelta_V%i_%s.eps',vi,filename)),'epsc2');
H = [H, h];