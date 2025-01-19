% Compare
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionDual';
load(fullfile(outdir, 'PPSync', 'SavedResult.mat'));
wEI_Sync = Ntwk.wEI;
InhbtMtrx_Sync = Ntwk.InhbtMtrx;
load(fullfile(outdir, 'PPAsync', 'SavedResult.mat'));
wEI_Async = Ntwk.wEI;
InhbtMtrx_Async = Ntwk.InhbtMtrx;


h = figure;
filename = sprintf('wEI_Sync_Async%1.0fmins', max(time));
imagesc((wEI_Sync - wEI_Async)');
caxis([-.25, .25]);
colormap(bluewhitered(256));
c = colorbar;
c.Location = 'northoutside';
c.Label.String = '\Deltaw_{R to G}';
xlabel('Inhbt channels');
ylabel('Exct channels');
mysavefig(h, filename, outdir, 14, [1.62, 3.6]);


h = figure;
filename = 'Meanfield-RcrtInhbt';
d = InhbtMtrx_Sync-InhbtMtrx_Async;
imagesc(d');
for i = 1:2
    for j = 1:2
        text(i-.3,j,sprintf('%.2e',d(i,j)), 'color','w');
    end
end
caxis([0, .007]);
colormap('gray');
c = colorbar;
c.Label.String = 'Recurrent Inhibition';
xticks([1,2]);
xticklabels({'Input 1', 'Input 2'});
yticks([1,2]);
yticklabels({'Input 1', 'Input 2'});
mysavefig(h, filename, outdir, 12, [3.4, 2.6]);

h = figure;
filename = 'Meanfield-RcrtInhbt_bars';
bar(InhbtMtrx_Sync-InhbtMtrx_Async);
xticks([1 2]) % Two groups of bars
xticklabels({'R1', 'R2'});
ylim([-.01,.01]);
xlabel('Origins');
ylabel('Inhibition');
legend('to R1', 'to R2');
mysavefig(h, filename, outdir, 12, [3.4, 2.6]);