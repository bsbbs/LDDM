%% visualization
if show
    h = figure;
    filename = sprintf('wEE_%1.0fmins', max(time));
    imagesc(Ntwk.wEE');
    caxis([0, .1]);
    % caxis([0, 1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'w_{R to R}';
    c.Location = 'northoutside';
    xlabel('Exct channels');
    ylabel('Exct channels');
    mysavefig(h, filename, testdir, 14, [3.2, 3.6]);

    h = figure;
    filename = sprintf('wEI_%1.0fmins', max(time));
    imagesc(Ntwk.wEI');
    caxis([0, 1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'w_{R to G}';
    c.Location = 'northoutside';
    xlabel('Inhbt channels');
    ylabel('Exct channels');
    mysavefig(h, filename, testdir, 14, [1.6, 3.6]);

    h = figure;
    filename = sprintf('wIE_%1.0fmins', max(time));
    imagesc(Ntwk.wIE');
    caxis([0, .2]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'w_{G to R}';
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, testdir, 14, [3.92, 1.46]);

    % synaptic weight change
    h = figure;
    filename = 'wEE_change';
    tmp = Ntwk.wEE - Ntwk.wEE_initial;
    imagesc(tmp');
    scale = max(abs(min(tmp(:))), max(tmp(:)));
    %caxis([-scale, scale]);
    caxis([-.1, .1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = '\Deltaw_{R to R}';
    c.Location = 'northoutside';
    xlabel('Exct channels');
    ylabel('Exct channels');
    mysavefig(h, filename, testdir, 14, [3.2, 3.6]);

    h = figure;
    filename = 'wEI_change';
    tmp = Ntwk.wEI - Ntwk.wEI_initial;
    imagesc(tmp');
    scale = ceil((max(abs(min(tmp(:))), max(tmp(:))))*10)/10;
    caxis([-scale, scale]);
    caxis([-.41, .41]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Location = 'northoutside';
    c.Label.String = '\Deltaw_{R to G}';
    xlabel('Inhbt channels');
    ylabel('Exct channels');
    mysavefig(h, filename, testdir, 14, [1.62, 3.6]);

    % synaptic weight change
    h = figure;
    filename = 'wIE_change';
    tmp = Ntwk.wIE - Ntwk.wIE_initial;
    imagesc(tmp');
    scale = ceil((max(abs(min(tmp(:))), max(tmp(:))))*10)/10;
    caxis([-scale, scale]);
    %caxis([-.05, .05]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = '\Deltaw_{G to R}';
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, testdir, 14, [3.92, 1.46]);

    % dynamics
    smpl_time = 5*60+[10000:30000]*dt; % unit in Secs
    segment = round(smpl_time/dt);
    % R
    h = figure;
    filename = sprintf('R_Dynamic%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    hold on;
    interval = 4;
    Amp = max(max(Rsection))/interval;
    for i = 1:interval:Ntwk.Exct.N
        plot(smpl_time, Rsection(i,:)+i*Amp, 'k-');
    end
    ylim([0, (i+interval)*Amp]);
    yticks([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]'));
    ylabel('Exct channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, testdir, 14, [3, 11.2]);
    % G
    h = figure;
    filename = sprintf('G_Dynamic%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    hold on;
    Amp = max(max(Gsection))/interval;
    for i = 1:interval:Ntwk.Inhbt.N
        plot(smpl_time, Gsection(i,:)+i*Amp, 'k-');
    end
    ylim([0, (i+interval)*Amp]);
    yticks([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]'));
    ylabel('Inhbt channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, testdir, 14, [3, 3.5]);
    
    % Example dynamics R and G
    h = figure;
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fs_%1.0fs',IdxE(1), IdxI(1), min(smpl_time), max(smpl_time));
    hold on;
    segment = round(smpl_time/dt);
    plot(smpl_time, Gsample(segment)/max(Gsample(segment)), '-', 'Color', "#659FBE", 'LineWidth',.5);
    plot(smpl_time, Rsample(segment)/max(Rsample(segment)), 'k-', 'LineWidth',.5);
    xlabel('Time (secs)');
    ylabel('Normalized activity');
    %ylabel('Firing rates (Hz)');
    mysavefig(h, filename, testdir, 14, [1.8, 1.2], 2);

    h = figure;
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fmins',IdxE(1), IdxI(1), max(time));
    hold on;
    lg2 = plot(time, Gsample, '-', 'Color', "#659FBE", 'LineWidth', .5);
    lg1 = plot(time, Rsample, 'k-', 'LineWidth', .5);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [47, 47, .1, .1, 47], '-k');
    xlim([0, 35]);
    %ylim([0, 12]);
    ylabel('Firing rates (Hz)');
    xlabel('Time (mins)');
    legend([lg1, lg2], {"R", "G"}, 'Location', "best", 'FontSize', 10);
    mysavefig(h, filename, testdir, 14, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wEE_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, EEdynamc(segment), 'k-', 'LineWidth',1);
    ylim([0,.02]);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [1.8, 1.2], 2);

    h = figure; hold on;
    filename = 'wEE_dynamic';
    plot(time, EEdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [.011, .011, .009, .009, 0.011], '-k');
    xlim([0, 35]);
    ylim([0,.02]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wEI_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, EIdynamc(segment), 'k-', 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [1.8, 1.2], 2);

    h = figure; hold on;
    filename = 'wEI_dynamic';
    plot(time, EIdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [.058, 0.058, .038, .038, .058]-.006, '-k');
    xlim([0, 35]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wIE_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, IEdynamc(segment), 'k-', 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [1.8, 1.2], 2);
    
    h = figure; hold on;
    filename = 'wIE_dynamic';
    plot(time, IEdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, ([.063, .063, .051, .051, .063]+.035)/10, '-k');
    xlim([0, 35]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, testdir, 12, [3, 2], 2);

    %% recurrent inhibition
    % tune, input connection weight .* EI and IE weights
    h = figure;
    filename = 'Tuning_Exct';
    imagesc(Ntwk.Exct.tuning');
    caxis([0, 1]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Exct channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Tuning';
    mysavefig(h, filename, testdir, 12, [3, 1]);
    
    h = figure;
    filename = 'Tuning_Inhbt';
    Ntwk.Inhbt.tuning = Ntwk.wEI * Ntwk.Exct.tuning;
    imagesc(Ntwk.Inhbt.tuning');
    caxis([0, 2.5]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Inhbt channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Tuning';
    mysavefig(h, filename, testdir, 12, [3, 1]);
    
    h = figure;
    filename = 'Tuning_Inhbt_IE';
    Ntwk.Inhbt.tuningIE = Ntwk.wIE' * Ntwk.Exct.tuning;
    imagesc(Ntwk.Inhbt.tuningIE');
    caxis([0, 1.4]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Inhbt channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Selective Inhibition';
    mysavefig(h, filename, testdir, 12, [3, 1]);
    % mean-field connectome
    InhbtMtrx(1,1) = sum((Ntwk.Inhbt.tuning(:,1)/sum(Ntwk.Inhbt.tuning(:,1), 1)).*(Ntwk.Inhbt.tuningIE(:,1)/sum(Ntwk.Inhbt.tuningIE(:,1),1)), 1);
    InhbtMtrx(1,2) = sum((Ntwk.Inhbt.tuning(:,1)/sum(Ntwk.Inhbt.tuning(:,1), 1)).*(Ntwk.Inhbt.tuningIE(:,2)/sum(Ntwk.Inhbt.tuningIE(:,2),1)), 1);
    InhbtMtrx(2,1) = sum((Ntwk.Inhbt.tuning(:,2)/sum(Ntwk.Inhbt.tuning(:,2), 1)).*(Ntwk.Inhbt.tuningIE(:,1)/sum(Ntwk.Inhbt.tuningIE(:,1), 1)), 1);
    InhbtMtrx(2,2) = sum((Ntwk.Inhbt.tuning(:,2)/sum(Ntwk.Inhbt.tuning(:,2), 1)) .*(Ntwk.Inhbt.tuningIE(:,2)/sum(Ntwk.Inhbt.tuningIE(:,2), 1)), 1);
    Ntwk.InhbtMtrx = InhbtMtrx;
    h = figure;
    filename = 'Meanfield-RcrtInhbt';
    imagesc(Ntwk.InhbtMtrx');
    for i = 1:2
        for j = 1:2
            text(i-.3,j,sprintf('%.3f',Ntwk.InhbtMtrx(i,j)));
        end
    end
    caxis([0, .03]);
    colormap('gray');
    %colormap(bluewhitered(556));
    c = colorbar;
    c.Label.String = 'Recurrent Inhibition';
    xticks([1,2]);
    xticklabels({'Input 1', 'Input 2'});
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    mysavefig(h, filename, testdir, 12, [3.4, 2.6]);

    h = figure;
    filename = 'Meanfield-RcrtInhbt_bars';
    bar(Ntwk.InhbtMtrx);
    xticks([1 2]) % Two groups of bars
    xticklabels({'R1', 'R2'});
    ylim([0,.03]);
    xlabel('Origins');
    ylabel('Inhibition');
    legend('to R1', 'to R2');
    mysavefig(h, filename, testdir, 12, [3.4, 2.6]);

    save(SavedResult, 'Ntwk', 'EEdynamc', 'EIdynamc', 'IEdynamc', ...
        'Rsample', 'Gsample', 'Rsection', 'Gsection', 'Seq', 'dt', 'time', 'smpl_time', 'ti');
end