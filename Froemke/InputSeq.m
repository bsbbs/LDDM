%% V input(s) sequence
if show
    h = figure;
    filename = 'InputDynamic';
    subplot(6,1,1);
    hold on;
    for ii = 1:Ntwk.Input.N
        plot(inpt_time, Seq(1:numel(inpt_time),ii)+1.5*(ii-1), '-', "Color", OKeeffe(ii,:), 'LineWidth',1);
    end
    ylim([-.5, Ntwk.Input.N*1.5]);
    yticks(1.5*([1:Ntwk.Input.N]-1));
    yticklabels(num2str([1:Ntwk.Input.N]'));
    ylabel('Input');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 14, [3, 11], 2);
    subplot(6,1,2);
    hold on;
    for ii = 1:Ntwk.Input.N
        plot(time, Seq(1:numel(time),ii)+1.5*(ii-1), '-', "Color", OKeeffe(ii,:), 'LineWidth',.5);
    end
    ylim([-.5, Ntwk.Input.N*1.5]);
    plot([min(inpt_time), min(inpt_time), max(inpt_time), max(inpt_time), min(inpt_time)]/60, [-.5, Ntwk.Input.N*1.5, Ntwk.Input.N*1.5, -.5, -.5], 'k-');
    yticks(1.5*([1:Ntwk.Input.N]-1));
    yticklabels(num2str([1:Ntwk.Input.N]'));
    xlim([0, 35]);
    ylabel('Input');
    xlabel('Time (mins)');
    mysavefig(h, filename, Inputdir, 14, [3, 11], 2);
end