%% Input features
Ntwk.Input.N = Ninput; % number of inputs
Inputx = linspace(-.5, .5, 1+Ntwk.Input.N*2);
Inputx = Inputx(2:2:(numel(Inputx)-1));
Ntwk.Input.Location = Ntwk.Scale*[Inputx', ones(Ntwk.Input.N, 1)*.25]; % the location of the input intervene
Ntwk.Input.Variance = 100; % the range of each input intervening
%% Tuning curve of the excitatory cells to the inputs, with each
% Exct cell receiving random connections from each input based on the
% distance from the center of each input
rng(2023);
Coupling = [];
for i = 1:Ntwk.Input.N
    Coupling(:,i) = normpdf(Ntwk.Exct.Location(:,1), Ntwk.Input.Location(i,1), Ntwk.Input.Variance); % connection strength
    % as the probability density in a Gaussian distribution with the mean as
    % the center of the input and variance
end
Ntwk.Exct.tuning = Coupling/max(Coupling(:));
Ntwk.Visualization.IdxE = find(Ntwk.Exct.tuning(:,1) == max(Ntwk.Exct.tuning(:,1)));
Ntwk.Visualization.IdxI = find(Ntwk.wIE(Ntwk.Visualization.IdxE,:) == max(Ntwk.wIE(Ntwk.Visualization.IdxE,:)));
save(fullfile(Inputdir, 'NtwkwithInput.mat'), "Ntwk");
if show
    h = figure;
    filename = 'InputTuning';
    subplot(2,1,1); hold on;
    % connection strength distribution from each input
    pseudolocation = linspace(-Ntwk.Scale/2, Ntwk.Scale/2, 400);
    legendLabels = cell(1, Ntwk.Input.N);
    for i = 1:Ntwk.Input.N
        y = normpdf(pseudolocation, Ntwk.Input.Location(i,1), Ntwk.Input.Variance)/normpdf(Ntwk.Input.Location(i,1), Ntwk.Input.Location(i,1), Ntwk.Input.Variance);
        plot(pseudolocation, y, 'LineWidth', 2, 'Color', OKeeffe(i,:));
        legendLabels{i} = sprintf('Input %d', i);
    end
    legend(legendLabels, 'Location', 'best');
    xlabel('\mum');
    ylabel('Coupling strength');
    mysavefig(h, filename, outdir, 12, [3, 3], 2);
    subplot(2,1,2); hold on;
    plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/5);
    for i = 1:Ntwk.Input.N
        for j = 1:Ntwk.Exct.N
            weight = Ntwk.Exct.tuning(j,i)/max(Ntwk.Exct.tuning(:));
            if weight > 1e-6
                plot([Ntwk.Input.Location(i,1), Ntwk.Exct.Location(j,1)], ...
                    [Ntwk.Input.Location(i,2), Ntwk.Exct.Location(j,2)],...
                    '-', 'Color', OKeeffe(i,:), 'LineWidth', weight*1);
            end
        end
    end
    xlabel('\mum');
    ylabel('\mum');
    xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
    ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
    mysavefig(h, filename, Inputdir, 14, [3, 3], 2);
end
