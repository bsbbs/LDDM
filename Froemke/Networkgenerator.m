Ntwkfile = fullfile(gnrloutdir, 'Ntkw.mat');
if ~exist(Ntwkfile, 'file')
    %% The structure of the network
    Ntwk.Scale = 8000; % scale of the micro structure to test, 2D square (length and width) in unit of um
    Ntwk.Exct.Location = []; % the physical location of excitatory cells
    Ntwk.Exct.N = 400; % number of the excitarory cells
    Ntwk.Exct.Properties.size =15; % physical size of the soma, in um
    Ntwk.Exct.Properties.FWHM = 1500; % the Full width at half maximum(FWHM) range of axons and dendrites connection
    Ntwk.Inhbt.Location = [];
    Ntwk.Inhbt.N = 100;
    Ntwk.Inhbt.Properties.size = 8;
    Ntwk.Inhbt.Properties.FWHM = 500;
    Ntwk.Noise.tauN = .002; % time constant for the OU process of noise
    Ntwk.Noise.sgm = 10; % amplitude of noise
    Ntwk.tauR = .1853; % time constant for the three types of neurons
    Ntwk.tauG = .2244;
    Ntwk.tauD = .3231;
    Ntwk.c = .01; % connection weight from input to R
    Ntwk.w = .01; % connection weight from R to G
    Ntwk.alpha = .01; % self-excitation weight from R to R
    % time constant of synaptic plasticity for pre-post and post-pre kernels
    Ntwk.tau_prepost = .02; % 20 ms
    Ntwk.tau_postpre = .02; % 20 ms
    Ntwk.A = .001; % the maximum amplitude change on the synapctic plasticity for each pair of spikes
    % the time constant of synaptic weight change for E and I neurons
    % according to Feild,... Freomke, 2020, Neuron, Figure 4
    % The time constants for homosynaptic plasticity are 1.0 min and 0.6 min for excitation and inhibition, respectively
    % The time constants for heterosynaptic plasticity are 4.8 and 9.1 min for excitation and inhibition, respectively
    Ntwk.tauhet_e = 4.8*60;
    Ntwk.tauhom_e = 1*60;
    Ntwk.tauhet_i = 9.1*60;
    Ntwk.tauhom_i = .6*60;
    Ntwk.taus = Ntwk.tauhet_i; % we assume that the decaying time constant for the synaptic weight follows the heterosynpatic decaying dynamic
    Ntwk.gamma = 1; % the change rate of synaptic weight

    %% Physical location of the cells, assuming located on the same layer, thus
    % spreading the 2D intersection slice
    %% locations of excitatory cells
    rng(2023);
    xE = -Ntwk.Scale/2 + rand(Ntwk.Exct.N,1)*Ntwk.Scale;
    [xE, I] = sort(xE);
    %yE = zeros(Ntwk.Exct.N,1);
    yE = randn(Ntwk.Exct.N,1)*Ntwk.Scale/24;
    Ntwk.Exct.Location = [xE, yE];
    %% locations of Inhibitory cells
    rng(2023);
    xI = -Ntwk.Scale/2 + rand(Ntwk.Inhbt.N,1)*Ntwk.Scale;
    [xI, Ip] = sort(xI);
    yI = -Ntwk.Scale/4 + rand(Ntwk.Inhbt.N,1)*Ntwk.Scale/2;
    Ntwk.Inhbt.Location = [xI, yI];
    if show
        h = figure;
        filename = 'Network structure';
        hold on;
        plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/3, 'MarkerFaceColor','auto');
        plot(Ntwk.Inhbt.Location(:,1), Ntwk.Inhbt.Location(:,2), '.', 'Color', [.5, .5, .5], 'MarkerSize', Ntwk.Inhbt.Properties.size);
        xlabel('\mum');
        ylabel('\mum');
        xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
        ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
    end
    % connection from E -> I
    [XE, XI] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1)); % rows represent Inhbt and columns represent Exct
    [YE, YI] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Inhbt.Location(:,2));
    DstcEI = sqrt((XE - XI).^2 + (YE - YI).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
    % p_EI = exp(-DstcEI/Ntwk.Exct.Properties.FWHM); % probability of connection from E to I based on distance
    p_EI = 1 - (1 - exp(-DstcEI/Ntwk.Exct.Properties.FWHM)).^4; % probability of connection from E to I based on distance
    % connection from I -> E
    [XI, XE] = meshgrid(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent Exct and columns represent Inhbt
    [YI, YE] = meshgrid(Ntwk.Inhbt.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcIE = sqrt((XI - XE).^2 + (YI - YE).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
    p_IE = 1 - (1 - exp(-DstcIE/Ntwk.Inhbt.Properties.FWHM)).^4; % probability of connection from I to E based on distance
    % p_IE = exp(-DstcIE/Ntwk.Inhbt.Properties.FWHM); % probability of connection from I to E based on distance
    % connection from E -> E
    [XE1, XE2] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1));
    [YE1, YE2] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcEE = sqrt((XE1 - XE2).^2 + (YE1 - YE2).^2); %  Euclidean distance between each pair of Exct neurons
    p_EE = 1 - (1 - exp(-DstcEE/Ntwk.Exct.Properties.FWHM)).^4; % probability of connection from E to E based on distance
    % p_EE = exp(-DstcEE/Ntwk.Exct.Properties.FWHM); % probability of connection from E to E based on distance
    Cnnct_EI = p_EI >= rand(size(p_EI)); % connections from E to I, 0 or 1
    Cnnct_IE = p_IE >= rand(size(p_IE)); % connections from I to E, 0 or 1
    Cnnct_EE = p_EE >= rand(size(p_EE)); % connections from E to E, 0 or 1
    Ntwk.Cnnct_EI = Cnnct_EI;
    Ntwk.Cnnct_IE = Cnnct_IE;
    Ntwk.Cnnct_EE = Cnnct_EE;
    if show
        exampleE = 280;
        EIs = find(Cnnct_EI(:, exampleE));
        for i = EIs'
            lgd1 = plot([Ntwk.Exct.Location(exampleE,1), Ntwk.Inhbt.Location(i,1)],...
                [Ntwk.Exct.Location(exampleE,2), Ntwk.Inhbt.Location(i,2)],...
                '-', 'Color', OKeeffe(3,:), 'LineWidth',1);
        end
        %plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        exampleI = 48;%14;
        IEs = find(Cnnct_IE(:, exampleI));
        for i = IEs'
            lgd2 = plot([Ntwk.Exct.Location(i,1), Ntwk.Inhbt.Location(exampleI,1)],...
                [Ntwk.Exct.Location(i,2), Ntwk.Inhbt.Location(exampleI,2)],...
                '-', 'Color', OKeeffe(4,:), 'LineWidth',1);
        end
        %plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        legend([lgd1, lgd2], {'R to G','G to R'});
        mysavefig(h, filename, gnrloutdir, 14, [5, 2.5], 2);
    end
    %% Synaptic connections
    % initial weights
    wEE_initial = .01*p_EE.*Ntwk.Cnnct_EE; % synaptic weight from E to E, weak initial connections of self-excitation
    wEI_initial = .01*p_EI.*Ntwk.Cnnct_EI; % synaptic weight from E to I, weak initial connections from E to I
    wIE_initial = .01*p_IE.*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST

    Ntwk.wEE_initial = wEE_initial;
    Ntwk.wEI_initial = wEI_initial;
    Ntwk.wIE_initial = wIE_initial;
    if show
        h = figure;
        filename = 'wEE_initial';
        imagesc(wEE_initial');
        caxis([0, .1]);
        colormap(bluewhitered(256));
        c = colorbar;
        c.Label.String = 'Initial synaptic weights';
        c.Location = 'northoutside';
        xlabel('Exct channels');
        ylabel('Exct channels');
        mysavefig(h, filename, gnrloutdir, 14, [3.2, 3.6]); % [4, 3]

        h = figure;
        filename = 'wEI_initial';
        imagesc(wEI_initial');
        caxis([0, .1]);
        colormap(bluewhitered(256));
        xlabel('Inhbt channels');
        ylabel('Exct channels');
        mysavefig(h, filename, gnrloutdir, 14, [1.6, 3]);

        h = figure;
        filename = 'wIE_initial';
        imagesc(wIE_initial');
        caxis([0, .1]);
        colormap(bluewhitered(256));
        ylabel('Inhbt channels');
        xlabel('Exct channels');
        mysavefig(h, filename, gnrloutdir, 14, [3.2, 1.46]);
    end
    save(Ntwkfile, 'Ntwk');
else
    load(Ntwkfile);
end
