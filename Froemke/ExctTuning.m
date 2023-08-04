% Excitatory neurons, tuning to the inputs, binary as an example
%% setup
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionPP';
addpath('./utils/');
addpath('./functions/');
addpath('../utils/bluewhitered/');
pastelPalette = [ ...
    0.8941, 0.4471, 0.2471;  % Pastel Red
    0.3686, 0.6118, 0.8275;  % Pastel Blue
    0.4314, 0.8039, 0.6314;  % Pastel Green
    0.9569, 0.7608, 0.5176;  % Pastel Yellow
    0.8157, 0.5412, 0.7647;  % Pastel Purple
    1.0000, 0.6471, 0.5098;  % Pastel Orange
    ];
OKeeffe = [
    255, 192, 203;  % Pink (Pale Violet Red)
    100, 149, 237;  % Cornflower Blue
    255, 127, 80;   % Coral
    144, 238, 144;  % Light Green (Pale Green)
    255, 228, 196;  % Bisque
    147, 112, 219;  % Medium Purple
    0, 206, 209;    % Dark Turquoise
    250, 128, 114;  % Salmon
    152, 251, 152;  % Pale Green (Light Green)
    218, 112, 214;  % Orchid
    ]/255;
h = figure;
hold on;
for i = 1:size(OKeeffe,1)
    plot(i, 0, '.', 'Color', OKeeffe(i,:), 'MarkerSize',180);
end
%% The structure of the network
Ntwk.Scale = 4000; % scale of the micro structure to test, 2D square (length and width) in unit of um
Ntwk.Input.N = 2; % number of inputs 
Inputx = linspace(-.5, .5, Ntwk.Input.N+2);
Ntwk.Input.Location = Ntwk.Scale*[Inputx(2:end-1)', ones(Ntwk.Input.N, 1)*.25]; % the location of the input intervene
Ntwk.Input.Variance = 200; % the range of each input intervening
Ntwk.Exct.Location = []; % the physical location of excitatory cells
Ntwk.Exct.N = 100; % number of the excitarory cells
Ntwk.Exct.Properties.size = 25; % physical size of the soma, in um
Ntwk.Exct.Properties.FWHM = 1500; % the Full width at half maximum(FWHM) range of axons and dendrites connection
Ntwk.Inhbt.Location = [];
Ntwk.Inhbt.N = 50;
Ntwk.Inhbt.Properties.size = 12;
Ntwk.Inhbt.Properties.FWHM = 500;
Ntwk.Noise.tauN = .002; % time constant for the OU process of noise
Ntwk.Noise.sgm = 3; % amplitude of noise
Ntwk.tauR = .1853; % time constant for the three types of neurons
Ntwk.tauG = .2244;
Ntwk.tauD = .3231;
Ntwk.w = 1; % connection weight from R to G
% the time constant of synaptic plasticity for E and I neurons
% according to Feild,... Freomke, 2020, Neuron, Figure 4
% The time constants for homosynaptic plasticity are 1.0 min and 0.6 min for excitation and inhibition, respectively
% The time constants for heterosynaptic plasticity are 4.8 and 9.1 min for excitation and inhibition, respectively
Ntwk.tauhet_e = 4.8*60;
Ntwk.tauhom_e = 1*60; 
Ntwk.tauhet_i = 9.1*60;
Ntwk.tauhom_i = .6*60; 
Ntwk.gamma = .001; % the changing rate of potentiation
Ntwk.show = 1; % control for visualization
%% physical location of the cells, assuming located on the same layer, thus
% spreading on the 2D surface
rng(2023);
xE = -Ntwk.Scale/2 + rand(Ntwk.Exct.N,1)*Ntwk.Scale;
yE = zeros(Ntwk.Exct.N,1);
Ntwk.Exct.Location = [xE, yE];

% tuning curve of the excitatory cells to the inputs 1 and 2, with each
% Exct cell receiving random connections from each input based on the
% distance from the center of each input
Coupling = [];
for i = 1:Ntwk.Input.N
    Coupling(:,i) = normpdf(Ntwk.Exct.Location(:,1), Ntwk.Input.Location(i,1), Ntwk.Input.Variance); % connection strength
    % as the probability density in a Gaussian distribution with the mean as
    % the center of the input and variance
end
Ntwk.Exct.tuning = Coupling/max(Coupling(:));
show = 1;
if show
    h = figure;
    filename = 'ExctTuning';
    subplot(2,1,1); hold on;
    % connection strength distribution from each input
    pseudolocation = linspace(-Ntwk.Scale/2, Ntwk.Scale/2, 100);
    legendLabels = cell(1, Ntwk.Input.N);
    for i = 1:Ntwk.Input.N
        plot(pseudolocation, normpdf(pseudolocation, Ntwk.Input.Location(i,1), Ntwk.Input.Variance)/normpdf(Ntwk.Input.Location(i,1), Ntwk.Input.Location(i,1), Ntwk.Input.Variance), 'LineWidth', 2, 'Color', OKeeffe(i,:));
        legendLabels{i} = sprintf('Input %d', i);
    end
    legend(legendLabels, 'Location', 'best');
    xlabel('\mum');
    ylabel('Connection strength');
    subplot(2,1,2); hold on;
    plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/3);
    for i = 1:Ntwk.Input.N
        for j = 1:Ntwk.Exct.N
            plot([Ntwk.Input.Location(i,1), Ntwk.Exct.Location(j,1)], ...
                [Ntwk.Input.Location(i,2), Ntwk.Exct.Location(j,2)],...
                '-', 'Color', OKeeffe(i,:), 'LineWidth', Ntwk.Exct.tuning(j,i)/max(Ntwk.Exct.tuning(:))*3);
        end
    end
    xlabel('\mum');
    ylabel('\mum');
    xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
    ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
    mysavefig(h, filename, outdir, 14, [6, 6]);
end
  
%% locations of Inhibitory cells
rng(2023);
xI = -Ntwk.Scale/2 + rand(Ntwk.Inhbt.N,1)*Ntwk.Scale;
yI = -Ntwk.Scale/4 + rand(Ntwk.Inhbt.N,1)*Ntwk.Scale/2;
Ntwk.Inhbt.Location = [xI, yI];
h = figure;
filename = 'Network structure';
hold on;
plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/3, 'MarkerFaceColor','auto');
plot(Ntwk.Inhbt.Location(:,1), Ntwk.Inhbt.Location(:,2), '.', 'Color', [.5, .5, .5], 'MarkerSize', Ntwk.Inhbt.Properties.size);
xlabel('\mum');
ylabel('\mum');
xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
% connection from E -> I and I -> E
[XE, XI] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1)); % rows represent Inhbt and columns represent Exct
[YE, YI] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Inhbt.Location(:,2));
Dstc = sqrt((XE - XI).^2 + (YE - YI).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
p_EI = exp(-Dstc/Ntwk.Exct.Properties.FWHM); % probability of connection from E to I based on distance
p_IE = exp(-Dstc/Ntwk.Inhbt.Properties.FWHM); % probability of connection from I to E based on distance
Cnnct_EI = p_EI >= rand(size(p_EI)); % connections from E to I, 0 or 1
Cnnct_IE = p_IE >= rand(size(p_IE)); % connections from I to E, 0 or 1
Ntwk.Cnnct_EI = Cnnct_EI;
Ntwk.Cnnct_IE = Cnnct_IE;
example = randi(Ntwk.Exct.N);
EIs = find(Cnnct_EI(:, example));
for i = EIs'
    lgd1 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
        [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
        '-', 'Color', OKeeffe(3,:), 'LineWidth',1.6);
end
plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);

IEs = find(Cnnct_IE(:, example));
for i = IEs'
    lgd2 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
        [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
        '-', 'Color', OKeeffe(4,:), 'LineWidth',1.6);
end
plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
legend([lgd1, lgd2], {'E->I','I->E'});
mysavefig(h, filename, outdir, 14, [6, 3]);

%% Input values and excitatory activities with Background white noise (OU process) input to the excitatory cells
dt = .001; % time precision for simulation, in unit of second
Ntrial = 258;
rng(2023);
Seq = CreateEvents(Ntrial, dt);
% time = [0:(size(Seq, 1)-1)]*dt/60; % unit in Mins
time = [0:20000]*dt; % unit in Secs
value = [20, 20];
R = zeros(Ntwk.Exct.N,numel(time));
G = zeros(Ntwk.Inhbt.N,numel(time));
wIE = .02*randn(Ntwk.Exct.N, Ntwk.Inhbt.N).*Ntwk.Cnnct_IE'; % synaptic weight from I to E, weak initial connections from the nearby SST 
wIE(wIE<0) = 0;
wIE_initial = wIE;
%%
if Ntwk.show
    % synaptic connection weights
    h = figure;
    filename = sprintf('SynapticW_at_time%1.2fs', 0*dt);
    imagesc(wIE);
    caxis([-1, 1]);
    colormap(bluewhitered(256)), colorbar;
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 14, [9, 6]);
end
noise = zeros(Ntwk.Exct.N, 1);
for ti = 1:(length(time)-1)
    if (mod(ti*dt, 1)==0)
        fprintf('%is.',ti*dt*1);
    end
    dR = (-R(:,ti) + Ntwk.Exct.tuning*(Seq(ti,:).*value)'./(1 + wIE*G(:,ti)))*dt/Ntwk.tauR;
    dG = (-G(:,ti) + Ntwk.w*Ntwk.Cnnct_EI*R(:,ti))*dt/Ntwk.tauG;
    preAP = G;
    postAP = R;
    threshR = 5;
    threshG = 50;
    tau_pos = .02; % 20 ms
    tau_neg = .02; % 20 ms
    % Calculate STDP kernel
    synapticIE_rollover = zeros([size(Cnnct_IE), ti]);
    for t = 1:ti % accumulation over time points
        % Potentiation effect
        prebinary = preAP(:,t) >= threshG;
        t_pos = t:min([t+3*tau_pos/dt,ti]);
        postbinary = (postAP(:,t_pos) > threshR)';
        synapticIE_rollover(:,:,t) = prebinary * (exp(-(t_pos - t)*dt / tau_pos) * postbinary);
        % Depression effect
        postbinary = (postAP(:,t) >= threshR)';
        t_neg = max([1, t-3*tau_neg/dt]):t;
        prebinary = (preAP(:,t_neg) > threshG);
        synapticIE_rollover(:,:, t) = synapticIE_rollover(:,:,t) + (prebinary * exp((t_neg - t)*dt / tau_neg)') * postbinary;
    end
    synapticIE_change = sum(synapticIE_rollover,3).*Ntwk.Cnnct_IE;
    dP = (-wIE/Ntwk.tauhet_i + Ntwk.gamma*(1-wIE).*synapticIE_change'/Ntwk.tauhom_i)*dt;
    wIE = wIE + dP;
    noise = noise + (-noise + randn(size(noise))*sqrt(dt)*Ntwk.Noise.sgm)*dt/Ntwk.Noise.tauN;
    R(:,ti+1) = R(:,ti) + dR + noise;
    G(:,ti+1) = G(:,ti) + dG;
    R(R(:, ti+1) < 0, ti+1) = 0;
    G(G(:, ti+1) < 0, ti+1) = 0;
    wIE(wIE<0) = 0;
end
save(fullfile(outdir, 'Simulation20s.mat'), 'Ntwk', 'R', 'G', 'wIE', 'wIE_initial', 'synapticIE_change', 'synapticIE_rollover', 'Seq');
%% visualization
if Ntwk.show
    % kernal
    h = figure; hold on;
    filename = 'Kernal';
    t_pos = -3*tau_pos/dt:0;
    y1 = exp(t_pos*dt / tau_pos);
    t_neg = 0:3*tau_neg/dt;
    y2 = exp(-t_neg*dt / tau_neg);
    plot([t_pos, t_neg], [y1, y2], '-', 'Color',OKeeffe(10,:), 'LineWidth',2);
    plot([0,0],[-1,1], 'k--');
    plot([-60,60],[0,0], 'k--');
    ylim([-1, 1]);
    ylabel('Connection change');
    xlabel('\Delta t_{pre} - t_{post}');
    text(-45,0.6,'Pre-Post');
    text(30,0.6,'Post-Pre');
    mysavefig(h, filename, outdir, 12, [5, 4]);
    % V
    h = figure;
    filename = 'InputDynamic';
    subplot(6,1,1);
    hold on;
    plot(time, Seq(1:numel(time),1)+1.5, '-', "Color", OKeeffe(1,:), 'LineWidth',2);
    plot(time, Seq(1:numel(time),2), '-', "Color", OKeeffe(2,:), 'LineWidth',2);
    ylim([-1, 3]);
    yticks([0,1.5]);
    yticklabels(num2str([1, 2]'));
    ylabel('Events');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [5, 11]);
    % R
    h = figure;
    filename = 'R_Dynamic';
    hold on;
    Amp = max(R(:));
    for i = 1:Ntwk.Exct.N
        plot(time, R(i,:)+i*Amp, 'k-');
    end
    ylim([0, (i+1)*Amp]);
    yticks([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]'));
    ylabel('Channels (Exct)');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [5, 11]);
    % G
    h = figure;
    filename = 'G_Dynamic';
    hold on;
    Amp = max(G(:));
    for i = 1:Ntwk.Inhbt.N
        plot(time, G(i,:)+i*Amp, 'k-');
    end
    ylim([0, (i+1)*Amp]);
    yticks([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]'));
    ylabel('Channels (Inhbt)');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [5, 11]);
    % synapctic change
    h = figure;
    filename = sprintf('PairingSum_at_time%1.2fs', ti*dt);
    imagesc(synapticIE_change);
    colormap(bluewhitered(256)), colorbar;
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 14, [9, 6]);
    % synaptic connection weights
    h = figure;
    filename = sprintf('SynapticW_at_time%1.2fs', ti*dt);
    imagesc(wIE);
    caxis([-1, 1]);
    colormap(bluewhitered(256)), colorbar;
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 14, [9, 6]);
    % synaptic weight change
    h = figure;
    filename = sprintf('SynapticW_change');
    imagesc(wIE - wIE_initial);
    caxis([-1, 1]);
    colormap(bluewhitered(256)), colorbar;
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 14, [9, 6]);
    % Example dynamics R and G
    h = figure;
    minValue = max(synapticIE_change(:));
    [IdxI, IdxE] = find(synapticIE_change == minValue);
    filename = sprintf('Exmplpair_at_E%iandI%i',IdxE(1), IdxI(1));
    hold on;
    plot(time, R(IdxE(1),:)/max(R(IdxE(1),:)), '-', "Color", OKeeffe(3,:), 'LineWidth',1);
    plot(time, G(IdxI(1),:)/max(G(IdxI(1),:)), '-', "Color", OKeeffe(4,:), 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Normalized activity');
    legend({"PYR", "SST"}, 'Location', "best");
    mysavefig(h, filename, outdir, 12, [5, 3]);
    %%
    rng(2023);
    h = figure;
    filename = 'Network structure';
    hold on;
    plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/3, 'MarkerFaceColor','auto');
    plot(Ntwk.Inhbt.Location(:,1), Ntwk.Inhbt.Location(:,2), '.', 'Color', [.5, .5, .5], 'MarkerSize', Ntwk.Inhbt.Properties.size);
    xlabel('\mum');
    ylabel('\mum');
    xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
    ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
    % trained connection from E -> I and I -> E
    EIs = find(Cnnct_EI(:, example));
    for i = EIs'
        lgd1 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
            [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
            '-', 'Color', OKeeffe(3,:), 'LineWidth',1.6);
    end
    plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);

    IEs = find(Cnnct_IE(:, example));
    for i = IEs'
        lgd2 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
            [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
            '-', 'Color', OKeeffe(4,:), 'LineWidth',1.6);
    end
    plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
    legend([lgd1, lgd2], {'E->I','I->E'});
    mysavefig(h, filename, outdir, 14, [6, 3]);
end