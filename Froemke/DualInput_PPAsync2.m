% Excitatory neurons, tuning to the inputs, binary as an example
%% setup
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionDual';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
Setup;
%% The structure of the network
Ntwk.Scale = 8000; % scale of the micro structure to test, 2D square (length and width) in unit of um
Ntwk.Input.N = 2; % number of inputs 
Inputx = linspace(-.5, .5, 1+Ntwk.Input.N*2);
Inputx = Inputx(2:2:(numel(Inputx)-1));
Ntwk.Input.Location = Ntwk.Scale*[Inputx', ones(Ntwk.Input.N, 1)*.25]; % the location of the input intervene
Ntwk.Input.Variance = 100; % the range of each input intervening
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
Ntwk.show = 1; % control for visualization
%% physical location of the cells, assuming located on the same layer, thus
% spreading on the 2D surface
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
h = figure;
filename = 'Network structure';
hold on;
plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', Ntwk.Exct.Properties.size/3, 'MarkerFaceColor','auto');
plot(Ntwk.Inhbt.Location(:,1), Ntwk.Inhbt.Location(:,2), '.', 'Color', [.5, .5, .5], 'MarkerSize', Ntwk.Inhbt.Properties.size);
xlabel('\mum');
ylabel('\mum');
xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
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
example = round(Ntwk.Exct.N/2);
EIs = find(Cnnct_EI(:, example));
for i = EIs'
    lgd1 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
        [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
        '-', 'Color', OKeeffe(3,:), 'LineWidth',1);
end
plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);

IEs = find(Cnnct_IE(example, :));
for i = IEs
    lgd2 = plot([Ntwk.Exct.Location(example,1), Ntwk.Inhbt.Location(i,1)],...
        [Ntwk.Exct.Location(example,2), Ntwk.Inhbt.Location(i,2)],...
        '-', 'Color', OKeeffe(4,:), 'LineWidth',1);
end
plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
legend([lgd1, lgd2], {'E to I','I to E'});
mysavefig(h, filename, outdir, 12, [5, 2.5], 2);
%% Synaptic connections
% initial weights
wEE_initial = .01*p_EE.*Ntwk.Cnnct_EE; % synaptic weight from E to E, weak initial connections of self-excitation
wEI_initial = .01*p_EI.*Ntwk.Cnnct_EI; % synaptic weight from E to I, weak initial connections from E to I
wIE_initial = .01*p_IE.*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST 

h = figure;
filename = 'wEE_initial';
imagesc(wEE_initial');
caxis([0, .03]);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Initial synaptic weights';
c.Location = 'northoutside';
xlabel('Exct channels');
ylabel('Exct channels');
mysavefig(h, filename, outdir, 12, [3.2, 3.6]); % [4, 3]

h = figure;
filename = 'wEI_initial';
imagesc(wEI_initial');
caxis([0, .03]);
colormap(bluewhitered(256));
xlabel('Inhbt channels');
ylabel('Exct channels');
mysavefig(h, filename, outdir, 12, [1.6, 3]);

h = figure;
filename = 'wIE_initial';
imagesc(wIE_initial');
caxis([0, .03]);
colormap(bluewhitered(256));
ylabel('Inhbt channels');
xlabel('Exct channels');
mysavefig(h, filename, outdir, 12, [3.2, 1.46]);

%% tuning curve of the excitatory cells to the inputs 1 and 2, with each
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
show = 1;
if show
    h = figure;
    filename = 'InputTuning';
    subplot(2,1,1); hold on;
    % connection strength distribution from each input
    pseudolocation = linspace(-Ntwk.Scale/2, Ntwk.Scale/2, 400);
    legendLabels = cell(1, Ntwk.Input.N);
    for i = 1:Ntwk.Input.N
        plot(pseudolocation, normpdf(pseudolocation, Ntwk.Input.Location(i,1), Ntwk.Input.Variance)/normpdf(Ntwk.Input.Location(i,1), Ntwk.Input.Location(i,1), Ntwk.Input.Variance), 'LineWidth', 2, 'Color', OKeeffe(i,:));
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
            if weight > 0
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
    mysavefig(h, filename, outdir, 12, [3, 3], 2);
end
%% Input values and seqeuences
outdir = fullfile(outdir, 'PPAsync');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
dt = .001; % time precision for simulation, in unit of second
Ntrial = 258;
rng(2023);
Seq = CreateEvents(Ntrial, dt);
value = 20*ones(1, Ntwk.Input.N);
time = [0:(size(Seq, 1)-1)]*dt/60; % unit in Mins
% time = [0:1000*60]*dt/60;
smpl_time = 5*60+[10000:30000]*dt; % unit in Secs
inpt_time = [0:40000]*dt;
% V input(s)
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
mysavefig(h, filename, outdir, 12, [3, 11], 2);
subplot(6,1,2);
hold on;
for ii = 1:Ntwk.Input.N
    plot(time, Seq(1:numel(time),ii)+1.5*(ii-1), '-', "Color", OKeeffe(ii,:), 'LineWidth',.5);
end
ylim([-.5, Ntwk.Input.N*1.5]);
plot([min(inpt_time), min(inpt_time), max(inpt_time), max(inpt_time), min(inpt_time)]/60, [-.5, Ntwk.Input.N*1.5, Ntwk.Input.N*1.5, -.5, -.5], 'r-');
yticks(1.5*([1:Ntwk.Input.N]-1));
yticklabels(num2str([1:Ntwk.Input.N]'));
xlim([0, 35]);
ylabel('Input');
xlabel('Time (mins)');
mysavefig(h, filename, outdir, 12, [3, 11], 2);

%% simulation begin ... excitatory activities with Background white noise (OU process) input to the excitatory cells
R = zeros(Ntwk.Exct.N,numel(time));
G = zeros(Ntwk.Inhbt.N,numel(time));
noise = zeros(Ntwk.Exct.N, 1);
Ntwk.wEE = wEE_initial;
Ntwk.wEI = wEI_initial;
Ntwk.wIE = wIE_initial;
IdxE = find(Ntwk.Exct.tuning(:,1) == max(Ntwk.Exct.tuning(:,1)));
IdxI = find(Ntwk.wIE(IdxE,:) == max(Ntwk.wIE(IdxE,:)));
EEdynamc = Ntwk.wEE(IdxE,IdxE);
EIdynamc = Ntwk.wEI(IdxI,IdxE);
IEdynamc = Ntwk.wIE(IdxE,IdxI);
a_pre_EE = zeros(1, Ntwk.Exct.N);
a_post_EE = zeros(Ntwk.Exct.N, 1);
a_pre_EI = zeros(1, Ntwk.Exct.N);
a_post_EI = zeros(Ntwk.Inhbt.N, 1);
a_pre_IE = zeros(1, Ntwk.Inhbt.N);
a_post_IE = zeros(Ntwk.Exct.N, 1);
for ti = 1:(length(time)-1)
    if (mod(ti*dt, 1)==0)
        fprintf('%is.',ti*dt*1);
    end
    if (mod(ti*dt, 60)==0)
        fprintf('\n');
    end
    % update R and G activities
    dR = (-R(:,ti) + (Ntwk.wEE*R(:,ti) + Ntwk.Exct.tuning*(Seq(ti,:).*value)' + noise)./(1 + Ntwk.wIE*G(:,ti)))*dt/Ntwk.tauR;
    dG = (-G(:,ti) + Ntwk.wEI*R(:,ti))*dt/Ntwk.tauG;
    R(:,ti+1) = R(:,ti) + dR;
    G(:,ti+1) = G(:,ti) + dG;
    R(R(:, ti+1) < 0, ti+1) = 0;
    G(G(:, ti+1) < 0, ti+1) = 0;
    % update OU noise
    noise = noise + (-noise + randn(size(noise))*sqrt(dt)*Ntwk.Noise.sgm)*dt/Ntwk.Noise.tauN;
    % apply E to E STDP kernel
    preAP = R(:,ti)';
    postAP = R(:,ti);
    a_pre_EE = a_pre_EE + (-a_pre_EE/Ntwk.tau_prepost)*dt + Ntwk.A*100*preAP*dt; % the intermiediate variables for pre
    a_post_EE = a_post_EE + (-a_post_EE/Ntwk.tau_postpre)*dt + Ntwk.A*100*postAP*dt; % and post action potential induced plasticity
    a = postAP*dt*a_pre_EE - a_post_EE*preAP*dt; % the AP induced overall plasticity change
    ds = Ntwk.gamma*(1 - Ntwk.wEE).*a.*Ntwk.Cnnct_EE; % a already contains the dt information, only applied to physical connections
    % (-Ntwk.wEE/Ntwk.taus)*dt + 
    Ntwk.wEE = Ntwk.wEE + ds;
    Ntwk.wEE(Ntwk.wEE<0) = 0;
    EEdynamc(ti+1) = Ntwk.wEE(IdxE,IdxE);
    % apply E to I STDP kernel
    preAP = R(:,ti)';
    postAP = G(:,ti);
    a_pre_EI = a_pre_EI + (-a_pre_EI/Ntwk.tau_prepost)*dt + Ntwk.A*100*preAP*dt; % the intermiediate variables for pre
    a_post_EI = a_post_EI + (-a_post_EI/Ntwk.tau_postpre)*dt + Ntwk.A*100*postAP*dt; % and post action potential induced plasticity
    a = postAP*dt*a_pre_EI - a_post_EI*preAP*dt; % the AP induced overall plasticity change
    ds = Ntwk.gamma*(1 - Ntwk.wEI).*a.*Ntwk.Cnnct_EI; % a already contains the dt information, only applied to physical connections
    % (-Ntwk.wEI/Ntwk.tauhet_e)*dt +
    Ntwk.wEI = Ntwk.wEI + ds;
    Ntwk.wEI(Ntwk.wEI<0) = 0;
    EIdynamc(ti+1) = Ntwk.wEI(IdxI,IdxE);
    % apply I to E STDP kernel to the course  of firing rates
    preAP = G(:,ti)';
    postAP = R(:,ti);
    a_pre_IE = a_pre_IE + (-a_pre_IE/Ntwk.tau_prepost)*dt + Ntwk.A*preAP*dt; % the intermiediate variables for pre
    a_post_IE = a_post_IE + (-a_post_IE/Ntwk.tau_postpre)*dt + Ntwk.A*postAP*dt; % and post action potential induced plasticity
    a = postAP*dt*a_pre_IE + a_post_IE*preAP*dt; % the AP induced overall plasticity change
    ds = Ntwk.gamma*(1 - Ntwk.wIE).*a.*Ntwk.Cnnct_IE; % a already contains the dt information, only applied to physical connections
    % (-Ntwk.wIE/Ntwk.taus)*dt + 
    Ntwk.wIE = Ntwk.wIE + ds;
    Ntwk.wIE(Ntwk.wIE<0) = 0;
    IEdynamc(ti+1) = Ntwk.wIE(IdxE,IdxI);
end
fprintf('\n');
save(fullfile(outdir, 'Simulation_full.mat'), 'Ntwk', 'wEE_initial', 'wEI_initial', 'wIE_initial', 'EEdynamc', 'EIdynamc', 'IEdynamc', 'Seq', 'ti', 'time', 'IdxE', 'IdxI');

%% visualization
if Ntwk.show
    % load(fullfile(outdir, 'Simulation_full.mat'));

    h = figure;
    filename = sprintf('wEE_%1.0fmins', max(time));
    imagesc(Ntwk.wEE');
    caxis([-.1, .1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'Trained synptic weights';
    c.Location = 'northoutside';
    xlabel('Exct channels');
    ylabel('Exct channels');
    mysavefig(h, filename, outdir, 12, [3.2, 3.6]); % [4, 3]

    h = figure;
    filename = sprintf('wEI_%1.0fmins', max(time));
    imagesc(Ntwk.wEI');
    caxis([-.35, .35]);
    colormap(bluewhitered(256));
    xlabel('Inhbt channels');
    ylabel('Exct channels');
    mysavefig(h, filename, outdir, 12, [1.6, 3]);

    h = figure;
    filename = sprintf('wIE_%1.0fmins', max(time));
    imagesc(Ntwk.wIE');
    caxis([-.2, .2]);
    colormap(bluewhitered(256));
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 12, [3.2, 1.46]);


    % synaptic weight change
    h = figure;
    filename = 'wEE_change';
    imagesc((Ntwk.wEE - wEE_initial)');
    caxis([-.1, .1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'Synaptic weight change';
    c.Location = 'northoutside';
    xlabel('Exct channels');
    ylabel('Exct channels');
    mysavefig(h, filename, outdir, 12, [3.2, 3.6]); % [4, 3]

    h = figure;
    filename = 'wEI_change';
    imagesc((Ntwk.wEI - wEI_initial)');
    caxis([-.35, .35]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'Weights change';
    xlabel('Inhbt channels');
    ylabel('Exct channels');
    mysavefig(h, filename, outdir, 12, [1.6, 3]);

    % synaptic weight change
    h = figure;
    filename = 'wIE_change';
    imagesc((Ntwk.wIE - wIE_initial)');
    caxis([-.2, .2]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'Weights change';
    xlabel('Exct channels');
    ylabel('Inhbt channels');
    mysavefig(h, filename, outdir, 12, [3.2, 1.46]);

    % dynamics
    smpl_time = 5*60+[10000:30000]*dt; % unit in Secs
    segment = round(smpl_time/dt);
    % R
    h = figure;
    filename = sprintf('R_Dynamic%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    hold on;
    interval = 4;
    Amp = max(max(R(:,segment)))/interval;
    for i = 1:interval:Ntwk.Exct.N
        plot(smpl_time, R(i,segment)+i*Amp, 'k-');
    end
    ylim([0, (i+interval)*Amp]);
    yticks([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]'));
    ylabel('Exct channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [3, 11.2]);
    % G
    h = figure;
    filename = sprintf('G_Dynamic%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    hold on;
    Amp = max(max(G(:,segment)))/interval;
    for i = 1:interval:Ntwk.Inhbt.N
        plot(smpl_time, G(i,segment)+i*Amp, 'k-');
    end
    ylim([0, (i+interval)*Amp]);
    yticks([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]'));
    ylabel('Inhbt channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [3, 3.5]);
    
    % Example dynamics R and G
    h = figure;
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fs_%1.0fs',IdxE(1), IdxI(1), min(smpl_time), max(smpl_time));
    hold on;
    % plot(smpl_time, R(IdxE(1),1:numel(smpl_time)), '-', "Color", OKeeffe(3,:), 'LineWidth',1);
    % plot(smpl_time, G(IdxI(1),1:numel(smpl_time)), '-', "Color", OKeeffe(4,:), 'LineWidth',1);
    segment = round(smpl_time/dt);
    plot(smpl_time, G(IdxI(1),segment)/max(G(IdxI(1),segment)), '-', "Color", OKeeffe(4,:), 'LineWidth',.5);
    plot(smpl_time, R(IdxE(1),segment)/max(R(IdxE(1),segment)), '-', "Color", OKeeffe(3,:), 'LineWidth',.5);
    xlabel('Time (secs)');
    ylabel('Normalized activity');
    %ylabel('Firing rates (Hz)');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);

    h = figure;
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fmins',IdxE(1), IdxI(1), max(time));
    hold on;
    lg2 = plot(time, G(IdxI(1),1:numel(time)), '-', "Color", OKeeffe(4,:), 'LineWidth', .5);
    lg1 = plot(time, R(IdxE(1),1:numel(time)), '-', "Color", OKeeffe(3,:), 'LineWidth', .5);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [6.1, 6.1, .1, .1, 6.1], '-k');
    xlim([0, 35]);
    ylabel('Firing rates (Hz)');
    xlabel('Time (mins)');
    legend([lg1, lg2], {"R", "G"}, 'Location', "best", 'FontSize', 10);
    mysavefig(h, filename, outdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wEE_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, EEdynamc(segment), 'k-', 'LineWidth',1);
    ylim([0,.02]);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);

    h = figure; hold on;
    filename = 'wEE_dynamic';
    plot(time, EEdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [.011, .011, .009, .009, 0.011], '-k');
    xlim([0, 35]);
    ylim([0,.02]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wEI_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, EIdynamc(segment), 'k-', 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);

    h = figure; hold on;
    filename = 'wEI_dynamic';
    plot(time, EIdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [.058, 0.058, .038, .038, .058], '-k');
    xlim([0, 35]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = sprintf('wIE_dynamic_%1.0fs_%1.0fs', min(smpl_time), max(smpl_time));
    plot(smpl_time, IEdynamc(segment), 'k-', 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);
    
    h = figure; hold on;
    filename = 'wIE_dynamic';
    plot(time, IEdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [.062, .062, .053, .053, .062]+.006, '-k');
    xlim([0, 35]);
    xlabel('Time (mins)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [3, 2], 2);


    %%
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
    mysavefig(h, filename, outdir, 12, [3, 1]);
    
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
    mysavefig(h, filename, outdir, 12, [3, 1]);
    
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
    mysavefig(h, filename, outdir, 12, [3, 1]);
    
    % h = figure;
    % filename = 'Tuning_Inhbt_IE_sort';
    % [IE1, I] = sort(Ntwk.Inhbt.tuningIE(:,1), 'descend');
    % IE2 = Ntwk.Inhbt.tuningIE(I,2);
    % imagesc([IE1, IE2]');
    % caxis([0, max(Ntwk.Inhbt.tuningIE(:))]);
    % colormap(bluewhitered(256));
    % c = colorbar;
    % xlabel('Inhbt channels');
    % yticks([1,2]);
    % yticklabels({'Input 1', 'Input 2'});
    % c.Label.String = 'Selective Inhibition';
    % mysavefig(h, filename, outdir, 12, [3, 1]);

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
            text(i-.3,j,sprintf('%.2e',Ntwk.InhbtMtrx(i,j)));
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
    %xlabel('E -> I');
    %ylabel('I -> E');
    mysavefig(h, filename, outdir, 12, [3.4, 2.6]);
    save(fullfile(outdir, 'Simulation_full.mat'), 'Ntwk', 'wEE_initial', 'wEI_initial', 'wIE_initial', 'EEdynamc', 'EIdynamc', 'IEdynamc', 'Seq', 'ti', 'IdxE', 'IdxI');

end


%% Coding in well-trained network
time = [0:10000]*dt; % unit in Secs
Seq = ones(numel(time),1);
values = 0:100;
avgR = zeros(Ntwk.Exct.N, numel(values));
for vi = 1:numel(values)
    fprintf('Value = %i: ',values(vi));
    R = zeros(Ntwk.Exct.N,numel(time));
    G = zeros(Ntwk.Inhbt.N,numel(time));
    noise = zeros(Ntwk.Exct.N, 1);
    for ti = 1:(length(time)-1)
        if (mod(ti*dt, 1)==0)
            fprintf('%is.',ti*dt*1);
        end
        % update R and G activities
        dR = (-R(:,ti) + (Ntwk.alpha*Ntwk.Cnnct_EE*R(:,ti) + Ntwk.Exct.tuning*(Seq(ti,:).*values(vi)) + noise)./(1 + s*G(:,ti)))*dt/Ntwk.tauR;
        dG = (-G(:,ti) + Ntwk.w*Ntwk.Cnnct_EI*R(:,ti))*dt/Ntwk.tauG;
        R(:,ti+1) = R(:,ti) + dR;
        G(:,ti+1) = G(:,ti) + dG;
        R(R(:, ti+1) < 0, ti+1) = 0;
        G(G(:, ti+1) < 0, ti+1) = 0;
        % update OU noise
        noise = noise + (-noise + randn(size(noise))*sqrt(dt)*Ntwk.Noise.sgm)*dt/Ntwk.Noise.tauN;
    end
    fprintf('\n');
    avgR(:, vi) = R(:, end);
end
%% 
h = figure;
filename = 'ResponseCurves';
IdxE = find(avgR(:,values == 20) == max(avgR(:,values == 20)));
subplot(2,1,1); hold on;
plot(values, avgR(IdxE,:), '-', 'Color',OKeeffe(3,:), 'LineWidth',2);
xlabel('Input value');
ylabel('Firing rates (Hz)');
mysavefig(h, filename, outdir, 12, [3, 4], 2);

subplot(2,1,2); hold on;
plot((values(1:end-1) + values(2:end))/2, diff(avgR(IdxE,:)), '-', 'Color',OKeeffe(3,:), 'LineWidth',2);
xlabel('Input value');
ylabel('1st derivative');
mysavefig(h, filename, outdir, 12, [3, 4], 2);