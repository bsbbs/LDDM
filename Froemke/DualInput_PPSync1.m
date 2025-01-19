% Excitatory neurons, tuning to the inputs, binary as an example
%% setup
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionDual';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
Setup;
%% The structure of the network
Ntwk.Scale = 4000; % scale of the micro structure to test, 2D square (length and width) in unit of um
Ntwk.Input.N = 2; % number of inputs 
Inputx = linspace(-.5, .5, Ntwk.Input.N+2);
Ntwk.Input.Location = Ntwk.Scale*[Inputx(2:end-1)', ones(Ntwk.Input.N, 1)*.25]; % the location of the input intervene
Ntwk.Input.Variance = 200; % the range of each input intervening
Ntwk.Exct.Location = []; % the physical location of excitatory cells
Ntwk.Exct.N = 200; % number of the excitarory cells
Ntwk.Exct.Properties.size = 15; % physical size of the soma, in um
Ntwk.Exct.Properties.FWHM = 1500; % the Full width at half maximum(FWHM) range of axons and dendrites connection
Ntwk.Inhbt.Location = [];
Ntwk.Inhbt.N = 50;
Ntwk.Inhbt.Properties.size = 8;
Ntwk.Inhbt.Properties.FWHM = 500;
Ntwk.Noise.tauN = .002; % time constant for the OU process of noise
Ntwk.Noise.sgm = 10; % amplitude of noise
Ntwk.tauR = .1853; % time constant for the three types of neurons
Ntwk.tauG = .2244;
Ntwk.tauD = .3231;
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
p_EI = 1 - (1 - exp(-DstcEI/Ntwk.Exct.Properties.FWHM)).^4; % probability of connection from E to I based on distance
% connection from I -> E
[XI, XE] = meshgrid(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent Exct and columns represent Inhbt
[YI, YE] = meshgrid(Ntwk.Inhbt.Location(:,2), Ntwk.Exct.Location(:,2));
DstcIE = sqrt((XI - XE).^2 + (YI - YE).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
p_IE = 1 - (1 - exp(-DstcIE/Ntwk.Inhbt.Properties.FWHM)).^4; % probability of connection from I to E based on distance
% connection from E -> E
[XE1, XE2] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1));
[YE1, YE2] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Exct.Location(:,2));
DstcEE = sqrt((XE1 - XE2).^2 + (YE1 - YE2).^2); %  Euclidean distance between each pair of Exct neurons
p_EE = 1 - (1 - exp(-DstcEE/Ntwk.Exct.Properties.FWHM)).^4; % probability of connection from E to E based on distance

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
h = figure;
filename = 'pEE';
imagesc(p_EE');
caxis([0, 1]);
c = colorbar;
c.Label.String = 'Probability of physical connections';
% colormap(bluewhitered(256)), colorbar;
xlabel('Exct channels');
ylabel('Exct channels');
mysavefig(h, filename, outdir, 12, [4, 3]);

h = figure;
filename = 'p_EI';
imagesc(p_EI');
caxis([0, 1]);
%c = colorbar;
%c.Label.String = 'Probability of physical connections';
% colormap(bluewhitered(256)), colorbar;
xlabel('Inhbt channels');
ylabel('Exct channels');
mysavefig(h, filename, outdir, 12, [1.6, 3]);

h = figure;
filename = 'p_IE';
imagesc(p_IE');
caxis([0, 1]);
%c = colorbar;
%c.Label.String = 'Probability of physical connections';
% colormap(bluewhitered(256)), colorbar;
ylabel('Inhbt channels');
xlabel('Exct channels');
mysavefig(h, filename, outdir, 12, [3.22, 1.46]);

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
    pseudolocation = linspace(-Ntwk.Scale/2, Ntwk.Scale/2, 100);
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
            plot([Ntwk.Input.Location(i,1), Ntwk.Exct.Location(j,1)], ...
                [Ntwk.Input.Location(i,2), Ntwk.Exct.Location(j,2)],...
                '-', 'Color', OKeeffe(i,:), 'LineWidth', Ntwk.Exct.tuning(j,i)/max(Ntwk.Exct.tuning(:))*1);
        end
    end
    xlabel('\mum');
    ylabel('\mum');
    xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
    ylim([-Ntwk.Scale/4, Ntwk.Scale/4]);
    mysavefig(h, filename, outdir, 12, [3, 3], 2);
end

%% initial I->E synaptic weights
rng(2023);
wIE_initial = .01*randn(Ntwk.Exct.N, Ntwk.Inhbt.N).*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST 
wIE_initial(wIE_initial<0) = 0;
s = wIE_initial;
if Ntwk.show
    % synaptic connection weights
    h = figure;
    filename = sprintf('SynapticW_at_time%1.0fs', 0);
    imagesc(wIE_initial');
    caxis([-.1, .1]);
    colormap(bluewhitered(256));
    c = colorbar;
    c.Label.String = 'Connection weights';
    ylabel('Inhbt channels');
    xlabel('Exct channels');
    mysavefig(h, filename, outdir, 12, [3.928, 1.46]);
end

%% Input values and seqeuences
dt = .001; % time precision for simulation, in unit of second
Ntrial = 258;
rng(2023);
Seq = CreateEvents(Ntrial, dt);
Seq(:,2) = Seq(:, 1);
value = 20*ones(1, Ntwk.Input.N);
time = [0:(size(Seq, 1)-1)]*dt/60; % unit in Mins
% time = [0:5000*60]*dt/60;
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
ylabel('Input');
xlabel('Time (mins)');
mysavefig(h, filename, outdir, 12, [3, 11], 2);

%% simulation begin ... excitatory activities with Background white noise (OU process) input to the excitatory cells
outdir = fullfile(outdir, 'PPSync');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
R = zeros(Ntwk.Exct.N,numel(time));
G = zeros(Ntwk.Inhbt.N,numel(time));
noise = zeros(Ntwk.Exct.N, 1);
sdynamc = s(98,32);
a_pre = zeros(1, Ntwk.Inhbt.N);
a_post = zeros(Ntwk.Exct.N, 1);
for ti = 1:(length(time)-1)
    if (mod(ti*dt, 1)==0)
        fprintf('%is.',ti*dt*1);
    end
    if (mod(ti*dt, 60)==0)
        fprintf('\n');
    end
    % update R and G activities
    dR = (-R(:,ti) + (Ntwk.alpha*Ntwk.Cnnct_EE*R(:,ti) + Ntwk.Exct.tuning*(Seq(ti,:).*value)' + noise)./(1 + s*G(:,ti)))*dt/Ntwk.tauR;
    dG = (-G(:,ti) + Ntwk.w*Ntwk.Cnnct_EI*R(:,ti))*dt/Ntwk.tauG;
    R(:,ti+1) = R(:,ti) + dR;
    G(:,ti+1) = G(:,ti) + dG;
    R(R(:, ti+1) < 0, ti+1) = 0;
    G(G(:, ti+1) < 0, ti+1) = 0;
    % update OU noise
    noise = noise + (-noise + randn(size(noise))*sqrt(dt)*Ntwk.Noise.sgm)*dt/Ntwk.Noise.tauN;
    % apply STDP kernel to the course  of firing rates
    preAP = G(:,ti)';
    postAP = R(:,ti);
    a_pre = a_pre + (-a_pre/Ntwk.tau_prepost)*dt + Ntwk.A*preAP*dt; % the intermiediate variables for pre
    a_post = a_post + (-a_post/Ntwk.tau_postpre)*dt + Ntwk.A*postAP*dt; % and post action potential induced plasticity
    a = postAP*dt*a_pre + a_post*preAP*dt; % the AP induced overall plasticity change
    % update on the synaptic weights
    ds = (-s/Ntwk.taus)*dt + Ntwk.gamma*(1 - s).*a.*Ntwk.Cnnct_IE; % a already contains the dt information, only applied to physical connections
    s = s + ds;
    s(s<0) = 0;
    sdynamc(ti+1) = s(98,32);
end
fprintf('\n');
wIE = s;
Ntwk.wIE = s;
save(fullfile(outdir, 'Simulation_full.mat'), 'Ntwk', 's', 'sdynamc', 'wIE_initial', 'wIE', 'Seq');
% synapctic connection weights
h = figure;
filename = sprintf('SynapticW_at_time%1.0fs', ti*dt);
imagesc(wIE');
caxis([-.1, .1]);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Connection weights';
ylabel('Inhbt channels');
xlabel('Exct channels');
mysavefig(h, filename, outdir, 12, [3.928, 1.46]);
% synaptic weight change
h = figure;
filename = sprintf('SynapticWChange_at_time%1.0fs', ti*dt);
imagesc((wIE - wIE_initial)');
caxis([-.1, .1]);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Weights change';
ylabel('Inhbt channels');
xlabel('Exct channels');
mysavefig(h, filename, outdir, 12, [3.928, 1.46]);
h = figure;
plot(time, sdynamc, 'k-', 'LineWidth',1);

%% visualization
if Ntwk.show
    % R
    h = figure;
    filename = sprintf('R_Dynamic%1.0fs', max(smpl_time));
    hold on;
    Amp = max(max(R(:,1:numel(smpl_time))));
    for i = 1:Ntwk.Exct.N
        plot(smpl_time, R(i,1:numel(smpl_time))+i*Amp, 'k-');
    end
    ylim([0, (i+1)*Amp]);
    yticks([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Exct.N/5):round(Ntwk.Exct.N/5):Ntwk.Exct.N]'));
    ylabel('Exct channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [3, 11.2]);
    % G
    h = figure;
    filename = sprintf('G_Dynamic%1.0fs', max(smpl_time));
    hold on;
    Amp = max(max(G(:,1:numel(smpl_time))));
    for i = 1:Ntwk.Inhbt.N
        plot(smpl_time, G(i,1:numel(smpl_time))+i*Amp, 'k-');
    end
    ylim([0, (i+1)*Amp]);
    yticks([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]*Amp);
    yticklabels(num2str([1,round(Ntwk.Inhbt.N/5):round(Ntwk.Inhbt.N/5):Ntwk.Inhbt.N]'));
    ylabel('Inhbt channels');
    xlabel('Time (secs)');
    mysavefig(h, filename, outdir, 12, [3, 3.5]);
    
    % Example dynamics R and G
    h = figure;
    maxValue = max(wIE(:));
    [IdxE, IdxI] = find(wIE == maxValue);
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fs',IdxE(1), IdxI(1), max(smpl_time));
    hold on;
    % plot(smpl_time, R(IdxE(1),1:numel(smpl_time)), '-', "Color", OKeeffe(3,:), 'LineWidth',1);
    % plot(smpl_time, G(IdxI(1),1:numel(smpl_time)), '-', "Color", OKeeffe(4,:), 'LineWidth',1);
    segment = round(smpl_time/dt);
    plot(smpl_time, R(IdxE(1),segment)/max(R(IdxE(1),segment)), '-', "Color", OKeeffe(3,:), 'LineWidth',.5);
    plot(smpl_time, G(IdxI(1),segment)/max(G(IdxI(1),segment)), '-', "Color", OKeeffe(4,:), 'LineWidth',.5);
    xlabel('Time (secs)');
    ylabel('Normalized activity');
    %ylabel('Firing rates (Hz)');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);

    h = figure;
    maxValue = max(wIE(:));
    [IdxE, IdxI] = find(wIE == maxValue);
    filename = sprintf('RGExmplpair_at_E%iandI%i_%1.0fmins',IdxE(1), IdxI(1), max(time));
    hold on;
    plot(time, R(IdxE(1),1:numel(time)), '-', "Color", OKeeffe(3,:), 'LineWidth', .5);
    plot(time, G(IdxI(1),1:numel(time)), '-', "Color", OKeeffe(4,:), 'LineWidth', .5);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [26, 26, 2, 2, 26], '-k');
    xlim([0, 35]);
    ylabel('Firing rates (Hz)');
    xlabel('Time (mins)');
    legend({"R", "G"}, 'Location', "best", 'FontSize', 10);
    mysavefig(h, filename, outdir, 12, [3, 2], 2);

    h = figure; hold on;
    filename = 's_dynamic_seg330';
    plot(smpl_time, sdynamc(segment), 'k-', 'LineWidth',1);
    xlabel('Time (secs)');
    ylabel('Synaptic weight');
    mysavefig(h, filename, outdir, 12, [1.8, 1.2], 2);

    h = figure; hold on;
    filename = 's_dynamic';
    plot(time, sdynamc, 'k-', 'LineWidth',1);
    plot([min(smpl_time), max(smpl_time), max(smpl_time), min(smpl_time), min(smpl_time)]/60, [4.2, 4.2, 3.0, 3.0, 4.2]/100, '-k');
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
    mysavefig(h, filename, outdir, 12, [5, 1]);
    
    h = figure;
    filename = 'Tuning_Inhbt';
    Ntwk.Inhbt.tuning = Ntwk.Cnnct_EI * Ntwk.Exct.tuning;
    imagesc(Ntwk.Inhbt.tuning');
    caxis([0, max(Ntwk.Inhbt.tuning(:))]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Inhbt channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Tuning';
    mysavefig(h, filename, outdir, 12, [5, 1]);
    
    h = figure;
    filename = 'Tuning_Inhbt_IE';
    Ntwk.Inhbt.tuningIE = Ntwk.wIE' * Ntwk.Exct.tuning;
    imagesc(Ntwk.Inhbt.tuningIE');
    caxis([0, max(Ntwk.Inhbt.tuningIE(:))]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Inhbt channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Selective Inhibition';
    mysavefig(h, filename, outdir, 12, [5, 1]);
    
    h = figure;
    filename = 'Tuning_Inhbt_IE_sort';
    [IE1, I] = sort(Ntwk.Inhbt.tuningIE(:,1), 'descend');
    IE2 = Ntwk.Inhbt.tuningIE(I,2);
    imagesc([IE1, IE2]');
    caxis([0, max(Ntwk.Inhbt.tuningIE(:))]);
    colormap(bluewhitered(256));
    c = colorbar;
    xlabel('Inhbt channels');
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    c.Label.String = 'Selective Inhibition';
    mysavefig(h, filename, outdir, 12, [5, 1]);

    % mean-field connectome
    InhbtMtrx(1,1) = sum((Ntwk.Inhbt.tuning(:,1)/sum(Ntwk.Inhbt.tuning(:,1), 1)).*(Ntwk.Inhbt.tuningIE(:,1)/sum(Ntwk.Inhbt.tuningIE(:,1),1)), 1);
    InhbtMtrx(1,2) = sum((Ntwk.Inhbt.tuning(:,1)/sum(Ntwk.Inhbt.tuning(:,1), 1)).*(Ntwk.Inhbt.tuningIE(:,2)/sum(Ntwk.Inhbt.tuningIE(:,2),1)), 1);
    InhbtMtrx(2,1) = sum((Ntwk.Inhbt.tuning(:,2)/sum(Ntwk.Inhbt.tuning(:,2), 1)).*(Ntwk.Inhbt.tuningIE(:,1)/sum(Ntwk.Inhbt.tuningIE(:,1), 1)), 1);
    InhbtMtrx(2,2) = sum((Ntwk.Inhbt.tuning(:,2)/sum(Ntwk.Inhbt.tuning(:,2), 1)) .*(Ntwk.Inhbt.tuningIE(:,2)/sum(Ntwk.Inhbt.tuningIE(:,2), 1)), 1);
    h = figure;
    filename = 'Meanfield-RcrtInhbt';
    imagesc(InhbtMtrx);
    for i = 1:2
        for j = 1:2
            text(i,j,num2str(InhbtMtrx(i,j)));
        end
    end
    caxis([0, max(InhbtMtrx(:))*1.2]);
    colormap('gray');
    %colormap(bluewhitered(556));
    c = colorbar;
    c.Label.String = 'Recurrent Inhibition';
    xticks([1,2]);
    xticklabels({'Input 1', 'Input 2'});
    yticks([1,2]);
    yticklabels({'Input 1', 'Input 2'});
    xlabel('E -> I');
    ylabel('I -> E');
    mysavefig(h, filename, outdir, 12, [5, 4]);
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