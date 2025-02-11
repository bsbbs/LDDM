%%
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersion';
addpath('./utils/');
%% Pairing the plasticity
% according to Feild,... Freomke, 2020, Neuron, Figure 4
% The time constants of homosynaptic plasticity are 1.0 min and 0.6 min for excitation and inhibition, respoectively
% The time constants of heterosynaptic plasticity are 4.8 and 9.1 min for excitation and inhibition, respoectively
dt = .1;
duration = 30*60; %s
twindow = 60;
gamma = .001;
tauhom_e = 1*60;
tauhom_i = .6*60;
tauhet_e = 4.8*60;
tauhet_i = 9.1*60;

tauH = 100*60;
lambda = 0.1;
h = figure;
filename = 'Plasticity_dynamic_simulation';
plti = 0;
for type = ['e','i']
    for pair = [1, 0]
        switch pair
            case 1
                mycol = 'r';
            case 0
                mycol = 'b';
        end
        plti = plti + 1;
        subplot(2,2,plti);
        hold on;
        H = 0; % plasticity mediator with delayed response to pre-post pairing
        P0 = .5;
        P = P0; % the connection plasticity from I to E
        output = zeros(ceil(duration/dt),3);
        x = -ceil((twindow+3)/dt):ceil(duration/dt);
        for i = 1:length(x)
            PrePost = 100 *(x(i)<=0 & x(i) >= -twindow/dt) * pair;
            dH = (-H/tauH + PrePost*lambda)*dt;
            H = H + dH;
            dP = (-P/eval(['tauhet_' type]) + (1-P)*gamma*H/eval(['tauhom_' type]))*dt;
            P = P + dP;
            output(i,:) = [PrePost, H, P];
        end
        plot(x, output(:,1)./max(output(:,1)) + 1, 'k-', 'LineWidth', 1);
        plot(x, output(:,2)./max(output(:,2)) + 1, 'k--', 'LineWidth', 1);
        plot(x, output(:,3)/P0, '-', 'Color', mycol, 'LineWidth', 2);
        xticks([0:5:30]*60/dt);
        xticklabels(num2str([0:5:30]'));
        xlabel('Time (min)');
        ylabel('Strength rescaled');
        if plti == 2
            legend({'Pre-Post pairing','Hidden process','Synaptic strength'}, 'Location', 'Best');
        end
    end
end
mysavefig(h, filename, outdir, 14, [8, 6]);
%% Dynamics of the synaptic connections
dt = .01;
tau = .3;
duration = 20; %s
twindow = 3;
gamma = .001;
tauR = .1;
tauG = .1;
tauD = .1;
tauP = 60;
tauH = 60;
w = 1;
R = 0;
G = 0;
D = 0;
H = 0; % plasticity mediator with delayed response to pre-post pairing
P = 0; % the connection plasticity from I to E

Input = 20;
output = zeros(ceil(duration/dt),5);
for ti = 1:ceil(duration/dt)
    dR = (-R + Input/(1 + P*G))*dt/tauR;
    dG = (-G + w*R)*dt/tauG;
    G = G + dG;
    R = R + dR;

    PrePost = R*G *(ti<twindow/dt);
    dH = (-H + PrePost)*dt/tauH;
    H = H + dH;
    dP = (-P + (1-P)*gamma*H)*dt/tauP;
    P = P + dP;
    output(ti,:) = [R, G, PrePost, H, P];
end

plot(output./max(output), 'LineWidth', 2);
legend({'R','G','PrePost','H','P'});

%% The network architecture
rng(2023);
r_e = 800; % radius of the pyramidal axons, um
r_i = 250;% radius of the inhibitory axons, um
h = figure;
filename = 'Schema_Before';
hold on;
D = 2000; % physical distance between the two pools of pyramidal neurons, um
x = [-D, D]/2;
y = [0, 0];
plot(x, y, '.k', 'MarkerSize', 100);
N = 2000;
xI = -D + rand(N,1)*D*2;
yI = -D/2 + rand(N,1)*D;
plot(xI, yI, '.b', 'MarkerSize',1);
Dstc = []; % distance between the pyramidal pool and the interneurons
PhysclCnnctProb_EI = []; % physical connection probability from the pyramidal
% pool to the interneurons, exponentially decay according to r_e
PhysclCnnct_EI = []; % physical connection from the pyramidal pool to the interneurons
PhysclCnnctProb_IE = []; % physical connection probability from the interneurons
% to the pyramidal pool, exponentially decay according to r_e
PhysclCnnct_IE = []; % physical connection from the interneurons to the pyramidal pool
for pyri = 1:2
    Dstc(:,pyri) = sqrt((xI - x(pyri)).^2 + (yI - y(pyri)).^2);
    PhysclCnnctProb_EI(:, pyri) = exp(-Dstc(:,pyri)/r_e);
    PhysclCnnctProb_IE(:, pyri) = exp(-Dstc(:,pyri)/r_i);
end
PhysclCnnct_EI = PhysclCnnctProb_EI >= rand(size(PhysclCnnctProb_EI));
PhysclCnnct_IE = PhysclCnnctProb_IE >= rand(size(PhysclCnnctProb_IE));
% for pyri = 1:2
%     for j = 1:numel(xI)
%         if PhysclCnnct_EI(j, pyri) == 1
%             plot([x(pyri), xI(j)], [y(pyri), yI(j)], 'k-', 'LineWidth', 1);
%         end
%     end
% end
mask1 = any(PhysclCnnct_EI, 2);
plot(xI(mask1), yI(mask1),'.b', 'MarkerSize',15);
plot(xI(~mask1), yI(~mask1),'.', 'Color', [.5,.5,.5], 'MarkerSize',15);
for j = 1:numel(xI)
    for pyri = 1:2
        if PhysclCnnct_IE(j, pyri) == 1 && PhysclCnnct_EI(j, pyri) == 1
            plot([xI(j), x(pyri)], [yI(j), y(pyri)], 'r-',  'LineWidth', 2);
        end
    end
end
mask2 = any(PhysclCnnct_IE, 2);
mask3 = all(PhysclCnnct_EI, 2);
plot(xI(mask2&mask3), yI(mask2&mask3),'.c', 'MarkerSize',25);
for j = 1:numel(xI)
    for pyri = 1:2
        if PhysclCnnct_IE(j, pyri) == 1 && all(PhysclCnnct_EI(j, :) == 1)
            plot([xI(j), x(pyri)], [yI(j), y(pyri)], 'r-',  'LineWidth', 2);
            % plot([x(3-pyri), xI(j)], [y(3-pyri), yI(j)], 'c-', 'LineWidth', 1);
            plot([x(3-pyri), xI(j)], [y(3-pyri), yI(j)], 'c-', 'LineWidth', 1);
        end
    end
end
xlabel('\mum');
ylabel('\mum');
mysavefig(h, filename, outdir, 14, [8, 4]);

%% STDP - singleton
Input = [20, 20];
R = 0;
G = zeros(size(PhysclCnnct_IE(:,1)));
H = G;
P = G;

for duration = [20, 30*60]
    switch duration
        case 20
            filename = 'Shortterm_dynamic';
            dt = 0.01;
            xt = 0:5:duration;
            xtl = xt;
            unit = 'Time (s)';
        case 30*60
            filename = "Longterm_dynamic";
            dt = 0.1;
            xt = 0:5*60:duration;
            xtl = xt/60;
            unit = 'Time (min)';
    end
    time = 0:ceil(duration/dt);
    Rt = [];
    Gt = [];
    Ht = [];
    Pt = [];
    type = 'i';
    for i = 1:length(time)
        G0 = G;
        dR = (-R + Input(1)/(1 + sum(P.*G)))*dt/tauR;
        dG = (-G + w*R*PhysclCnnct_EI(:,1))*dt/tauG;
        R = R + dR;
        G = G + dG;
        PrePost = G0*R;
        dH = (-H/tauH + lambda*PrePost.*PhysclCnnct_IE(:,1))*dt;
        H = H + dH;
        dP = (-P/eval(['tauhet_' type]) + gamma*(1-P).*H/eval(['tauhom_' type]))*dt;
        P = P + dP;
        Rt(:,i) = R;
        Gt(:,i) = G;
        Ht(:,i) = H;
        Pt(:,i) = P;
    end
    % visualization
    h = figure;
    subplot(2,1,1);
    hold on;
    plot(Rt, 'k-', 'LineWidth', 2);
    plot(mean(Gt(~PhysclCnnct_EI(:,1),:), 1), '-', 'Color', [.5, .5, .5], 'LineWidth', 2);
    plot(mean(Gt(PhysclCnnct_EI(:,1),:), 1), 'b-', 'LineWidth', 1);
    plot(mean(Gt(PhysclCnnct_EI(:,1) & PhysclCnnct_IE(:,1),:), 1), 'r--', 'LineWidth', 1);
    legend({'R','G-unconnected','G-connected', 'G-reciprocal'}, "Location", "best");
    xticks(xt/dt);
    xticklabels(num2str(xtl'));
    xlabel(unit);
    ylabel('Firing rates (Hz)');
    mysavefig(h, filename, outdir, 14, [6,4]);

    subplot(2,1,2);
    hold on;
    plot(mean(Pt(~PhysclCnnct_EI(:,1) & PhysclCnnct_IE(:,1),:), 1), 'b-', 'LineWidth', 2);
    plot(mean(Pt(PhysclCnnct_EI(:,1) & PhysclCnnct_IE(:,1),:), 1), 'r-', 'LineWidth', 2);
    legend({'Synaptic strength - uni-directional','Synaptic strength - reciprocal'}, "Location", "best");
    xticks(xt/dt);
    xticklabels(num2str(xtl'));
    xlabel(unit);
    ylabel('Synaptic strength');
    mysavefig(h, filename, outdir, 14, [6,8]);
end


%% temporal_kernal(Pre, Post)
function dW = temporal_kernal(Pre_history, Post_history)
tau_plus = 5; %ms
tau_minus = 5; %ms
% granger causality

% adaptation on the connection weigths, i.e., dW

end

%% Orenstein-Ulenbeck process
function noise = OU(noise, sgm, dt, size) 
tauN = .002; % s
noise = noise + (-noise + randn(size).*sqrt(dt).*sgm)/tauN*dt;
end