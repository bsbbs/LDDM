% Excitatory neurons, tuning to the inputs, binary as an example
%% setup
gnrloutdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/General';
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionSngl2';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
Setup;





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