% Trajectory of trinary choice
%% SNR over training
%% Signal-to-Noise ratio over training 
N = 3;
cp = 6.4/100;
b = eye(N)*b0;
a = eye(N)*a0;
predur = 0;
presentt = 0;
stimdur = 5;
triggert = 0;
thresh = 70;
dur = 5; % second
stoprule = 1;
sgmG = 0;
BR = 0;
BG = 0;
v3vec = 0:0.05:1;
names = ["Early noise", "Late noise"];
mycols = [1,0,0; 0,0,1];
CellTypes = ["SST", "PV"]; 
%% example dynamics
sgmInput = 66*0.62;
sgmR = 0;
Vprior = zeros(1,N);
% winput = squeeze(winputc(end,:));
% wrg = squeeze(wrgc(end,:,:));
% wgr = squeeze(wgrc(end,:));
winput = ones(1,N)*3.1250;
wrg = ones(N,N)*.4419;
wgr = ones(1,N); %*.1927;

R0 = 20*ones(1,N);
D0 = 0*R0;
G0 = (wrg*R0')';
initialvals = [R0; G0; D0];
cp = 0.06;
Vinput = 100*[1+cp, 1-cp, 1];

[choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv3(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
    sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
h = figure; hold on;
filename = "ExmplDynmcs_TrinayCircuit";
plot(R(:,1));
plot(R(:,2));
plot(R(:,3));
xlabel('Time (ms)');
ylabel('Neural Activity(Option 2)');
%xlim([0, thresh]);
ylim([0, thresh]);
%% 
% --- Modified Plotting Section for Plane Projection ---
h = figure;
filename = "ExmplDynmcs_TrinayCircuit_Projected";

% 1. Define the transformation to 2D triangular coordinates
% x_2d = (R2 - R1) / sqrt(2)
% y_2d = (2*R3 - R1 - R2) / sqrt(6)
R1 = R(:,1); R2 = R(:,2); R3 = R(:,3);
x_proj = (R2 - R1) / sqrt(2);
y_proj = (2*R3 - R1 - R2) / sqrt(6);

% 2. Plot the projected trajectory
plot(x_proj, y_proj, 'LineWidth', 1.5, 'Color', [0.2, 0.2, 0.2]);
hold on;

% 3. Draw the decision boundary triangle (where R1+R2+R3 = thresh)
% Corners of the triangle in 2D
c1 = [-thresh/sqrt(2), -thresh/sqrt(6)]; % Option 1 corner (thresh, 0, 0)
c2 = [thresh/sqrt(2), -thresh/sqrt(6)];  % Option 2 corner (0, thresh, 0)
c3 = [0, 2*thresh/sqrt(6)];              % Option 3 corner (0, 0, thresh)
triangle = [c1; c2; c3; c1];
plot(triangle(:,1), triangle(:,2), 'k--', 'LineWidth', 1);

% 4. Add labels for the options at the corners
text(c1(1), c1(2)-5, 'Option 1', 'HorizontalAlignment', 'center');
text(c2(1), c2(2)-5, 'Option 2', 'HorizontalAlignment', 'center');
text(c3(1), c3(2)+5, 'Option 3', 'HorizontalAlignment', 'center');

% 5. Formatting
axis equal;
grid on;
title('Neural Trajectory Projected onto Decision Plane');
xlabel('Relative Evidence (R2 vs R1)');
ylabel('Relative Evidence (R3 vs others)');
% --- End of Modified Section ---

%% 3D to 2D projection
h = figure;
filename = "ExmplDynmcs_TrinayCircuit_Projected";

% 1. 计算轨迹到平面 R1 + R2 + R3 = thresh 的正交投影
% 平面的法向量为 [1, 1, 1]。对于空间中的点 R，其到平面的投影 R_proj 为：
% R_proj = R - k * [1, 1, 1]，其中 k = (R1 + R2 + R3 - thresh) / 3
k = (sum(R, 2) - thresh) / 3; 
R_proj = R - repmat(k, 1, 3); % 确保矩阵维度匹配相减

% 2. 绘制原始的 3D 轨迹 (用浅灰色虚线作为参考)
plot3(R(:,1), R(:,2), R(:,3), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;

% 3. 绘制投影到决策平面上的 3D 轨迹
plot3(R_proj(:,1), R_proj(:,2), R_proj(:,3), 'b-', 'LineWidth', 2);

% 4. 绘制切割 70Hz 阈值的决策平面（三角形）
corners = [thresh, 0, 0; 
           0, thresh, 0; 
           0, 0, thresh; 
           thresh, 0, 0];
% 画出三角形的边
plot3(corners(:,1), corners(:,2), corners(:,3), 'k-', 'LineWidth', 1.2);
% 给平面加上半透明的填充色以便于观察
patch('Vertices', corners(1:3,:), 'Faces', [1 2 3], ...
      'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% 5. 坐标轴及视角设置
xlabel('Neural Activity (Option 1)');
ylabel('Neural Activity (Option 2)');
zlabel('Neural Activity (Option 3)');
xlim([0, thresh]);
ylim([0, thresh]);
zlim([0, thresh]);
grid on;

% 调整视角，使其能够正视该投影平面
view(135, 35);
%%
for type = 2
    task = sprintf('SNR_V3_N%i_%s',N, CellTypes(type));
    h = figure;
    set(h, 'Units', 'inches','Position', [1,1,4,1.8]);
    filename = task;
    ldg = [];
                sgmInput = 66*0.62;
                sgmR = 0;
        Vprior = zeros(1,N);
        % winput = squeeze(winputc(end,:));
        % wrg = squeeze(wrgc(end,:,:));
        % wgr = squeeze(wgrc(end,:));
        winput = ones(1,N)*3.1250;
        wrg = ones(N,N)*.4419;
        wgr = ones(1,N)*.1927;

        R0 = 20*ones(1,N);
        D0 = 0*R0;
        G0 = (wrg*R0')';
        initialvals = [R0; G0; D0];
        rept = 40; %10;
        simfile = fullfile(Simdir, sprintf('%s_%s_%2.1f_%2.1f_%i.mat', task, names{testi}, sgmInput, sgmR, rept));
        if ~exist(simfile, 'file')
            MU = nan(length(v3vec),2);
            COV = nan(length(v3vec),3);
            SNR = nan(length(v3vec),2);
            r = nan(length(v3vec), 1);
            fprintf('%s',names{testi});
            for v3i = 1:length(v3vec)
                fprintf('.');
                Vinput = 100*[1+cp, 1-cp, v3vec(v3i)];
                
                mu = [];
                Sigma2 = [];
                parfor i = 1:rept
                    data = [];
                    [choice, rt, R, G, D, Vcourse] = LDDM_RndInputrv3(Vprior, Vinput, BR, BG, winput, wrg, wgr, a, b,...
                        sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
                    data(:,1) = R(10000:end,1);
                    data(:,2) = R(10000:end,2);
                    mu(i,:) = mean(data);
                    Sigma2(i,:,:) = cov(data);
                end
                % Sigma2 = cov(data);
                MU(v3i,:) = mean(mu, 1);
                % COV(v3i,:) = [Sigma2(1,1), Sigma2(2,2), Sigma2(1,2)];
                COV(v3i,:) = [mean(Sigma2(:,1,1)), mean(Sigma2(:,2,2)), mean(Sigma2(:,1,2))];
                SNR(v3i,:) = MU(v3i,:)./([COV(v3i,1), COV(v3i,2)]);
                r(v3i) = COV(v3i,3)/sqrt(COV(v3i,1)*COV(v3i,2));
            end
            fprintf('\n');
            dp = abs(MU(:,1) -  MU(:,2))./(COV(:,1) + COV(:,2) - 2*COV(:,3));
            save(simfile, 'MU','COV','SNR','r', 'dp');
        else
            load(simfile);
        end
        subplot(1,2,1); hold on;
        plot(v3vec, SNR(:,1), '-', 'Color', mycols(testi,:), 'LineWidth', lwd);
        plot(v3vec, SNR(:,2), '--', 'Color', mycols(testi,:), 'LineWidth', lwd);
        ylabel('SNR (V1 & V2)');
        xlabel('V3');

        subplot(1,2,2); hold on;
        ldg(testi) = plot(v3vec, dp, '-', 'Color', mycols(testi,:), 'LineWidth', lwd);
        ylabel('Discriminability (V1 & V2)');
        xlabel('V3');

    legend(ldg, {'Early noise', 'Late noise'}, "Location","east");
    mysavefig(h, filename, plotdir, fontsize, [5, 1.8]);
    % exportgraphics(h, fullfile(plotdir, [filename, '.pdf']), 'ContentType','vector');
end