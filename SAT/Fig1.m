% Fig 2 - Example dynamics of different models, and their prediction about
% the shape of reaction time



%% panel b, Example dynamics of the LDDM
nColors = 7;
cols_warm = [ones(nColors, 1), linspace(0,1,nColors)'.*[1,1]]; %autumn(7);
cols_cool = [linspace(0,1,nColors)'.*[1,1], ones(nColors, 1)]; %winter(7);
cols_gray = linspace(0,1,nColors)'.*[1,1, 1];

w = ones(2);
cp = .128; %.512;
initialvals = [10,10;20,20;0,0]*2;
predur = 0;
presentt = 0;
dur = 1.5;
stimdur = dur;
sgm = 3;
sgmInput = 1*scale0;
sgmtest = 4;
triggert = 0;
thresh = 70;
stoprule = 1;

h = figure;
filename = 'Fig2b';
Vprior = [1, 1]*0;
Vinput = [1 + cp, 1 - cp]*scale0;
sets = [40, .7; % adjust self-excitation
    80, .7;
    0, 1.265; % adjust disinhibition
    0, 1.76];
for i = 1:2
    if i == 1
        opt2thick = cols_gray(2,:);
        opt2thin = cols_gray(5,:);
        opt1thin = cols_warm(5,:);
        opt1thick = cols_warm(1,:);
    elseif i == 2
        opt2thick = cols_gray(2,:);
        opt2thin = cols_gray(5,:);
        opt1thin = cols_cool(5,:);
        opt1thick = cols_cool(1,:);
    end
    subplot(1,2,i); hold on;
    for j = 1:2
        si = 2*(i-1)+j;
        a = sets(si,1)*eye(2);
        b = sets(si,2)*eye(2);
        % switch si
        %     case 1
        %         mycol = cols_warm(1,:);
        %         colgray = cols_gray(1,:);
        %     case 2
        %         mycol = cols_warm(1,:);
        %         colgray = cols_gray(1,:);
        %     case 3
        %         mycol = cols_cool(1,:);
        %         colgray = cols_gray(1,:);
        %     case 4
        %         mycol = cols_cool(1,:);
        %         colgray = cols_gray(1,:);
        % end
        % colgray = cols_gray(5,:);
        % [~, rtbar, Rbar, Gbar, Dbar] = LDDM(Vprior, Vinput, w, a, b,...
        %     sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        [~, rtbar, Rbar, Gbar, Dbar, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
            0, 0, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        for ti = 1:10
            % [choice, rt, R, G, D] = LDDM(Vprior, Vinput, w, a, b,...
            %     sgmtest, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            [choice, rt, R, G, D, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
                sgm, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            plot(R(:,2), '--', 'Color', opt2thin, 'LineWidth', .5);
            plot(R(:,1), '-', 'Color', opt1thin, 'LineWidth', .5);
        end
        lgd2 = plot(Rbar(:,2), '--','Color', opt2thick, 'LineWidth',lwd);
        lgd1 = plot(Rbar(:,1), '-','Color', opt1thick, 'LineWidth',lwd);
        plot((predur + triggert + rtbar + [-.14, .14])/dt,  [thresh, thresh],'k-','LineWidth',.5);
        %plot(diff(Rbar(:,1)), '-','Color', mycol, 'LineWidth',lwd);
    end
    xlabel('Time (ms)');
    ylim([0, 70]);
    yticks([0, 70]);
    yticklabels({'0','70'});
    ylabel('Firing rates (Hz)');
    legend([lgd1, lgd2],{'Option 1','Option 2'},...
        'Location','Best','FontSize',fontsize-5, 'FontName','Times New Roman', ...
        'FontAngle','italic','Box','off');
    % ylabel('Change rate (Hz/ms)');
    mysavefig(h, filename, plotdir, fontsize, [6, 3]);
end

%% Simulating RT distributions
w = ones(2);
cp = .128; %.512;
Vprior = [1, 1]*0;
Vinput = [1 + cp, 1 - cp]*scale0;
initialvals = [10,10;20,20;0,0]*2;
predur = 0;
presentt = 0;
dur = 5;
stimdur = dur;
sgm = 3;
sgmInput = 1*scale0;
sgmtest = 4;
triggert = 0;
thresh = 70;
stoprule = 1;
sims = 2048;

[rt, choice, ~] = LDDM_Rndinput_GPU(Vprior, Vinput, w, a, b,...
    sgm, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);