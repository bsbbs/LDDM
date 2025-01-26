% Fig 2 - Example dynamics of different models, and their prediction about
% the shape of reaction time
nColors = 7;
red = [1, 0, 0];
orange = [1, .5, 0]; %[255, 140, 0]/255;
cols_warm = [linspace(orange(1), red(1), nColors)', ...
    linspace(orange(2), red(2), nColors)', ...
    linspace(orange(3), red(3), nColors)'];
% [ones(nColors, 1), linspace(0,1,nColors)'.*[1,1]]; % graded red
cyan = [0, 1, 1];
blue = [0, 0, 1];
cols_cool = [linspace(cyan(1), blue(1), nColors)', ...
    linspace(cyan(2), blue(2), nColors)', ...
    linspace(cyan(3), blue(3), nColors)'];
% [linspace(0,1,nColors)'.*[1,1], ones(nColors, 1)]; % graded blue
cols_gray = linspace(0,1,nColors)'.*[1,1, 1];
%% panel b, Example dynamics of the LDDM
w = ones(2);
cp = .128; %.512;
initialvals = [1,1;2,2;0,0]*30;
predur = 0;
presentt = 0;
dur = 10;
stimdur = dur;
sgm = 2;
sgmInput = 7*scale0;
triggert = 0;
thresh = 70;
stoprule = 1;
sims = 102400;
h = figure;
plotname = 'RTChoice_Demo2b';
Vprior = [1, 1]*0;
Vinput = [1 + cp, 1 - cp]*scale0;
sets = [30, .7; % adjust self-excitation
    80, .7;
    0, 1.176; % adjust disinhibition
    0, 1.71];
opt2col = cols_gray(2,:);
for i = 1:2
    subplot(1,2,i); 
    hold on;
    for j = 1:2
        si = 2*(i-1)+j;
        if si == 1
            opt1col = orange;
        elseif si == 2
            opt1col = red;
        elseif si == 3
            opt1col = cyan;
        elseif si == 4
            opt1col = blue;
        end
        a = sets(si,1)*eye(2);
        b = sets(si,2)*eye(2);
        for ti = 1:5
            [choice, rt, R, G, D, ~] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
                sgm/2, sgmInput/3, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R(R>70) = 70;
            plot([0:length(R)-1]*dt, R(:,2), '--', 'Color', [mean(opt1col)*[1, 1, 1], .2], 'LineWidth', .5);
            plot([0:length(R)-1]*dt, R(:,1), '-', 'Color', [opt1col, .2], 'LineWidth', .5);
        end
        [~, rtbar, Rbar, Gbar, Dbar, Vcourse] = LDDM_RndInput(Vprior, Vinput, w, a, b,...
            0, 0, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        lgd2 = plot([0:length(Rbar)-1]*dt, Rbar(:,2), '--','Color', mean(opt1col)*[1, 1, 1], 'LineWidth',lwd);
        lgd1 = plot([0:length(Rbar)-1]*dt, Rbar(:,1), '-','Color', opt1col, 'LineWidth',lwd);
        filename = sprintf('RTChoice_a%2.1f_b%1.3f_sgmInput%1.1f_sgm%1.1f_cp%0.3f', a(1), b(1), sgmInput/scale0, sgm, cp);
        simmat = fullfile(datadir, [filename, '.mat']);
        if ~exist(simmat,"file")
            [rt, choice, ~] = LDDM_Rndinput_GPU(Vprior, Vinput, w, a, b,...
                sgm, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            rt = squeeze(rt);
            choice = squeeze(choice);
            save(simmat, 'rt','choice', '-mat');
        else
            load(simmat);
        end
        nanmask = ~isnan(rt);
        rt = rt(nanmask);
        choice = choice(nanmask);
        [kde, x] = ksdensity(rt, "Kernel", 'normal', 'NumPoints', 100);
        kde_smth = smooth(x, kde, .07, 'loess');
        kde_smth(kde_smth<0) = 0;
        rate = 10*(j^.5)/max(kde_smth);
        plot(x, kde_smth*rate+71, '-','Color', opt1col, 'LineWidth',2);
        text(mean(rt), 71+max(kde_smth*rate), sprintf("%1.2f", mean(2 - choice)));
    end
    xlim([0,3]);
    plot(xlim,  [thresh, thresh],'k-','LineWidth',.5);
    xlabel('Time (s)');
    ylim([0, 90]);
    yticks([0, 70]);
    yticklabels({'0','70'});
    ylabel('Firing rates (Hz)');
    legend([lgd1, lgd2],{'Option 1','Option 2'},...
        'Location','Best','FontSize',fontsize-5, 'FontName','Arial', ...
        'FontAngle','italic','Box','off');
    mysavefig(h, plotname, plotdir, fontsize, [6, 3]);
end

%% Simulating RT distributions
Rslts = table('Size', [0 8], 'VariableTypes', {'uint8', 'uint8', 'double', 'double', 'double', 'double', 'double', 'double'},...
    'VariableNames', {'testi', 'i', 'a', 'b', 'accuracy', 'rtmean', 'rtstd', 'rtskew'});
cp = .128;
Vinput = [1 + cp, 1 - cp]*scale0;
sgm = 2;
sgmInput = 7*scale0;
% adapting alpha
test = {'Adjust alpha', 'Adjust beta'};
valtext  = {'\alpha', '\beta'};
ntest = 7;
avec = linspace(30, 80, ntest);
bvec = linspace(1.176, 1.71, ntest);
ratios = nan(2, ntest);
meanrts = nan(2,ntest);
stdrts = nan(2,ntest);
skewrts = nan(2,ntest);
vals = nan(2,ntest);
for testi = 1:2
    for i = 1:ntest
        if testi == 1
            a = avec(i)*eye(2);
            b = .7*eye(2);
            cols = cols_warm;
            vals(testi,i) = avec(i);
        elseif testi == 2
            a = 0*eye(2);
            b = bvec(i)*eye(2);
            cols = cols_cool;
            vals(testi,i) = bvec(i);
        end
        filename = sprintf('RTChoice_a%2.1f_b%1.2f_sgmInput%1.1f_sgm%1.1f_cp%0.3f', a(1), b(1), sgmInput/scale0, sgm, cp);
        simmat = fullfile(datadir, [filename, '.mat']);
        if ~exist(simmat,"file")
            [rt, choice, ~] = LDDM_Rndinput_GPU(Vprior, Vinput, w, a, b,...
                sgm, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
            rt = squeeze(rt);
            choice = squeeze(choice);
            save(simmat, 'rt','choice', '-mat');
        else
            load(simmat);
        end
        nanmask = ~isnan(rt);
        rt = rt(nanmask);
        choice = choice(nanmask);
%         h = figure;
%         MyH = histogram(rt, 20, 'FaceColor',[.7, .7, .7], 'Normalization','pdf');
%         hold on;
%         [kde, x] = ksdensity(rt, "Kernel", 'normal');
%         plot(x, kde, '-', 'Color', cols(i,:), 'LineWidth',2);
%         xlabel('Reaction time (s)');
%         ylabel('Density');
%         title(sprintf('a = %2.1f, b = %1.2f', a(1), b(1)));
        ratios(testi,i) = mean(2-choice);
        meanrts(testi,i) = mean(rt);
        stdrts(testi,i) = std(rt);
        skewrts(testi,i) = skewness(rt);
        new_row = table(testi, i, a(1), b(1), mean(2-choice), mean(rt), std(rt), skewness(rt), 'VariableNames', Rslts.Properties.VariableNames);
        Rslts = [Rslts; new_row];
        writetable(Rslts, fullfile(datadir, 'SimRslts.txt'), 'Delimiter', '\t');
    end
    h = figure;
    filename = sprintf('RTChoice_%s', test{testi});
    subplot(1,4,1);
    scatter(vals(testi,:), ratios(testi,:)*100, 28, cols, 'filled');
    xlabel(valtext{testi});
    ylabel('Accuracy (%)');
    mysavefig(h, filename, plotdir, fontsize, [6, 1.8]);
    subplot(1,4,2);
    scatter(vals(testi,:), meanrts(testi,:), 28, cols, 'filled');
    xlabel(valtext{testi});
    ylabel('RT Mean (s)');
    mysavefig(h, filename, plotdir, fontsize, [8, 1.8]);
    subplot(1,4,3);
    scatter(vals(testi,:), stdrts(testi,:), 28, cols, 'filled');
    xlabel(valtext{testi});
    ylabel('RT S.D. (s)');
    mysavefig(h, filename, plotdir, fontsize, [6, 1.8]);
    subplot(1,4,4);
    scatter(vals(testi,:), skewrts(testi,:), 28, cols, 'filled');
    xlabel(valtext{testi});
    ylabel('RT skewness');
    mysavefig(h, filename, plotdir, fontsize, [8, 1.8]);
end
%% Analysis changing rates
% plot(diff(Rbar(:,1)), '-','Color', mycol, 'LineWidth',lwd);
% ylabel('Change rate (Hz/ms)');