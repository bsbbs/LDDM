% sgmi = 1;
% name = 'LcD';
% name = 'AsymW';
% name = 'XJ';
%% simulation
if strcmp(name, 'LcD')
    paramspecify_LcDsInhbt;
    thresh = 50;
    a = eye(2)*5;
    b = eye(2)*1.4;
    presentt = .2;
    triggert = .9;
elseif strcmp(name, 'AsymW')
    paramspecify_LcDsInhbt;
    thresh = 50;
    a = eye(2)*5;
    w = [.1, 1;
        1, .1];
    presentt = .2;
    triggert = .9;
elseif strcmp(name, 'XJ')
    paramspecify_WongWang;
    thresh = 15;
    presentt = .2;
    triggert = presentt;
end
Vinput = [320, 192];
dt = .001;
dur = 1.5;

stoprule = 1;
stimdur = dur - presentt;%.7;
switch sgmi
    case 1
        sgm = 0;
    case 2
        if strcmp(name, 'XJ')
            sgm = .008;
        else
            sgm = .5;
        end
end
h = figure;
hold on;
mycl = jet(length(GABA)*10);
for GABAi = 1:length(GABA)
    if strcmp(name,'LcD')
        [R, G, I, rt, choice] = LcDsInhbt_GABAblock(Vinput, GABA(GABAi), w,...
            a, b, sgm, Tau, dur, dt, presentt, triggert, thresh*15, initialvals, stimdur, stoprule);
        R(R>thresh) = NaN;
        slide_wind = dt;
    elseif strcmp(name,'AsymW')
        [R, G, I, rt, choice] = AsymW_GABA(Vinput, GABA(GABAi), w, a, sgm, Tau, dur,...
            dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    elseif strcmp(name,'XJ')
        [nu_wind, s_wind, rt, choice, H, S] = wong06_Gaba(Vinput,GABA(GABAi), miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals_dflt, stoprule);
        R = nu_wind;
        slide_wind = dt;
    end
    plot(R(:,1), 'LineWidth', lwd, 'Color',mycl(5+(GABAi-1)*10,:)); % , 'Color', mycl(GABAi,:)
    plot(R(:,2), '--', 'LineWidth', lwd, 'Color',mycl(5+(GABAi-1)*10,:));
end
plot(triggert/slide_wind:dur/slide_wind,ones(1,length(triggert/slide_wind:dur/slide_wind))*thresh, '--k');
if strcmp(name,'LcD')
    text(length(R) - .4/slide_wind,thresh*1.05,'threshold');
else
    text(dur/slide_wind - .3/slide_wind,thresh*1.1,'threshold');
end
xlabel('Time (a.u.)');
if strcmp(name,'XJ') || strcmp(name, 'AsymW')
    xticks([presentt/dt]);
    xticklabels({'on'});
else
    xticks([presentt/dt, triggert/dt]);
    xticklabels({'stimuli','action'}); % ,'3.5','4.0','4.5'
end
ylabel('Firing Rates (Hz)');
% if strcmp(name,'LcD')
    autoylim = [0,thresh*1.1];
% elseif 
%     autoylim = [0,20];
% end
ylim(autoylim*1.1);
% title(modelname);
lgd = legend('R1 Placebo', 'R2 Placebo', 'R1 Agonist', 'R2 Agonist', 'Location','Best','NumColumns',2, 'FontSize', fontsize-5);
if ~exist('./graphics/timeCourse','dir')
    mkdir('./graphics/timeCourse');
end
filename = sprintf('timeCourse/%s_GABA%1.1f_%1.1f_FD_sgm%2.3f.eps',name,GABA, sgm);
savefig(h, filename, outdir, fontsize, aspect3);