outdir = '/Volumes/GoogleDrive/My Drive/LDDM/Churchland2008_v2';
addpath('/Users/bs3667/Documents/LDDM/CoreFunctions');
%% multi-alternative choice - dynamic no noise
clevel = 18;
mygray = gray(clevel);
myblue = (repmat([255,255,255],clevel,1)-[linspace(70,255,clevel)',linspace(70,255,clevel)',linspace(40,127,clevel)'])/255;
myred = (repmat([255,255,255],clevel,1)-[linspace(0,25,clevel)',linspace(40,255,clevel)',linspace(40,255,clevel)'])/255;
mygreen = (repmat([255,255,255],clevel,1)-[linspace(70,255,clevel)',linspace(40,127,clevel)',linspace(70,255,clevel)'])/255;

fontsize = 10;
paramspecify_LcDsInhbt;
% alpha beta2, sgm, beta4
samebeta = 1.5;
params = [0, samebeta, 7, samebeta];
% a, b, noise, scale, b2
Tau = [.1, .1, .1];
% Tau = [.1853, .2244, .3231];
% Tau = [.1, .1, .3231];
predur = .6;
dur = 3;
delay = .19;
% Cohr = [3.2, 9, 25.6]/100;
Cohr = [0 32 64 128 256 512 768]/1000;
dt = .001;
sigma = 0;
Buildup = NaN(320,6,2);
Suppress = NaN(320,6,2);
x = [1:130]';
h = figure;
hold on;
for i = 1:2
    items = i*2;
    alpha = params(1) * eye(items);
    beta = params(2) * eye(items);
    if items == 4
        beta = params(4) * eye(items);
    end
    sgm = params(3);
    w = 1*ones(items);
    initialvals = 5*[ones(1,items); items*ones(1,items); zeros(1, items)];
    Rstar = 32;
    Scl = ((1-alpha(1,1))*Rstar + 2*Rstar^2);
    Vprior = Scl*ones(1,items);
    presentt = dt;
    stimdur = dur - presentt;
    stoprule = 1;
    [R0, G0, I0, ~, ~] = LcDsInhbt(Vprior, w, alpha, beta, sigma, Tau, predur, dt, presentt, Inf, Inf, initialvals, stimdur, stoprule);

    presentt = 0;
    stimdur = dur - presentt;
    triggert = presentt;
    stoprule = 1;
    cohi = 0;
    for cohr = Cohr
        cohi = cohi + 1;
        s = cohi + 1;
        if length(Cohr) <= 3
            if i == 1
                mycol = mygray(18-s*4-1,:);
            else
                mycol = myred(s*4+2,:);
            end
        end
        initialvals = [R0(end,:);G0(end,:);I0(end,:)];
        Vinput = Scl/1.5 * [1, ones(1,items-1)];
        [Rd, Gd, Id, ~, ~] = LcDsInhbt(Vinput, w, alpha, beta, sigma, Tau, delay, dt, presentt, presentt+delay, thresh, initialvals, delay, 0);
        initialvals = [Rd(end,:);Gd(end,:);Id(end,:)];
        Vinput = Scl * [(1+cohr), (1-cohr)*ones(1,items-1)];
        [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sigma, Tau, dur-delay, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        Rout = [R0;Rd;R];
        if numel(Cohr) <= 3
            if items == 2
                lp(items) = plot(Rout(:,1),'-','Color',mycol,'LineWidth',1.5);
            else
                lp(items) = plot(Rout(:,1),'-','Color',mycol,'LineWidth',1.5);
            end
            plot(Rout(:,2),'--','Color',lp(items).Color,'LineWidth',1.5);
            if 0 %items == 4
                plot(Rout(:,3),'--','Color',lp(items).Color,'LineWidth',2.5);
                plot(Rout(:,4),'--','Color',lp(items).Color,'LineWidth',2.5);
            end
        end
        endt = min(floor(rt/dt), 480);
        Buildup(1:endt,cohi,i) = R(1:endt,1);
        Suppress(1:endt,cohi,i) = R(1:endt,2);
    end
end
plot([0, 0], [0,70],'k-');
plot(ones(1,2)+(predur)/dt, [0,70],'k-');
ylabel('Activity (H.z.)');
xlabel('Time (ms)');
xlim([-200,2000]);
xticks([0,predur/dt+1,predur/dt+190,predur/dt+320]);
xticklabels({'','','',''});
set(gca,'FontSize',fontsize);
% legend(lp,{'two items','four items'},'Location','North');
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 4.8 3.0];
saveas(h,fullfile(outdir,sprintf('Churchland_4items_%icoh.eps',length(Cohr))),'epsc2');
% build up rate
clear up down;
for i = 1:2
    items = i*2;
    switch i
        case 1
            mycol = 'b';
        case 2
            mycol = 'r';
    end
    for cohi = 1:length(Cohr)
        % min(Buildup(:,cohi,i))
        y_in = Buildup(x,cohi,i);
        p = polyfit(x(~isnan(y_in)),y_in(~isnan(y_in)),1);
        up(cohi,i) = p(1);
        
        y_out = Suppress(x,cohi,i);
        p = polyfit(x(~isnan(y_out)),y_out(~isnan(y_out)),1);
        down(cohi,i) = p(1);
    end
end
%
h = figure; hold on;
plot(Cohr'*100, up(:,1)*1000,'k.-','LineWidth',1,'MarkerSize',14);
plot(Cohr'*100, down(:,1)*1000,'ko--','LineWidth',1,'MarkerSize',5);
plot(Cohr'*100, up(:,2)*1000,'r.-','LineWidth',1,'MarkerSize',14);
plot(Cohr'*100, down(:,2)*1000,'ro--','LineWidth',1,'MarkerSize',5);
xlabel('Motion strength (% coh)');
ylabel('Build up rate (Hz/s)');
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
hold off;
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 3.2];
saveas(h,fullfile(outdir,sprintf('MultiChoice_BuildupRate_%icoh.eps',length(Cohr))),'epsc2');
%% multi-alternative choice, with noise
%samebeta = 2.5;
%params = [0, samebeta, 15, samebeta];
sims = 1000;
gap = 0;
Cohr = [0 32 64 128 256 512 768]/1000; % percent of coherence
name = sprintf('sims%ik_a%1.2f_b%1.2f_noise%2.2f_b%1.2f',sims/1000,params);
file = fullfile(outdir,sprintf('MultiChoice_%s.mat',name));
if exist(file,'file')
    load(file);
else
    % rndseed = rng(2);
    rtmat = [];
    choicemat = [];
    Buildup = NaN(320,sims,7,2); % 131
    Suppress = NaN(320,sims,7,2); % 131
    for i = 1:2
        items = i*2;
        fprintf('items: %i',items);
        alpha = params(1) * eye(items);
        beta = params(2) * eye(items);
        if items == 4
            beta = params(4) * eye(items);
        end
        sgm = params(3);
        w = 1*ones(items);
        initialvals = 8*[ones(1,items); items*ones(1,items); zeros(1, items)];
        Scl = ((1-alpha(1,1))*Rstar + 2*Rstar^2);
        Vprior = Scl*ones(1,items);
        presentt = dt;
        stimdur = dur - presentt;
        triggert = Inf;
        stoprule = 0;
        sigma = 0;
        [R0, G0, I0, rt, choice] = LcDsInhbt(Vprior, w, alpha, beta, sigma, Tau, predur, dt, presentt, triggert, Inf, initialvals, stimdur, stoprule);
        for cohi = 1:length(Cohr)
            fprintf('.');
            cohr = Cohr(cohi);
            dur = 5;
            presentt = 0;
            stimdur = dur - presentt;
            triggert = presentt;
            stoprule = 1;
            for rep = 1:sims
                initialvals = [R0(end,:);G0(end,:);I0(end,:)];
                Vinput = Scl/1.5 * [1, ones(1,items-1)];
                [Rd, Gd, Id, ~, ~] = LcDsInhbt(Vinput, w, alpha, beta, sgm, Tau, delay, dt, presentt, presentt+delay, thresh, initialvals, delay, 0);
                initialvals = [Rd(end,:);Gd(end,:);Id(end,:)];
                Vinput = Scl*[(1+cohr), (1-cohr)*ones(1,items-1)];
                [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
                if isnan(rt)
                    warning('rt is NaN, decision not made');
                else 
                    rt = rt + delay;
                end
                rtmat(rep,cohi,i) = rt;
                choicemat(rep,cohi,i) = choice;
                if rt >= .06+delay
                    endt = min(320-delay/dt, round((rt-delay)/dt));
                    Buildup(1:endt,rep,cohi,i) = R(1:endt,1);
                    Suppress(1:endt,rep,cohi,i) = R(1:endt,2);
                end
            end
        end
        fprintf('\n');
    end
    save(fullfile(outdir,sprintf('MultiChoice_%s.mat',name)),...
        'rtmat','choicemat','Buildup','Suppress','params');
end
choicemat = choicemat == 1;
% ACC & RT, build-up rate
fontsize = 14;
head = .012;
h = figure;
subplot(2,1,1); hold on;
meanchoice = mean(choicemat, 1, 'omitnan');
plot([head, Cohr(2:end)], meanchoice(:,:,1),'k.-','LineWidth',1,'MarkerSize',14);
plot([head, Cohr(2:end)], meanchoice(:,:,2),'r.-','LineWidth',1,'MarkerSize',14);
ylim([0,1]);
xlim([head,1]);
yticks([0:.25:1]);
xticks([head, .1, 1]);
xticklabels({'0','10','100'});
ylabel('Probability correct');
set(gca,'XScale','log');
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
hold off;

subplot(2,1,2); hold on;
meanrt = mean(rtmat,1,'omitnan')/dt;
plot([head, Cohr(2:end)], meanrt(:,:,1),'k.-','LineWidth',1,'MarkerSize',14);
plot([head, Cohr(2:end)], meanrt(:,:,2),'r.-','LineWidth',1,'MarkerSize',14);
xlim([head,1]);
yticks([200:200:1800]);
xticks([head, .1, 1]);
xticklabels({'0','10','100'});
xlabel('Motion strength (% coh)');
ylabel('Reaction time (ms)');
set(gca,'XScale','log');
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
hold off;
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 3.2];
saveas(h,fullfile(outdir,sprintf('MultiChoice_RT&ACC_%s.eps',name)),'epsc2');

% build up rate
x = [1:130]';
h = figure; hold on;
for i = 1:2
    items = i*2;
    switch i
        case 1
            mycol = 'b';
        case 2
            mycol = 'r';
    end
    for cohi = 1:length(Cohr)
        in = choicemat(:,cohi,i) == 1;
        longt = 1;%rtmat(:,cohi,i) > .45;
        block_in = Buildup(1:130,in&longt,cohi,i);
        y_in = mean(block_in,2,'omitnan');
        p = polyfit(x,y_in,1);
        up(cohi,i) = p(1);
        mdl = fitlm(x,y_in,'linear')
        plot(x,y_in,'.');
        %plot(x,p(1)*x+p(2),'-','Color',mycol);
        
        out = choicemat(:,cohi,i) == 1;
        block_out = Suppress(1:130,out&longt,cohi,i);
        y_out = mean(block_out,2,'omitnan');
        p = polyfit(x,y_out,1);
        down(cohi,i) = p(1);
        mdl = fitlm(x,y_out,'linear')
        plot(x,y_out,'.');
        %plot(x,p(1)*x+p(2),'--','Color',mycol);
    end
end
%
h = figure; hold on;
plot(Cohr'*100, up(:,1)*1000,'k.-','LineWidth',1,'MarkerSize',14);
plot(Cohr'*100, down(:,1)*1000,'ko--','LineWidth',1,'MarkerSize',5);
plot(Cohr'*100, up(:,2)*1000,'r.-','LineWidth',1,'MarkerSize',14);
plot(Cohr'*100, down(:,2)*1000,'ro--','LineWidth',1,'MarkerSize',5);
xlabel('Motion strength (% coh)');
ylabel('Build up rate (Hz/s)');
set(gca,'FontSize',fontsize-2);
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
hold off;
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 3.2];
saveas(h,fullfile(outdir,sprintf('MultiChoice_BuildupRate_%s.eps',name)),'epsc2');
%% averaged trace, Tin and Tout
Cohr = [0, 9, 25.6, 51.2]/100;
absCohr = [1,2,3,4];
rtmat = [];
choicemat = [];
R1mat = []; %NaN(round((predur+dur)/dt),1000,5,2);
R2mat = []; %NaN(round((predur+dur)/dt),1000,5,2);
name = sprintf('sims%ik_a%1.2f_b%1.2f_noise%2.2f_b%1.2f',sims/1000,params);
file = fullfile(outdir,sprintf('MultiChoiceTrace_%s.mat',name));
if exist(file,'file')
    load(file);
else
    % rndseed = rng(2);
    for i = 1:2
        items = i*2;
        fprintf('items: %i\n',items);
        alpha = params(1) * eye(items);
        beta = params(2) * eye(items);
        if items == 4
            beta = params(4) * eye(items);
        end
        sgm = params(3);
        w = 1*ones(items);
        initialvals = 8*[ones(1,items); items*ones(1,items); zeros(1, items)];
        Scl = ((1-alpha(1,1))*Rstar + 2*Rstar^2);
        Vprior = Scl*ones(1,items);
        %Vprior = params(4)*512*ones(1,items);
        predur = .6;
        presentt = dt;
        stimdur = predur - presentt;
        triggert = Inf;
        stoprule = 0;
        sigma = 0;
        [R0, G0, I0, rt, choice] = LcDsInhbt(Vprior, w, alpha, beta, sigma, Tau, predur+delay, dt, presentt, triggert, Inf, initialvals, stimdur+delay, stoprule);
        for cohi = 1:length(Cohr)
            cohr = Cohr(cohi);
            fprintf('coherence %2.1f\n', cohr*100);
            % Vinput = params(4)*[256*(1+cohr), 256*(1-cohr)*ones(1,items-1)];
            Vinput = Scl*[(1+cohr), (1-cohr)*ones(1,items-1)];
            initialvals = [R0(end,:);G0(end,:);I0(end,:)];
            dur = 5;
            presentt = 0;
            stimdur = dur - presentt;
            triggert = presentt;
            stoprule = 0;
            for rep = 1:sims
                [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sgm, Tau, dur-delay, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
                if isnan(rt)
                    warning('rt is NaN, decision not made');
                end
                rtmat(rep,cohi,i) = rt+delay;
                choicemat(rep,cohi,i) = choice;
                % Rvec(rep,)) = [R0(:,1); R(:,1)]';
                R1mat(1:round((predur+dur)/dt),rep,cohi,i) = [R0(:,1); R(:,1)];
                R2mat(1:round((predur+dur)/dt),rep,cohi,i) = [R0(:,2); R(:,2)];
            end
        end
    end
    save(fullfile(outdir,sprintf('MultiChoiceTrace_%s.mat',name)),...
        'rtmat','choicemat','R1mat','R2mat','params');
end
%
clevel = 18;
mygray = gray(clevel);
myblue = (repmat([255,255,255],clevel,1)-[linspace(70,255,clevel)',linspace(70,255,clevel)',linspace(40,127,clevel)'])/255;
myred = (repmat([255,255,255],clevel,1)-[linspace(0,25,clevel)',linspace(40,255,clevel)',linspace(40,255,clevel)'])/255;
mygreen = (repmat([255,255,255],clevel,1)-[linspace(70,255,clevel)',linspace(40,127,clevel)',linspace(70,255,clevel)'])/255;

fontsize = 14;
h = figure; hold on;
R1sac = [];
R2sac = [];
R1stim = R1mat;
R1tmp = R1mat;
R2stim = R2mat;
R2tmp = R2mat;
for i = 1:2
    for cohi = 1:length(Cohr)
        cohr = Cohr(cohi);
        s = absCohr(cohi);
        if i == 1
            mycol = mygray(18-s*4-1,:);
        else
            mycol = myred(s*4+2,:);
        end
        for rep = 1:sims
            rtloc = round(rtmat(rep,cohi,i)/dt);
            R1tmp(1:(predur/dt+200),rep,cohi,i) = NaN;
            R2tmp(1:(predur/dt+200),rep,cohi,i) = NaN;
            R1sac(:,rep,cohi,i) = R1tmp(round(predur/dt)+[(rtloc-250):(rtloc+50)],rep,cohi,i);
            R2sac(:,rep,cohi,i) = R2tmp(round(predur/dt)+[(rtloc-250):(rtloc+50)],rep,cohi,i);
            R1stim((round((predur+rtmat(rep,cohi,i))/dt)-60):end,rep,cohi,i) = NaN;
            R2stim((round((predur+rtmat(rep,cohi,i))/dt)-60):end,rep,cohi,i) = NaN;
        end
        %         if cohr <= 0
        %             out = choicemat(:,cohi,i) ~= 1;
        %             % longt = rtmat(:,cohi,i) > .45;
        %             Rout = Rstim([1:round((predur+.33)/dt)], out, cohi, i);
        %             Rout_sac = Rsac(:, out, cohi, i);
        %             plot([mean(Rout,2,'omitnan'); NaN(200,1); mean(Rout_sac,2,'omitnan')],'--','Color',mycol,'LineWidth',1.5);
        %         end
        if cohr >= 0
            in = choicemat(:,cohi,i) == 1;
            longt = 1;%rtmat(:,cohi,i) >= (.45-delay);
            Rin = R1stim([1:round((predur+.33)/dt)], in&longt, cohi, i);
            Rin_sac = R1sac(:, in&longt, cohi, i);
            lp(cohi+(i-1)*length(Cohr)) = plot([mean(Rin,2,'omitnan'); NaN(100,1); mean(Rin_sac,2,'omitnan')],'-','Color',mycol,'LineWidth',1.5);
            
            Rout = R2stim([1:round((predur+.33)/dt)], in&longt, cohi, i);
            Rout_sac = R2sac(:, in&longt, cohi, i);
            plot([mean(Rout,2,'omitnan'); NaN(100,1); mean(Rout_sac,2,'omitnan')],'--','Color',mycol,'LineWidth',1.5);
        end
    end
end
plot(ones(1,2)*0+predur/dt, [0,80],'k--');
plot(ones(1,2)*190+predur/dt, [0,80],'k--');
plot(ones(1,2)*320+predur/dt, [0,80],'k--');
patch([190+predur/dt ones(1,2)*320+predur/dt 190+predur/dt], [0, 0, 80, 80], [0.8 0.8 0.8], 'FaceAlpha',.3);
xlim([-50, 1400]);
xticks([0,.2,.4,.6,.8,1.08,1.28]*1000);
xticklabels({'cues','-.4','-.2','motion','.2','-.2','action'});
ylim([0,71]);
xlabel('Time (secs)');
ylabel('Firing rates (Hz)');
legend(lp([4,8]),{'Two choice','Four Choice'},'Location','NorthWest');
legend('boxoff');
set(gca,'TickDir','out');
H = gca;
H.LineWidth = 1;
set(gca,'FontSize',fontsize);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 4.8 4];
saveas(h,fullfile(outdir,sprintf('MultiChoice_Trace_%s.eps',name)),'epsc2');