%% panel a, dynamic of neural firing rates
outdir = './rslts/LDDM_IIA_Dynamic';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
fontsize = 14;
mksz = 5;
lwd = 2;

eqlbvec = 42; %[5:6:70]; % 11 e
alen = 9;
amat = nan([alen, numel(eqlbvec)]); % 9a * 11e
for ei = 1:numel(eqlbvec)
    maxa = 3*eqlbvec(ei) + 1;
    afit = (3*4.7*eqlbvec(ei) + 4.7)/5.7; % based on experience, when alpha*R* and input scale in a ratio of 4.7, the model works the best
    amat(:,ei) = [flip(afit - logspace(0,log10(afit-1),ceil(alen/2))) afit + logspace(0,log10(maxa-afit-1),floor(alen/2))];
end
eqlbmat = repmat(eqlbvec,alen,1);% 9a * 11e

dt = .001; % second
tauR = .1; % second
tauG = .1; % second
tauI = .1; % second
Tau = [tauR, tauG, tauI];
b = .2*eye(3);
w = ones(3);
prepresentt = dt;
predur = 1;%5.6;
prestimdur = predur - prepresentt;
presentt = dt;
dur = 4;
triggert = presentt;
stimdur = dur-presentt;
thresh = 70;
stoprule = 1;
c0 = [1,1,1];
%c = [3.2 12.8, 25.6, 38.4 51.2]'/100;
c3 = [0:.2:1]';
c1 = 1+.256;
c2 = 1;
cp = [ones(size(c3))*c1, c2*ones(size(c3)), c3];
initialset = [1,1,1;3,3,3;0,0,0];
mygray = flip(gray(length(c3) + 1));
mycol2 = colormap(winter(length(c3)));
aspect = [3,3]*ceil(sqrt(alen));
task = 'PrepSgnlABS';

%% in space of R1-R2
sgm = 0;
simname = sprintf('LDDM_DynamicR1R2_%s_c1%1.2f_c2%1.2f_%ic3_%ialevel_b%1.2f_w%1.1f_%ieqlbvals_init%1.2f_sgm%2.1f',...
    task,c1,c2,length(c3),alen,b(1,1),w(1,1),numel(eqlbvec),initialset(1,1),sgm);
lim = 0;
for ei = 1:numel(eqlbvec) % the same equilibirum level and starting points but different alpha level
    h = figure; 
    filename = sprintf('%s_eqlb%2.1f',simname, eqlbvec(ei));
    for ai = 1:alen
        subplot(ceil(sqrt(alen)),ceil(sqrt(alen)),ai); hold on;
        scale = 3*eqlbvec(ei).^2 + (1-amat(ai,ei)).*eqlbvec(ei); % 9a * 11e
        a = amat(ai,ei)*eye(3);
        V0 = c0*scale;
        lgd3 = [];
        for vi = 1:numel(c3)
            initialvals = initialset*eqlbvec(ei);
            [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
            Vinput = cp(vi,:)*scale;
            % initialvals = [R0(end,:);G0(end,:);I0(end,:)];
            initialvals = initialset*eqlbvec(ei);
            [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            %R0(end-100:end,:) = NaN;
            R = [R0; R];
            lgd3(vi) = plot(R(:,1), R(:,2), '-', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
            plot([0,R(end,1)],[0,R(end,2)],'.-', 'Color', mycol2(vi,:), 'LineWidth',lwd/2,'MarkerSize',mksz*2);
            lim = max([lim, max(R(:,1:2))]);
        end
        if 1
            plot([0, 2*thresh],[thresh,thresh], 'k--');
            plot([thresh,thresh],[0, 2*thresh], 'k--');
        end
        plot([0,thresh],[0,thresh], 'k--');
        ylim([0,lim*1.2]);
        ylabel('R2 Activity (a.u.)');
        xlim([0,lim*1.2]);
        xlabel('R1 Activity (a.u.)');
        if ai == 1
            title({sprintf('Eqlbrm %1.1f, slfexctlevel %i',eqlbvec(ei),ai),sprintf('alpha %3.1f, scale %2.0f',a(1,1),scale)},'FontSize',fontsize-7,'FontWeight','normal');
            lgd = legend(lgd3,cellstr(num2str(c3)),...
                'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
                'FontAngle','italic','NumColumns',1,'Box','off');
            title(lgd, 'V_3');
        else
            title({sprintf('slfexctlevel %i',ai), sprintf('alpha %3.1f, scale %2.0f',a(1,1),scale)},'FontSize',fontsize-7,'FontWeight','normal');
            
        end
        savefigs(h, filename, outdir, fontsize, aspect);
    end
end

%% in the space of time - R
simname = sprintf('LDDM_Dynamic_%s_c1%1.2f_c2%1.2f_%ic3_%ialevel_b%1.2f_w%1.1f_%ieqlbvals_init%1.2f_sgm%2.1f',...
    task,c1,c2,length(c3),alen,b(1,1),w(1,1),numel(eqlbvec),initialset(1,1),sgm);
for ei = 1:numel(eqlbvec) % the same equilibirum level and starting points but different alpha level
    h = figure; 
    filename = sprintf('%s_eqlb%2.1f',simname, eqlbvec(ei));
    for ai = 1:alen
        subplot(ceil(sqrt(alen)),ceil(sqrt(alen)),ai); hold on;
        scale = 3*eqlbvec(ei).^2 + (1-amat(ai,ei)).*eqlbvec(ei); % 9a * 11e
        a = amat(ai,ei)*eye(3);
        V0 = c0*scale;
        for vi = 1:numel(c3)
            initialvals = initialset*eqlbvec(ei);
            [R0, G0, I0, ~, ~] = LcDsInhbt(V0, w, a, b, sgm, Tau, predur, dt, prepresentt, Inf, Inf, initialvals, prestimdur, stoprule);
            Vinput = cp(vi,:)*scale;
            % initialvals = [R0(end,:);G0(end,:);I0(end,:)];
            initialvals = initialset*eqlbvec(ei);
            [R, ~, ~, ~, ~] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            %R0(end-100:end,:) = NaN;
            R = [R0; R];
            lgd3(vi) = plot(R(:,3), 'k-.', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
            lgd2(vi) = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
            lgd1(vi) = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
        end
        title(sprintf('alevel %i, Eqlbrm %1.1f',ai,eqlbvec(ei)),'FontSize',fontsize-7,'FontWeight','normal');
        plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
        % text(1200,thresh*1.1,'threshold');
        
        plot([1, prepresentt/dt prepresentt/dt],...
            [0,0,2]-7,'k-','LineWidth',lwd/1.5);
        plot([prepresentt/dt, (prepresentt+prestimdur)/dt],...
            [2,2]-7,'k--','LineWidth',lwd/1.5); % input, target, pre-motion
        plot([(predur)/dt, (predur+presentt+dur)/dt],...
            [2,2]-7,'k-','LineWidth',lwd/1.5); % inputs, stimuli, motion
        %plot([1, (prepresentt)/dt, (prepresentt)/dt, (predur)/dt, (predur)/dt, (predur+presentt+dur)/dt],[2,2,2,2,0,0]-11.5,'k-','LineWidth',lwd/1.5); % fixation point (go signal)
        plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-11.5,'k-','LineWidth',lwd/1.5); % disinhibition
        ylim([-16.5,70]);
        yticks([0]);
        yticklabels({'0'});
        ylabel('Activity (a.u.)');
        xticks([]);
        drawaxis(gca, 'x', 0, 'movelabel', 1);
        xlim([1, (prepresentt+predur+presentt+dur)/dt]);
        if ai == 1
            title({sprintf('Eqlbrm %1.1f, slfexctlevel %i',eqlbvec(ei),ai),sprintf('alpha %3.1f, scale %2.0f',a(1,1),scale)},'FontSize',fontsize-7,'FontWeight','normal');
            lgd = legend(lgd3,cellstr(num2str(c3)),...
                'Location','best','FontSize',fontsize-5, 'FontName','Times New Roman',...
                'FontAngle','italic','NumColumns',1,'Box','off');
            title(lgd, 'V_3');
        else
            title({sprintf('slfexctlevel %i',ai), sprintf('alpha %3.1f, scale %2.0f',a(1,1),scale)},'FontSize',fontsize-7,'FontWeight','normal');
        end
        savefigs(h, filename, outdir, fontsize, aspect);
    end
end
