%% Fig 8 - Persistent activity with and without local disinhibition

%% phase plane analysis and vector field. i = 1: b = 0; i = 2: 0<b<w; i = 3: b>w
lgdvec = [];
for bi = 1:2
    switch bi
        case 1
            filename = 'Fig8b';
        case 2
            filename = 'Fig8e';
    end
    Taup = [1,1,1];
    alpha = eye(2)*a0;
    a = alpha(1,1);
    w = [1,1;1,1];
    if bi == 1
        b = 0;
    elseif bi == 2
        b = .4;
    elseif bi == 3
        b = 1.2;
    end
    Dt = 1;
    if w(1,1) == w(1,2)
        if b > 0 && b < 1
            axismax = 25;
        elseif b == 0
            axismax = 20;
        elseif b >= 1
            axismax = 30;
        end
    else
        if w(1,1) > w(1,2)
            axismax = (alpha(1,1) - 1)/w(1,2)*1.1;
        elseif w(1,1) < w(1,2)
            axismax = (alpha(1,1) - 1)/w(1,1)*1.1;
        end
    end
    [R10,R20] = meshgrid(linspace(0,axismax,15),linspace(0,axismax,15));
    R10(1,1) = R10(1,2)/3;
    R20(1,1) = R20(2,1)/3;
    G1 = R10*w(1,1) + R20*w(1,2) - b*R10;
    G2 = R10*w(2,1) + R20*w(2,2) - b*R20;
    G1(G1<0) = 0;
    G2(G2<0) = 0;
    dR1 = (-R10 + (R10.*alpha(1,1))./(1+G1)) * Dt/Tau(1);
    dR2 = (-R20 + (R20.*alpha(2,2))./(1+G2)) * Dt/Tau(2);
    norm = sqrt(dR1.^2 + dR2.^2);
    rate = (norm.^.6)./(norm);
    dR1 = dR1.*rate;
    dR2 = dR2.*rate;
    h = figure;
    hold on;
    if w(2,2) - b > 0 && w(1,2) == w(1,1)
        R2Line = [0,(alpha(2,2) - 1)/(w(2,2) - b)];
    elseif w(2,2) - b <= 0 && w(1,2) == w(1,1)
        R2Line = [0,axismax];
    elseif b == 0 && w(1,2) ~= w(1,1)
        R2Line = [0,(alpha(1,1) - 1)/w(1,2)];
    end
    R1Line = -w(1,2)/(w(1,1)-b).*R2Line+(alpha(1,1) - 1)/(w(1,1)-b);
    lg1 = plot(R2Line,R1Line,'-','Color',colorpalette{4},'LineWidth',2);% dR1/dt = 0
    if w(2,2) - b > 0 && w(1,2) == w(1,1)
        R1Line = [0,(alpha(1,1) - 1)/(w(1,1) - b)];
    elseif w(2,2) - b <= 0 && w(1,2) == w(1,1)
        R1Line = [0,axismax];
    elseif b == 0 && w(1,2) ~= w(1,1)
        R1Line = [0,(alpha(2,2) - 1)/w(2,1)];
    end
    R2Line = -w(2,1)/(w(2,2)-b).*R2Line+(alpha(2,2) - 1)/(w(2,2)-b);
    lg2 = plot(R2Line,R1Line,'--','Color',colorpalette{1},'LineWidth',2);% dR2/dt = 0
    syms R1 R2;
    vars = [R1 R2];
    eqns = [(w(1,1)-b)*R1 + w(1,2)*R2 == alpha(1,1) - 1, ... % dR1/dt = 0
        w(2,1)*R1 + (w(2,2)-b)*R1 == alpha(2,2) - 1];% dR2/dt = 0
    [AnswR1,AnswR2] = solve(eqns, vars);
    AnswI1 = b*AnswR1;
    AnswI2 = b*AnswR2;
    AnswG1 = w(1,1)*AnswR1 + w(1,2)*AnswR2 - AnswI1;
    AnswG2 = w(2,1)*AnswR1 + w(2,2)*AnswR2 - AnswI2;
    PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
    Npoints = sum(PositiveRealAnsw);
    for i = find(PositiveRealAnsw)'
        R1star = double(AnswR1(i));
        R2star = double(AnswR2(i));
        G1star = double(AnswG1(i));
        G2star = double(AnswG2(i));
        JMat = [-1 + a/(1+G1star), -(a*R1star)/(1+G1star)^2, 0, 0, 0, 0
            w(1,1), -1, -1, w(1,2), 0, 0
            b, 0, -1, 0, 0, 0
            0, 0, 0, -1+a/(1+G2star),   -(a*R2star)/(1+G2star)^2, 0
            w(2,1), 0, 0, w(2,2), -1, -1
            0, 0, 0, b, 0, -1];
        A = eig(JMat);
        if all(real(A) < -1e-10)
            lgdvec(end +1) = plot(R2star, R1star, '.', 'Color','#06d6a0', 'MarkerSize',mksz/1.5);
            plot(R1star, R2star, 'o', 'Color','k', 'MarkerSize',mksz/4);
        elseif sum(real(A) > 1e-10)
            lgdvec(end +1) = plot(R2star,R1star,  '.', 'Color', '#ffd166', 'MarkerSize',mksz/1.5);
            plot(R1star, R2star, 'o', 'Color','k', 'MarkerSize',mksz/4);
        end
    end
    if w(1,1) == w(1,2)
        if b > 0 && b < 1
            xlim([-.6,axismax]);ylim([-.6,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1)+w(1,2)-b), (alpha(1,1)-1)/(w(1,1) - b)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2)+w(2,1)-b), (alpha(2,2)-1)/(w(2,2) - b)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$','$\frac{\alpha-1-G_0}{\omega-\beta}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$','$\frac{\alpha-1-G_0}{\omega-\beta}$'});
        elseif b == 0
            xlim([-.4,axismax]);ylim([-.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) - b)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) - b)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{\omega}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{\omega}$'});
        elseif b >= 1
            xlim([-1.4,axismax]);ylim([-1.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2) - b)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1) - b)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$'});
        end
    elseif w(1,1) ~= w(1,2) && b == 0
        if w(1,1) > w(1,2) % convergent
            xlim([-.8,axismax]);ylim([-.8,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2)),(alpha(1,1)-1)/w(1,2)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1)),(alpha(2,2)-1)/w(2,1)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{v}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{v}$'});
        elseif w(1,1) < w(1,2) % divergent
            xlim([-.4,axismax]);ylim([-.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2)),(alpha(1,1)-1)/w(1,1)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1)),(alpha(2,2)-1)/w(2,2)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{w}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{w}$'});
        end
    end
    set(gca,'TickLabelInterpreter', 'latex');
    xlabel('R_2 activity (a.u.)', 'FontAngle','italic');
    ylabel('R_1 activity (a.u.)', 'FontAngle','italic');
    lg3 = quiver(R10,R20,dR1,dR2,1.3,'k','LineWidth',1);
    ratio = aspect5;
    
    if bi == 1
        legend([lg1,lg2,lg3], {'\color[rgb]{.0667,.5412,.6980}dR_1/dt = 0',...
            '\color[rgb]{.9373,.2784,.4353}dR_2/dt = 0','\color{black}Change rate'},...
            'FontName','Times New Roman', ...
            'FontAngle','Italic','FontSize',fontsize-4, 'Location','northeastoutside','Box','off');
    elseif bi == 2
        legend([lg1,lg2,lg3,lgdvec(1)], {'\color[rgb]{.0667,.5412,.6980}dR_1/dt = 0',...
            '\color[rgb]{.9373,.2784,.4353}dR_2/dt = 0','\color{black}Change rate',...
            'Unstable point'}, 'FontName','Times New Roman', ...
            'FontAngle','Italic','FontSize',fontsize-4, 'Location','northeastoutside','Box','off');
    end
    savefigs(h, filename, plotdir, 11, ratio);
end

%% example dynamic of R1 and R2. i = 1: b = 0; i = 2: 0<b<w; i = 3: b>w
mygrayred = mygray;
mygrayred(:,1) = 1;
mygrayblue = mygray;
mygrayblue(:,3) = 1;
VR = [];
WM = [];
for bi = 1:2
    switch bi
        case 1
            filename = 'Fig8a';
        case 2
            filename = 'Fig8d';
    end
    
    a = a0*eye(2);
    w = [1,1;1,1];
    if bi == 1
        b = 0*eye(2);
    elseif bi == 2
        b = .4*eye(2);
    elseif bi == 3
        b = 1.2*eye(2);
    end
    initialvals = zeros(3,2);
    predur = 0;
    presentt = .1;
    dur = 3.2;
    stimdur = .9-presentt;
    sgm = 0;
    triggert = .9;
    thresh = Inf;
    stoprule = 0;
    Vprior = [1, 1]*scale0;
    h = figure; hold on;
    for vi = 2:5
        Vinput = [1+c(vi), 1-c(vi)]*scale0;
        [choice,rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
        lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd/2);
        VR(bi,vi,:) = R((predur+presentt+stimdur)/dt-1,:);
        WM(bi,vi,:) = R((predur+dur)/dt-1,:);
    end
    ylim([-.4,70]);
    yticks([0]);
    xlim([1, dur/dt]);
    xticks([presentt/dt, (presentt+stimdur)/dt]);
    xticklabels({'on','off'});
    ylabel('Activity (a.u.)');
    xlabel('Time (a.u.)');
    if bi == 1
    legend([lgd1, lgd2],{'R_1', 'R_2'},...
        'Location','North','FontSize',fontsize-5, 'FontName','Times New Roman', ...
        'FontAngle','italic','NumColumns',1,'Box','off');
    end
    mysavefig(h, filename, plotdir, fontsize - 2, aspect10);
end
%% coding two items
mygrayred = mygray;
mygrayred(:,1) = 1;
mygrayblue = mygray;
mygrayblue(:,3) = 1;
cp = [-flip(c(2:5)); c(2:5)];
mygray10 = flip(gray(length(cp) + 2));
VR = [];
WM = [];
h = figure;
filename = 'Fig8_CodingTwoItems';
for bi = 1:3
    a = a0*eye(2);
    w = ones(2);
    if bi == 1
        b = 0*eye(2);
    elseif bi == 2
        b = .4*eye(2);
    elseif bi == 3
        b = 1.2*eye(2);
    end
    initialvals = zeros(3,2);
    predur = 0;
    presentt = .1;
    dur = 3.2;
    stimdur = .9-presentt;
    sgm = 0;
    triggert = .9;
    thresh = Inf;
    stoprule = 0;
    Vprior = [1, 1]*scale0;
    for vi = 1:numel(cp)
        Vinput = [1+cp(vi), 1-cp(vi)]*scale0;
        [choice,rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        VR(bi,vi,:) = R((predur+presentt+stimdur)/dt-1,:);
        WM(bi,vi,:) = R((predur+dur)/dt-1,:);
    end
    
    si = 1;
    subplot(3,2,2*(bi-1)+si); hold on;
    for vi = 1:numel(cp)/2
        plot([0],[VR(bi,vi,1)],'.', 'Color', mygrayred(8-vi,:),'LineWidth',2,'MarkerSize',18);
        plot([VR(bi,vi,2)],[0],'.', 'Color', mygrayred(8-vi,:),'LineWidth',2,'MarkerSize',18);
        plot([0,VR(bi,vi,2),VR(bi,vi,2)],[VR(bi,vi,1),VR(bi,vi,1),0],'-','LineWidth',.25,'Color',mygray(2,:));
        plot([0,VR(bi,vi,2)],[0,VR(bi,vi,1)],'-','LineWidth',.25,'Color',mygray(3,:));
    end
    for vi = 5:numel(cp)
        plot([0],[VR(bi,vi,1)],'.', 'Color', mygrayred(vi-2,:),'LineWidth',2,'MarkerSize',18);
        plot([VR(bi,vi,2)],[0],'.', 'Color', mygrayred(vi-2,:),'LineWidth',2,'MarkerSize',18);
        plot([0,VR(bi,vi,2),VR(bi,vi,2)],[VR(bi,vi,1),VR(bi,vi,1),0],'-','LineWidth',.25,'Color',mygray(2,:));
        plot([0,VR(bi,vi,2)],[0,VR(bi,vi,1)],'-','LineWidth',.25,'Color',mygray(3,:));
    end
    % plot([0,a0-1],[a0-1,0],'-','LineWidth',.25,'Color',mygray(2,:));
    xlim([0,25]);
    ylim([0,25]);
    yticks([0,20]);
    xticks([0,20]);
    drawaxis(gca, 'x', 0, 'movelabel', 1);
    drawaxis(gca, 'y', 0, 'movelabel', 1);
    xlabel('R_2');
    ylabel('R_1');
    
    si = 2;
    subplot(3,2,2*(bi-1)+si); hold on;
    for vi = 1:numel(cp)/2
        plot([0],[WM(bi,vi,1)],'.', 'Color', mygrayblue(8-vi,:), 'LineWidth',2,'MarkerSize',12); %
        plot([WM(bi,vi,2)],[0],'.', 'Color', mygrayblue(8-vi,:),'LineWidth',2,'MarkerSize',12); % 
        plot([0,WM(bi,vi,2),WM(bi,vi,2)],[WM(bi,vi,1),WM(bi,vi,1),0],'-','LineWidth',.25,'Color',mygray(2,:));
        plot([0,WM(bi,vi,2)],[0,WM(bi,vi,1)],'-','LineWidth',.25,'Color',mygray(3,:));
    end
    for vi = 5:numel(cp)
        plot([0],[WM(bi,vi,1)],'.', 'Color', mygrayblue(vi-2,:), 'LineWidth',2,'MarkerSize',12); %
        plot([WM(bi,vi,2)],[0],'.', 'Color', mygrayblue(vi-2,:),'LineWidth',2,'MarkerSize',12); % 
        plot([0,WM(bi,vi,2),WM(bi,vi,2)],[WM(bi,vi,1),WM(bi,vi,1),0],'-','LineWidth',.25,'Color',mygray(2,:));
        plot([0,WM(bi,vi,2)],[0,WM(bi,vi,1)],'-','LineWidth',.25,'Color',mygray(3,:));
    end
    
    if bi == 1
        plot([0,a0-1],[a0-1,0],'-','LineWidth',.25,'Color',mygray(2,:));
    end
    xlim([0,25]);
    ylim([0,25]);
    yticks([0,20]);
    xticks([0,20]);
    drawaxis(gca, 'x', 0, 'movelabel', 1);
    drawaxis(gca, 'y', 0, 'movelabel', 1);
    xlabel('R_2');
    ylabel('R_1');
end
savefigs(h, filename, plotdir, fontsize - 2, [4.5,6.5]);
%% Coding multiple items
Vprior = [1,1,1,1,1]*scale0/5;
VR = [];
WM = [];
for bi = 1:2
    a = a0*eye(5)*2.5;
    w = ones(5);
    if bi == 1
        b = 0*eye(5);
    elseif bi == 2
        b = .1*eye(5);
    elseif bi == 3
        b = 1.2*eye(5);
    end
    initialvals = zeros(3,5);
    predur = 0;
    presentt = .1;
    dur = 11;
    stimdur = .9-presentt;
    sgm = 0;
    triggert = .9;
    thresh = Inf;
    stoprule = 0;
    for vi = 2:5
        Vinput = [1+c(vi), 1-c(vi), 1-c(vi), 1-c(vi), 1-c(vi)]*scale0/5;
        [choice,rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
%         lgd1 = plot(R(:,1), 'k-', 'Color', mygrayred(vi+1,:), 'LineWidth',lwd/2);
%         lgd2 = plot(R(:,2), 'k--', 'Color', mygrayblue(vi+1,:), 'LineWidth',lwd/2);
        VR(bi,vi,:) = R((predur+presentt+stimdur)/dt-1,:);
        WM(bi,vi,:) = R((predur+dur)/dt-1,:);
    end
    
    % plot
    maxy = 45; %max([max(max(WM(bi,2:5,:))),max(max(VR(bi,2:5,:)))]);
    si = 1;
    h = figure;
    filename = sprintf('Fig8cf_CodingFiveItems_VR_b%i_s%i',bi,si);
    s = spider_plot_class(squeeze(VR(bi,2:5,:)), 'LegendLabels',cellstr(num2str(c(2:5)))',...
        'AxesLabels',{'R_1','R_2','R_3','R_4','R_5'},'AxesLimits',repmat([0;maxy],1,5),...
        'AxesLabelsEdge','none','Color',mygrayred(3:end,:),'AxesOffset',0,...
        'AxesInterval',2,'LineWidth',1,'LabelFontSize',12,'AxesFontSize',8);
    h.PaperUnits = 'inches';
    h.PaperPosition = [0,0,1,1]*6;
    saveas(h, fullfile(plotdir,[filename, '.eps']),'epsc2');
   
    si = 2;
%     if si == 2 && bi == 2
%         maxy = 45;
%     end
    h = figure;
    filename = sprintf('Fig8cf_CodingFiveItems_WM_b%i_s%i',bi,si);
    s = spider_plot_class(squeeze(WM(bi,2:5,:)), 'LegendLabels',cellstr(num2str(c(2:5)))',...
        'AxesLabels',{'R_1','R_2','R_3','R_4','R_5'},'AxesLimits',repmat([0;maxy],1,5),...
        'AxesLabelsEdge','none','Color',mygrayblue(3:end,:),'AxesOffset',0,...
        'AxesInterval',2,'LineWidth',1,'LabelFontSize',12,'AxesFontSize',8);
    h.PaperUnits = 'inches';
    h.PaperPosition = [0,0,1,1]*6;
    saveas(h, fullfile(plotdir,[filename, '.eps']),'epsc2');
    
end