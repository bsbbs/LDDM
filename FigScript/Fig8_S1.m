%% Fig S3 - Persistent activity under Asymmetric gain control weights

%% phase plane analysis and vector field. i = 1: v<w; i = 2: v=w; i = 3: v>w
lgdvec = [];
for wi = 1:3
    switch wi
        case 1
        filename = 'Fig8-S1a';
        case 2
        filename = 'Fig8-S1b';
        case 3
        filename = 'Fig8-S1c';
    end
    Taup = [1,1,1];
    alpha = eye(2)*10;
    a = alpha(1,1);
    if wi == 1
        w = [1,.7;.7,1];
    elseif wi == 2
        w = [1,1;1,1];
    elseif wi == 3
        w = [1,3;3,1];
    end
    beta = 0; % 0, .4, 1.2
    Dt = 1;
    if w(1,1) == w(1,2)
        if beta > 0 && beta < 1
            axismax = 25;
        elseif beta == 0
            axismax = 10;
        elseif beta >= 1
            axismax = 50;
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
    G1 = R10*w(1,1) + R20*w(1,2) - beta*R10;
    G2 = R10*w(2,1) + R20*w(2,2) - beta*R20;
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
    if w(2,2) - beta > 0 && w(1,2) == w(1,1)
        R2Line = [0,(alpha(2,2) - 1)/(w(2,2) - beta)];
    elseif w(2,2) - beta <= 0 && w(1,2) == w(1,1)
        R2Line = [0,axismax];
    elseif beta == 0 && w(1,2) ~= w(1,1)
        R2Line = [0,(alpha(1,1)-1)/w(1,2)];
    end
    R1Line = -w(1,2)/(w(1,1)-beta).*R2Line+(alpha(1,1) - 1)/(w(1,1)-beta);
    lg1 = plot(R2Line,R1Line,'-','Color',colorpalette{4},'LineWidth',lwd*2);% dR1/dt = 0
    if w(2,2) - beta > 0 && w(1,2) == w(1,1)
        R1Line = [0,(alpha(1,1) - 1)/(w(1,1) - beta)];
    elseif w(2,2) - beta <= 0 && w(1,2) == w(1,1)
        R1Line = [0,axismax];
    elseif beta == 0 && w(1,2) ~= w(1,1)
        R1Line = [0,(alpha(2,2)-1)/w(2,1)];
    end
    R2Line = -w(2,1)/(w(2,2)-beta).*R2Line+(alpha(2,2) - 1)/(w(2,2)-beta);
    lg2 = plot(R2Line,R1Line,'--','Color',colorpalette{1},'LineWidth',lwd*2);% dR2/dt = 0
    syms R1 R2;
    vars = [R1 R2];
    eqns = [(w(1,1)-beta)*R1 + w(1,2)*R2 == alpha(1,1) - 1, ... % dR1/dt = 0
        w(2,1)*R1 + (w(2,2)-beta)*R1 == alpha(2,2) - 1];% dR2/dt = 0
    [AnswR1,AnswR2] = solve(eqns, vars);
    AnswI1 = beta*AnswR1;
    AnswI2 = beta*AnswR2;
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
            beta, 0, -1, 0, 0, 0
            0, 0, 0, -1+a/(1+G2star),   -(a*R2star)/(1+G2star)^2, 0
            w(2,1), 0, 0, w(2,2), -1, -1
            0, 0, 0, beta, 0, -1];
        A = eig(JMat);
        if all(real(A) < -1e-10)
            lgdvec(end +1) = plot(R2star, R1star, '.', 'Color','#06d6a0', 'MarkerSize',mksz);
            plot(R2star, R1star, 'o', 'Color','k', 'MarkerSize',mksz/3);
        elseif sum(real(A) > 1e-10)
            lgdvec(end +1) = plot(R2star, R1star, '.', 'Color', '#ffd166', 'MarkerSize',mksz);
            plot(R2star, R1star, 'o', 'Color','k', 'MarkerSize',mksz/3);
        end
    end
    if w(1,1) == w(1,2)
        if beta > 0 && beta < 1
            plot([0, (alpha(2,2)-1)/(w(2,2) - beta)],[(alpha(1,1)-1)/(w(1,1) - beta),(alpha(1,1)-1)/(w(1,1) - beta)],'k--','LineWidth',1); % constraints
            plot([(alpha(2,2)-1)/(w(2,2) - beta),(alpha(2,2)-1)/(w(2,2) - beta)],[0, (alpha(1,1)-1)/(w(1,1) - beta)],'k--','LineWidth',1); % constraints
            xlim([-.4,axismax]);ylim([-.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1)+w(1,2)-beta), (alpha(1,1)-1)/(w(1,1) - beta)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2)+w(2,1)-beta), (alpha(2,2)-1)/(w(2,2) - beta)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$','$\frac{\alpha-1-G_0}{\omega-\beta}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$','$\frac{\alpha-1-G_0}{\omega-\beta}$'});
        elseif beta == 0
            xlim([-.4,axismax]);ylim([-.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) - beta)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) - beta)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{\omega}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{\omega}$'});
        elseif beta >= 1
            xlim([-.4,axismax]);ylim([-.4,axismax]);
            yticks([0, (alpha(1,1)-1)/(w(1,1) + w(1,2) - beta)]);
            xticks([0, (alpha(2,2)-1)/(w(2,2) + w(2,1) - beta)]);
            yticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$'});
            xticklabels({'0','$\frac{\alpha-1-G_0}{2\omega-\beta}$'});
        end
    elseif w(1,1) ~= w(1,2) && beta == 0
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
            xticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{w}$'});
            yticklabels({'0','$\frac{\alpha-1-G_0}{w+v}$','$\frac{\alpha-1-G_0}{w}$'});
        end
    end
    set(gca,'TickLabelInterpreter', 'latex');
    xlabel('R_2 activity (a.u.)', 'FontAngle', 'italic');
    ylabel('R_1 activity (a.u.)', 'FontAngle', 'italic');
    lg3 = quiver(R20,R10,dR2,dR1,1.3,'k','LineWidth',1);
    ratio = [2.9, 2.54];%aspect3;
    if wi == 3
        legend([lg1,lg2,lg3,lgdvec], {'\color[rgb]{.0667,.5412,.6980}dR_1/dt = 0',...
            '\color[rgb]{.9373,.2784,.4353}dR_2/dt = 0','\color{black}Change rate',...
            'Stable point','Unstable point'}, 'FontName','Times New Roman', ...
            'FontAngle','Italic','FontSize',fontsize-4, 'Location','northeastoutside','Box','off');
        ratio = aspect5;
    end
    mysavefig(h, filename, plotdir, 11, ratio);
end

%% example dynamic of R1 and R2. i = 1: v<w; i = 2: v=w; i = 3: v>w
c = [3.2 12.8, 25.6, 38.4 51.2]'/100; % percentage of coherence
VmatDiag = 30*[1+c, 1-c];
for wi = 1:3
    switch wi
        case 1
            filename = 'Fig8-S1d';
        case 2
            filename = 'Fig8-S1e';
        case 3
            filename = 'Fig8-S1f';
    end
    if wi == 1
        w = [1,.7;.7,1];
    elseif wi == 2
        w = [1,1;1,1];
    elseif wi == 3
        w = [1,3;3,1];
    end
    
    a = 10*eye(2,2);
    b = 0*eye(2,2);
    initialvals = zeros(3,2);
    presentt = .12;
    dur = 1.5;
    stimdur = .45-presentt;
    sgm = 0;
    triggert = Inf;
    thresh = Inf;
    stoprule = 0;
    h = figure; hold on;
    %Vmat = [[1 15 45 60]' [30,30,30,30]'];
    for vi = 2:5
        Vinput = VmatDiag(vi,:);
        [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        %subplot(2,1,1); hold on;
        lgd1 = plot(R(:,1), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
        %subplot(2,1,2); hold on;
        lgd2 = plot(R(:,2), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
    end
%     %subplot(2,1,1);
%     ylim([-.4,max(ylim)*1.1]);
%     yticks([0, 30]);
%     xticks([presentt/dt, (presentt+stimdur)/dt]);
%     xticklabels({});
%     savefig(h, filename, plotdir, fontsize, aspect2);
%     subplot(2,1,2);
    ylim([-.4,max(ylim)*1.1]);
    yticks([0, 30]);
    xticks([presentt/dt, (presentt+stimdur)/dt]);
    xticklabels({'on','off'});
    ylabel('Activity (a.u.)');
    xlabel('Time (a.u.)');
    if wi == 1
    legend([lgd1, lgd2],{'R_1', 'R_2'},...
        'Location','North','FontSize',fontsize-5, 'FontName','Times New Roman', ...
        'FontAngle','italic','NumColumns',1,'Box','off');
    end
    mysavefig(h, filename, plotdir, fontsize, aspect10);
end

%% Input
h = figure; hold on;
for vi = 2:5
    Vinput = VmatDiag(vi,:);
    plot(Vinput(2),Vinput(1),'k.','MarkerSize',mksz/2,'Color',mygray(vi+1,:));
end
xlabel('V_2 (a.u.)');ylabel('V_1 (a.u.)');
xlim([0,60]);ylim([0,60]);
xticks([0,60]);yticks([0,60]);
set(gca,'FontSize',fontsize-6);
set(gca,'TickDir','out');
set(gca,'LineWidth',1);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 3.1 2.8]/3;
saveas(h,fullfile(plotdir,sprintf('ValueInputMatrix_Fig8_S1.eps')),'epsc2');