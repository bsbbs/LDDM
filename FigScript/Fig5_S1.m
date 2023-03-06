% Fig5-S1 - detailed phase plane analysis on the LDDM

%% colormap of a and beta/w space
V = [1, 1]*scale0;
w = 1;
v = 1;
V1 = V(1);
V2 = V(2);
%avec = linspace(0,60,601);
avec = 10.^[-1:.01:3];
rvec = linspace(0,4,401); % ratio of beta/w
Nb = length(rvec);
filename = sprintf('Stability_Numeric_V%1.0f_%1.0f_w%1.1f_v%1.1f_a%i_%1.2f_%1.2fr%i_%1.2f_%1.2f',...
    [V, w, v, length(avec),min(avec),max(avec),length(rvec),min(rvec),max(rvec)]);
output = fullfile(datadir,[filename, '.mat']);
if ~exist(output, 'file')
    Npoints = NaN(length(avec),length(rvec));
    Eigenvalue = -ones(length(avec),length(rvec),4);
    Stability = NaN(length(avec),length(rvec));
    for ai = 1:length(avec)
        a = avec(ai);
        fprintf('%2.1f.',a);
        syms R1 R2;
        vars = [R1 R2];
        parfor bi = 1:Nb
            fprintf('.');
            b = rvec(bi)*w;
            eqns = [(V1/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
                (V2/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
            [AnswR1,AnswR2] = solve(eqns, vars);
            AnswI1 = b*AnswR1;
            AnswI2 = b*AnswR2;
            AnswG1 = w*AnswR1 + v*AnswR2 - AnswI1;
            AnswG2 = v*AnswR1 + w*AnswR2 - AnswI2;
            PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
            Npoints(ai,bi) = sum(PositiveRealAnsw);
            NegEigen = -ones(1,4);
            for i = find(PositiveRealAnsw)'
                R1star = double(AnswR1(i));
                R2star = double(AnswR2(i));
                G1star = double(AnswG1(i));
                G2star = double(AnswG2(i));
                JMat = [-1 + a/(1+G1star), -(V1+a*R1star)/(1+G1star)^2, 0, 0, 0, 0
                    w, -1, -1, v, 0, 0
                    b, 0, -1, 0, 0, 0
                    0, 0, 0, -1+a/(1+G2star),   -(V2+a*R2star)/(1+G2star)^2, 0
                    v, 0, 0, w, -1, -1
                    0, 0, 0, b, 0, -1];
                A = eig(JMat);
                NegEigen(i) = all(real(A) < 0); % attractive 1, diversive 0
            end
            Stability(ai,bi)  = all(NegEigen>0);
            Eigenvalue(ai,bi,:) = NegEigen;
        end
        fprintf('\n');
    end
    visualize = sum(Eigenvalue,3);
    visualize(visualize == -2) = 2;
    save(output,'avec','rvec','V','w','v','Npoints','Eigenvalue','Stability','visualize');
else
    load(output);
end
% plot
visualize(visualize == -4) = -2;
visualize(visualize == -3) = -1;
h = figure;colormap(colorpalettergb([1,2,3,4,5],:));
filename = 'Fig5S1a';
s = surf(visualize,'EdgeColor','none');
ylim([1,length(avec)]);
xlim([1,length(rvec)]);
xticks(linspace(1,length(rvec),5));
xticklabels(linspace(0,4,5));
yticks(linspace(1,length(avec),5));
yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel('\beta');
ylabel('\alpha');
view(0,90);
mysavefig(h, filename, plotdir, fontsize, [3, 2.54]*1.5);

%% nullclines for R1 and R2 under unequal inputs, manually change parameters
c = 0;
V = scale0*[1+c, 1-c];
Vprior = scale0*[1, 1];
for spots = 3
    switch spots
        case 1
            a = 5;
            b = .3;
            filename = 'Fig5S1b';
        case 2
            a = .05;
            b = 1.003;
            filename = 'Fig5S1c';
        case 3
            a = 30;
            b = .9;
            filename = 'Fig5S1d';
        case 4
            a = 0;
            b = 1.5;
            filename = 'Fig5S1e';
        case 5
            a = 0;
            b = 2.5;
            filename = 'Fig5S1f';
    end
    w = 1;
    v = 1;
    % - time couse R1 & R2
    sgm = .4; dur = 5; predur = 0;
    presentt = dt; triggert = presentt; thresh = Inf;
    initialvals = [4,4;8,8;0,0]/15; stimdur = dur; stoprule = 0;
    % - Nullclines R1*-R2* space
    h = figure; hold on;
    R2 = [linspace(.2,50,400), 51:5:1000];
    R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
    lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
    plot([min(R2), V(2)/(1-a)],[V(1)/(1-a), V(1)/(1-a)],'--k','LineWidth',1);
    Line1 = [R1' R2'];
    R1 = [linspace(.2,50,400), 51:5:1000];
    R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
    lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
    plot([V(2)/(1-a), V(2)/(1-a)],[min(R1), V(1)/(1-a)],'--k','LineWidth',1);
    Line2 = [R1' R2'];
    if 1
        rng(7,'twister');
        [choice, rt, R, G, D] = LDDM(Vprior, V, [w v; v w], a*eye(2), b*eye(2),...
            sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        ldgtrc(1) = plot(R(round(presentt/dt):5:end,2), R(round(presentt/dt):5:end,1),'-','Color','#ef476f','LineWidth',lwd/3);
        [choice, rt, R, G, D] = LDDM(Vprior, V, [w v; v w], a*eye(2), b*eye(2),...
            sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        ldgtrc(2) = plot(R(round(presentt/dt):5:end,2), R(round(presentt/dt):5:end,1),'-','Color','#118ab2','LineWidth',lwd/3);
    end
    if 1 %b == 0 || V(1) ~= V(2)
        plot([0.1,1000],[0.1,1000],'--k');
    end
    b = b*w;
    syms R1 R2
    eqns = [(V(1)/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
        (V(2)/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
    vars = [R1 R2];
    [AnswR1,AnswR2] = solve(eqns, vars);
    AnswI1 = b*AnswR1;
    AnswI2 = b*AnswR2;
    AnswG1 = w*AnswR1 + v*AnswR2 - AnswI1;
    AnswG2 = v*AnswR1 + w*AnswR2 - AnswI2;
    PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
    Npoints = sum(PositiveRealAnsw);
    lgdvec = [];
    for i = [find(PositiveRealAnsw)]'
        R1star = double(AnswR1(i));
        R2star = double(AnswR2(i));
        G1star = double(AnswG1(i));
        G2star = double(AnswG2(i));
        JMat = [-1 + a/(1+G1star), -(V(1)+a*R1star)/(1+G1star)^2, 0, 0, 0, 0
            w, -1, -1, v, 0, 0
            b, 0, -1, 0, 0, 0
            0, 0, 0, -1+a/(1+G2star),   -(V(2)+a*R2star)/(1+G2star)^2, 0
            v, 0, 0, w, -1, -1
            0, 0, 0, b, 0, -1];
        A = eig(JMat);
        Stability(i) = all(real(A) < 0); % attractive 1, diversive 0
        if Stability(i)
            lgdvec(end +1) = plot(R2star, R1star, '.', 'Color','#06d6a0', 'MarkerSize',mksz/1.5);
            plot(R2star, R1star, 'o', 'Color','k', 'MarkerSize',mksz/4);
        else
            lgdvec(end +1) = plot(R2star, R1star, '.', 'Color', '#ffd166', 'MarkerSize',mksz/1.5);
            plot(R2star, R1star, 'o', 'Color','k', 'MarkerSize',mksz/4);
        end
    end
    set(gca, 'XScale','log');
    set(gca, 'YScale','log');
    xlim([.2,10^3]);ylim([.2,10^3]);
    if max(V)/(1-a) < 1000 && max(V)/(1-a) > 100
        xticks([1,10,100,V(1)/(1-a),1000]);yticks([1,10,100,V(1)/(1-a),1000]);
        xticklabels({'$10^0$','$10^1$','','$\frac{V_1+B_R}{1+B_G-\alpha}$',''});
        yticklabels({'$10^0$','$10^1$','','$\frac{V_2+B_R}{1+B_G-\alpha}$','$10^3$'});
        set(gca,'TickLabelInterpreter', 'latex');
        ratio = [3.1, 2.54];
    else
        xticks([1,10,100,1000]);yticks([1,10,100,1000]);
%         xticklabels({'$10^0$','$10^1$','$10^2$','$\frac{V_1+B_R}{1+B_G-\alpha}$','$10^3$'});
%         yticklabels({'$10^0$','$10^1$','$10^2$','$\frac{V_2+B_R}{1+B_G-\alpha}$','$10^3$'});
%         xticklabels({'$10^0$','$10^1$','$\frac{V_1+B_R}{1+B_G-\alpha}$','$10^2$'});
%         yticklabels({'$10^0$','$10^1$','$\frac{V_2+B_R}{1+B_G-\alpha}$','$10^2$'});
        set(gca,'TickLabelInterpreter', 'latex');
        ratio = [2.8, 2.54];
    end
    xlabel('R_2 activity (a.u.)', 'FontAngle','italic');
    ylabel('R_1 activity (a.u.)', 'FontAngle','italic');
    if spots == 1
        legend([lgd1, lgd2],{'dR_1/dt = 0','dR_2/dt = 0'}, 'FontName','Times New Roman', ...
            'FontAngle','Italic','FontSize',fontsize-5, 'Location','east','Box','off');
    end
    mysavefig(h, filename, plotdir, fontsize, ratio);
end