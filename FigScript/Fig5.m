% Fig 5. Winner-take-all competition
%% panel a, input structure
c = [3.2 12.8, 25.6, 38.4 51.2]'/100; % percentage of coherence
V = [1+c, 1-c]*scale0;
h = figure; hold on;
filename = 'Fig5a_InputMatrix';
% mygray = gray(length(c) + 3);
for vi = 1:length(V)
    Vinput = V(vi,:);
    plot(Vinput(2),Vinput(1),'k.','MarkerSize',mksz/2,'Color',mygray(vi+2,:));
end
ylabel('V_1 (a.u.)');xlabel('V_2 (a.u.)');
xlim([0,500]);ylim([0,500]);
xticks([0,500]);yticks([0,500]);
set(gca,'FontSize',fontsize-7);
set(gca,'TickDir','out');
set(gca,'LineWidth',1);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 3.1 2.8]/2.4;
saveas(h,fullfile(plotdir,[filename,'.eps']),'epsc2');
%% panel a & b, dynamic of neural firing rates
a = a0*eye(2);
b = b0*eye(2);
w = ones(2);
predur = .8;  
presentt = dt;
dur = 1.1;
stimdur = dur-presentt;
triggert =  presentt;
sgm = 0;
thresh = 70;
stoprule = 1;
Vprior = [1, 1]*scale0 + 0;
Tau = [.1,.1,.1];
initialvals = zeros(3,2);
for Nclass = ['R','G','D']
    h = figure; hold on;
    filename = sprintf('Fig5a_%s_latebeta',Nclass);
    for vi = 1:5
        Vinput = [1+c(vi), 1-c(vi)]*scale0 + 0;
        [choice, rt, R, G, D] = LDDM(Vprior, Vinput, w, a, b,...
            sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        lgd2(vi) = plot(eval([Nclass,'(:,2)']), 'k--', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
        lgd1(vi) = plot(eval([Nclass, '(:,1)']), 'k-', 'Color', mygray(vi+1,:), 'LineWidth',lwd);
    end
    plot([1.2, dur+predur]/dt,[thresh,thresh], 'k-');
    text(1200,thresh*1.05,'threshold');
    
    plot([1, 1, (predur)/dt],...
        [0, 2, 2]-7,'k--','LineWidth',lwd/2); % input, target, pre-motion
    plot([(predur)/dt, (predur+presentt+dur)/dt],...
        [2,2]-7,'k-','LineWidth',lwd/2); % inputs, stimuli, motion
    plot([1, (predur+triggert)/dt, (predur+triggert)/dt, (predur+presentt+dur)/dt],[0,0,2,2]-11.5,'k-','LineWidth',lwd/2); % disinhibition
    ylim([-16.5,71]);
    yticks([0, 70]);
    yticklabels({'0', '70'});
    ylabel('          Activity (a.u.)');
    xticks([1, predur/dt]);
    xticklabels({'',''});
    drawaxis(gca, 'x', 0, 'movelabel', 1);
    xlim([-50, (predur+presentt+dur)/dt]);
    legend([lgd1(5), lgd2(5)],{[Nclass,'_1'], [Nclass,'_2']},...
        'Location','NorthWest','FontSize',fontsize-5, 'FontName','Times New Roman', ...
        'FontAngle','italic','NumColumns',1,'Box','off');
    mysavefig(h, filename, plotdir, fontsize - 2, [2.8 2.54]);
end

%% panel c_left, nullclines for R1 and R2 under equal inputs
a = a0;
w = 1;
v = 1;
sgm = .02; dur = 12;
predur = 0;
presentt = dt; triggert = presentt; thresh = Inf;
initialvals = [4,4;8,8;0,0]/15; stimdur = dur; stoprule = 0;
filenamelist = {'Fig5cL','Fig5cM','Fig5cR'};
blist = [.9];
cplist = [0, .512, 1];
Vprior = [1, 1]*scale0 + 0;
for bi = 1
    b = blist(bi);
    for ci = 1:3
        cp = cplist(ci);
        rng('default'); rng(8);
        V = [1+cp 1-cp]*scale0 + 0;
        h = figure; hold on;
        filename = filenamelist{bi,ci};
        % - Nullclines R1*-R2* space
        R2 = 10.^(linspace(-1,3,300)); %linspace(.1,1000,300);
        R1 = (V(2)./R2 - (w - b)*R2 - (1-a))/v; % dR2/dt = 0
        lgd2 = plot(R2,R1,'k--','LineWidth',lwd/2); % dR2/dt = 0
        plot([min(R2), V(2)/(1-a)],[V(1)/(1-a), V(1)/(1-a)],'--k','LineWidth',1);
        Line1 = [R1' R2'];
        R1 = 10.^(linspace(-1,3,300)); %linspace(.1,1000,300);
        R2 = (V(1)./R1 - (w - b)*R1 - (1-a))/v; % dR1/dt = 0
        lgd1 = plot(R2,R1,'k-','LineWidth',lwd/2); % dR1/dt = 0
        plot([V(2)/(1-a), V(2)/(1-a)],[min(R1), V(1)/(1-a)],'--k','LineWidth',1);
        Line2 = [R1' R2'];
        % - time couse R1 & R2
        if 1
            [choice, rt, R, G, D] = LDDM(Vprior, V, [w v; v w], a*eye(2), b*eye(2),...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            ldgtrc(1) = plot(R(round(presentt/dt):5:end,2), R(round(presentt/dt):5:end,1),'-','Color','#ef476f','LineWidth',lwd/3);
            [choice, rt, R, G, D] = LDDM(Vprior, V, [w v; v w], a*eye(2), b*eye(2),...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            ldgtrc(2) = plot(R(round(presentt/dt):5:end,2), R(round(presentt/dt):5:end,1),'-','Color','#118ab2','LineWidth',lwd/3);
        end
        if 1 %b == 0 || V(1) ~= V(2)
            plot([0.5,1000],[0.5,1000],'--k');
        end
        % - Stability around fixed point(s)
        b = b*w;
        syms R1 R2
        eqns = [(V(1)/R1 - (w - b)*R1 - (1-a)) == R2*v, ... % dR1/dt = 0
            (V(2)/R2 - (w - b)*R2 - (1-a)) == R1*v];% dR2/dt = 0
        vars = [R1 R2];
        [AnswR1, AnswR2] = solve(eqns, vars);
        AnswD1 = b*AnswR1;
        AnswD2 = b*AnswR2;
        AnswG1 = w*AnswR1 + v*AnswR2 - AnswD1;
        AnswG2 = v*AnswR1 + w*AnswR2 - AnswD2;
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
        xticks([1,10,100,1000]);yticks([1,10,100,1000]);
        xlabel('R_2 activity (a.u.)');
        ylabel('R_1 activity (a.u.)');
        mysavefig(h, filename, plotdir, fontsize - 2, [2.8 2.54]*.8);
    end
end
%% panel d, equilibrium point for a noiseless system
a = 0; %a0;
b = blist(1);
w = 1;
v = 1;
cplist = linspace(-1,1,40);
Vprior = [1, 1]*scale0 + 0;
CodeRatio = nan(size(cplist));
name = sprintf('CodedRatio_WTA_LDDM_Sim%i_a%1.1f_b%1.1f',length(cplist), a, b);
output = fullfile(datadir,[name '.mat']);
if ~exist(output,'file')
    for ii = 1:length(cplist)
        fprintf('cp %3.1f',cplist(ii));
        fprintf('.');
        cp = [1 + cplist(ii), 1 - cplist(ii)];
        V = cp*scale0 + 0;
        syms R1 R2
        eqns = [(V(1)/R1 - (w - b)*R1 - (1-a))/v == R2, ... % dR1/dt = 0
            (V(2)/R2 - (w - b)*R2 - (1-a))/v == R1];% dR2/dt = 0
        vars = [R1 R2];
        [AnswR1,AnswR2] = solve(eqns, vars);
        AnswD1 = b*AnswR1;
        AnswD2 = b*AnswR2;
        AnswG1 = w*AnswR1 + v*AnswR2 - AnswD1;
        AnswG2 = v*AnswR1 + w*AnswR2 - AnswD2;
        PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0;
        Npoints = sum(PositiveRealAnsw);
        Ratio = double((AnswR1(PositiveRealAnsw)./(AnswR1(PositiveRealAnsw)+AnswR2(PositiveRealAnsw))) - .5);
        if sum(V == 0)
            Ratio = [Ratio,V(1)/sum(V)];
        end
        Cntrst = abs(Ratio);
        CodeRatio(ii) = Ratio(Cntrst == max(Cntrst))+.5;
        if V(1)/(V(1)+V(2)) == 1 || V(1)/(V(1)+V(2)) == 0
            CodeRatio(ii) = NaN;
        end
    end
    save(output,'CodeRatio','cplist','a','b','w','v');
else
    load(output);
end
h = figure;
filename = 'Fig5d'; hold on;
name = sprintf('CodedRatio_VR_LDDM_DNM_RNM');
output = fullfile(datadir,[name '.mat']);
load(output);
lgd2 = plot((1+cplist)/2, R1rp(:,1)./(R1rp(:,1)+R2rp(:,1)),'.-','Color',colorpalette{5},'MarkerSize',mksz/2);
lgd1 = plot((1+cplist)/2, CodeRatio, '.-','Color',colorpalette{3},'MarkerSize',mksz/2);
xlabel('Input ratio V_1/(V_1 + V_2)');ylabel('Coded ratio R_1/(R_1 + R_2)');
xticks([0,.25,.5,.75,1]); yticks([0,.25,.5,.75,1]);
legend([lgd1, lgd2],{'\color[rgb]{0.0235,0.8392,0.6275}Choice','\color[rgb]{0.0275,0.2314,0.2980}Representation'},'Box','off','Location','NorthWest','FontSize',fontsize-5);
mysavefig(h, filename, plotdir, fontsize - 2, [2.8 2.54]*.95);
%% panel e, parameter space for choice/representation under equal inputs
V = [1, 1]*scale0 + BR;
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
            AnswD1 = b*AnswR1;
            AnswD2 = b*AnswR2;
            AnswG1 = w*AnswR1 + v*AnswR2 - AnswD1;
            AnswG2 = v*AnswR1 + w*AnswR2 - AnswD2;
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
dyadic = visualize == 2 | visualize == 0;
h = figure;colormap(colorpalettergb([3,5],:));
filename = 'Fig5e';
s = surf(dyadic+1,'EdgeColor','none');
ylim([1,length(avec)]);
xlim([1,length(rvec)]);
xticks(linspace(1,length(rvec),5));
xticklabels(linspace(0,4,5));
yticks(linspace(1,length(avec),5));
% yticklabels(linspace(0,60,5));
yticklabels({'10^{-1}','10^0','10^1','10^2','10^3'});
xlabel('\beta');
ylabel('\alpha');
view(0,90);
mysavefig(h, filename, plotdir, fontsize - 2, [2.8 2.54]*.95);

