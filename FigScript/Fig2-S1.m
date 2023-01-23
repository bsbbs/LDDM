% Fig1-S2. Testing models
V = [30, 30];
w = 1;
a = 10;
b = .75;
gamma = 0;
eta = 0;
syms R1 R2
eqns = [(V(1)/R1 - (w - b)*R1 - (1-a)) == R2*(w + gamma + eta/R1), ... % dR1/dt = 0
    (V(2)/R2 - (w - b)*R2 - (1-a)) == R1*(w + gamma + eta/R2)];% dR2/dt = 0
vars = [R1 R2];
[AnswR1,AnswR2] = solve(eqns, vars);
PositiveRealAnsw = double(AnswR1) >0 & double(imag(AnswR1)) == 0 ...
    & double(AnswR2) >0 & double(imag(AnswR2)) == 0;
Npoints = sum(PositiveRealAnsw);
if ~exist(fullfile(plotdir, 'Fig2-S1Series'),'dir')
    mkdir(fullfile(plotdir, 'Fig2-S1Series'));
end
    
%% time course
dt = .001;
a = [1 0;
    0 1]*2;
eta = [0 1;
    1 0]*10;
w = [1 1;
    1 1];
r = [0 1;
    1 0]*10;
b = [1 0;
    0 1]*2.2;
tau = .1;
dur = 1.3/dt;
thresh = 12*6;
V = [35;25];
for acntrl = [0,1]
    for bcntrl = [0,1]
        h = figure;
        filename = fullfile('Fig2-S1Series',sprintf('a%ib%i.eps',acntrl,bcntrl));
        i = 0;
        for rcntrl = [0,1]
            for etacntrl = [0,1]
                i = i + 1;
                subplot(2,2,i); hold on;
                R0 = [0;0];
                G0 = [0;0];
                D0 = [0;0];
                I0 = [0;0];
                E0 = [0;0];
                R = R0;
                G = G0;
                D = D0;
                I = I0;
                E = E0;
                Rpool = [];
                for ti = 2:dur
                    dR = (-R + (V + acntrl*a*R - I)./(1 + G))*dt/tau;
                    dG = (-G + w*R + E - D)*dt/tau;
                    dD = (-D + bcntrl*b*R)*dt/tau;
                    dI = (-I + etacntrl*eta*R)*dt/tau;
                    dE = (-E + rcntrl*r*R)*dt/tau;
                    R = R + dR;
                    G = G + dG;
                    D = D + dD;
                    I = I + dI;
                    E = E + dE;
                    R(R<0) = 0;
                    G(G<0) = 0;
                    D(D<0) = 0;
                    I(I<0) = 0;
                    E(E<0) = 0;
                    Rpool(:,ti) = R;
                    if max(R) > thresh
                        break;
                    end
                end
                plot(Rpool(1,:),'LineWidth',lwd, 'Color',colorpalette{1});
                plot(Rpool(2,:),'LineWidth',lwd, 'Color',colorpalette{3});
                ylim([0,12]);
                xlim([-100,dur]);
                xticks([0,1000]);
                xticklabels({'0','1.0'});
                if i >= 3
                    xlabel('Time (a.u.)');
                end
                if mod(i,2) == 1
                    ylabel('Activity (a.u.)');
                end
                if acntrl == 0 && bcntrl == 0 && rcntrl == 0 && etacntrl == 0
                    legend({'R_1','R_2'},'FontName','Times New Roam',...
                    'FontAngle','italic','Location','Best', 'Box','off', 'FontSize',fontsize-5);
                end
                mysavefig(h, filename, plotdir, fontsize, aspect2);
            end
        end
    end
end