% Fig4. 
cplist = linspace(-1,1,121);
V1Iterp = (1 + cplist)*scale0;
V2Iterp = (1 + cplist)*scale0;
%% panel a, colormap of R1 activities as a function of V1 and V2 inputs
sgm = 0;
w = [1,1; 1,1];
predur = 0;
dur = 3.5;
presentt = 0;
stimdur = 2;
triggert = dur;
initialvals = [2,2;4,4]*0;
Vprior = [1, 1]*scale0 + B0;

% DNM
name = sprintf('DNMVR_dur%1.1f_W%1.2f%1.2f_sgm%1.2f',dur, w(1,1), w(1,2), sgm);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        for kk = 1:length(V2Iterp)
            fprintf('.');
            Vinput = [V1Iterp(ii), V2Iterp(kk)] + B0;
            [R, G] = DNM(Vinput, w, sgm, Tau, presentt, dur, dt, initialvals, stimdur);
            R1rp(kk,ii) = R(round((predur+presentt+stimdur)/dt) - 1,1);
            R2rp(kk,ii) = R(round((predur+presentt+stimdur)/dt) - 1,2);
            R1wm(kk, ii) = R(round((predur+triggert)/dt) - 1,1);
            R2wm(kk, ii) = R(round((predur+triggert)/dt) - 1,2);
        end
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm');
else
    load(output);
end
% - heatmap for value representation
h=figure;   colormap(jet);
filename = 'Fig4aM';
imagesc(V2Iterp,V1Iterp,R1rp'/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
caxis([0,1]);
xlabel('V_2 (a.u.)');ylabel('V_1 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefigs(h, filename, plotdir, fontsize - 2, [3.2, 2.3]);

% the LDDM model
a = a0*eye(2);
b = zeros(2);
name = sprintf('LDDMVR_dur%1.1f_W%1.2f%1.2f_alpha%1.1f_sgm%1.2f',dur, w(1,1), w(1,2), a(1), sgm);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        for kk = 1:length(V2Iterp)
            fprintf('.');
            Vinput = [V1Iterp(ii), V2Iterp(kk)] + B0;
            [choice, rt, R, G, I] = LDDM(Vprior, Vinput, w, a, b,...
                sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R1rp(kk,ii) = R(round((predur+presentt+stimdur)/dt) - 1,1);
            R2rp(kk,ii) = R(round((predur+presentt+stimdur)/dt) - 1,2);
            R1wm(kk, ii) = R(round((predur+triggert)/dt) - 1,1);
            R2wm(kk, ii) = R(round((predur+triggert)/dt) - 1,2);
        end
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm');
else
    load(output);
end
% - heatmap for value representation
h=figure;   colormap(jet);
filename = 'Fig4aL';
imagesc(V2Iterp,V1Iterp,R1rp'/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
caxis([0,1]);
xlabel('V_2 (a.u.)');ylabel('V_1 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefigs(h, filename, plotdir, fontsize - 2, [3.2, 2.3]);

% RNM (Wong and wang 2006)
presentt = 0;
dur = 12;
triggert = dur;
stimdur = 8;
miu0 = 30;
sgm = 0 ;
tauS = .1;
tauAMPA = .002;
I0 = .3255; % nA
gamma = .641;
JN = [.2609 -.0497
    -.0497   .2609]; % nA
initialvals = [2 2;.1 .1];
name = sprintf('WW06VR_dur%1.1f_stimdur%1.1f', dur, stimdur);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        for kk = 1:length(V2Iterp)
            fprintf('.');
            Vinput = [V1Iterp(ii), V2Iterp(kk)];
            cp = Vinput/250;
            [choice, rt, nu_wind, s_wind, H, S] = wong06(cp,miu0,sgm,I0,JN,...
                gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
            R1rp(kk,ii) = nu_wind(round((presentt+stimdur)/dt) - 1,1);
            R2rp(kk,ii) = nu_wind(round((presentt+stimdur)/dt) - 1,2);
            R1wm(kk, ii) = nu_wind(round((triggert)/dt) - 1,1);
            R2wm(kk, ii) = nu_wind(round((triggert)/dt) - 1,2);
        end
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm');
else
    load(output);
end
% - heatmap for value representation
h=figure;   colormap(jet);
filename = 'Fig4aR';
imagesc(V2Iterp,V1Iterp,R1rp'/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
caxis([0,1]);
xlabel('V_2 (a.u.)');ylabel('V_1 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefigs(h, filename, plotdir, fontsize - 2, [3.2, 2.3]);

%% Fit Louie et al., 2011
% You need to do the following adaptation based on your own environment in order to run the
% 1. re-arrange the related directories based on your own 
% 2. download the optimization algorithm "bads" from https://github.com/acerbilab/bads
cd ../Fit;
FitLouie2011_DNM;
FitLouie2011_LDDM;
FitLouie2011_RNM;
