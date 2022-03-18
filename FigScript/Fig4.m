% Fig4. 
%% panel c, colormap of R1 activities as a function of V1 and V2 inputs
sigma = 0;
w = [1,1; 1,1];
dur = 3.5;
presentt = 0;
stimdur = 3;
triggert = dur;
% Lofaro's
alpha = eye(2,2)*0;
name = sprintf('Lofaro_dur%1.1f_W%1.2f%1.2f_alpha%1.1f_sgm%1.2f',dur, w(1,1), w(1,2), alpha(1), sigma);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        for kk = 1:length(V2Iterp)
            fprintf('.');
            Vinput = [V1Iterp(ii), V2Iterp(kk)];
            [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sigma, ...
                Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R1rp(kk,ii) = R(round((presentt+stimdur)/dt) - 1,1);
            R2rp(kk,ii) = R(round((presentt+stimdur)/dt) - 1,2);
            R1wm(kk, ii) = R(round((triggert)/dt) - 1,1);
            R2wm(kk, ii) = R(round((triggert)/dt) - 1,2);
        end
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm');
else
    load(output);
end
% - heatmap for value representation
h=figure;   colormap(jet);
filename = 'Fig3cL';
imagesc(V1Iterp,V2Iterp,R1rp/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
xlabel('V_1 (a.u.)');ylabel('V_2 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefig(h, filename, plotdir, fontsize, aspect4);
% the hybrid model
alpha = eye(2,2)*15;
name = sprintf('LcDVR_dur%1.1f_W%1.2f%1.2f_alpha%1.1f_sgm%1.2f',dur, w(1,1), w(1,2), alpha(1), sigma);
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];R1wm = [];R2wm = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        for kk = 1:length(V2Iterp)
            fprintf('.');
            Vinput = [V1Iterp(ii), V2Iterp(kk)];
            [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sigma, ...
                Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R1rp(kk,ii) = R(round((presentt+stimdur)/dt) - 1,1);
            R2rp(kk,ii) = R(round((presentt+stimdur)/dt) - 1,2);
            R1wm(kk, ii) = R(round((triggert)/dt) - 1,1);
            R2wm(kk, ii) = R(round((triggert)/dt) - 1,2);
        end
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','R1wm','R2wm');
else
    load(output);
end
% - heatmap for value representation
h=figure;   colormap(jet);
filename = 'Fig3cM';
imagesc(V1Iterp,V2Iterp,R1rp/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
xlabel('V_1 (a.u.)');ylabel('V_2 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefig(h, filename, plotdir, fontsize, aspect4);
%% Wong and wang 2006 model
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
            [nu_wind, s_wind, rt, choice, H, S] = wong06(Vinput,miu0,sgm,I0,JN,...
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
filename = 'Fig3cR';
imagesc(V1Iterp,V2Iterp,R1rp/max(max(R1rp)));
xticks([0,500]);yticks([0,500]);
set(gca,'YDir','normal');
xlabel('V_1 (a.u.)');ylabel('V_2 (a.u.)');
c = colorbar;
ylabel(c, 'Rescaled R_1 activity');
xh = get(gca, 'xlabel'); p = get(xh, 'position'); p(2) = p(2)+abs(p(2))*.9; set(xh, 'position',p);
yh = get(gca, 'ylabel'); p = get(yh, 'position'); p(1) = p(1)+abs(p(1))*.8; set(yh, 'position',p);
savefig(h, filename, plotdir, fontsize, aspect4);

%% panel d, the coded ratio and input ratio under different connection weights
sigma = 0;alpha = eye(2)*15;beta = zeros(2);
dur = 10;presentt = dt;stimdur = 10;triggert = dur;stoprule = 0;
name = sprintf('CodedRatio_LDDM_DNM_RNM');
output = fullfile(datadir,[name '.mat']);
R1rp = [];R2rp = [];
if ~exist(output,'file')
    for ii = 1:length(V1Iterp)
        fprintf('V1 %3.1f',V1Iterp(ii));
        fprintf('.');
        initialvals = zeros(3,2);
        Vinput = [V1Iterp(ii), 512 - V1Iterp(ii)];
        [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sigma, ...
            Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        R1rp(ii,1) = R(round((presentt+stimdur)/dt) - 1,1);
        R2rp(ii,1) = R(round((presentt+stimdur)/dt) - 1,2);
        [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha*0, beta, sigma, ...
            Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        R1rp(ii,2) = R(round((presentt+stimdur)/dt) - 1,1);
        R2rp(ii,2) = R(round((presentt+stimdur)/dt) - 1,2);

        initialvals = [2 2;.1 .1];
        [nu_wind, s_wind, rt, choice, H, S] = wong06(Vinput,miu0,sgm,I0,JN,...
            gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
        R1rp(ii,3) = nu_wind(round((presentt+stimdur)/dt) - 1,1);
        R2rp(ii,3) = nu_wind(round((presentt+stimdur)/dt) - 1,2);
        fprintf('\n');
    end
    save(output,'R1rp','R2rp','V1Iterp');
else
    load(output);
end
%%
h = figure;
hold on;
filename = 'Fig3d';
i = 1; % LDDM
interval = 3;
lgd(i)=plot(V1Iterp(1:interval:end)/512, R1rp(1:interval:end,i)./(R1rp(1:interval:end,i)+R2rp(1:interval:end,i)),'.-','Color',colorpalette{5},'MarkerSize',mksz/2);
i = 2; % DNM, skip
i = 3;
lgd(i)=plot(V1Iterp(1:interval:end)/512, R1rp(1:interval:end,i)./(R1rp(1:interval:end,i)+R2rp(1:interval:end,i)),'.-','Color',colorpalette{1},'MarkerSize',mksz/2);
legend([lgd(1),lgd(3)],{'\color[rgb]{.0275,.2314,.298} Normalized coding',...
    '\color[rgb]{.9373,.2784,.4353}WTA competition'}, 'Box','off', 'Location','Southeast',...
    'FontSize', fontsize-7, 'FontName','Times New Roman','FontAngle','Italic',...
    'FontWeight','bold');
xlabel('Input ratio V_1/(V_1 + V_2)');ylabel('Coded ratio R_1/(R_1 + R_2)');
xticks([0,.25,.5,.75,1]); yticks([0,.25,.5,.75,1]);
savefig(h, filename, plotdir, fontsize, aspect3);
%% old panel d
sigma = 0;alpha = eye(2)*15;beta = zeros(2);
dur = 4;presentt = dt;stimdur = 2;triggert = dur;stoprule = 0;
h = figure;
hold on;
filename = 'Fig3d';
for i = 1:5
    if i == 1
        w = [1,.1;.1,1]; %w = w/sqrt(sum(sum(w.^2))/4);
    elseif i == 2
        w = [1,.8;.8,1]; %w = w/sqrt(sum(sum(w.^2))/4);
    elseif i == 5
        w = [1,1;1,1]; %w = w/sqrt(sum(sum(w.^2))/4);
    elseif i == 3
        w = [.8,1;1,.8]; %w = w/sqrt(sum(sum(w.^2))/4);
    elseif i == 4
        w = [.3,1;1,.3]; %w = w/sqrt(sum(sum(w.^2))/4);
    end
    name = sprintf('CodedRatioSim%i_dur%1.1f_W%1.2f%1.2f_alpha%1.1f_sgm%1.2f',length(V1Iter), dur, w(1,1), w(1,2), alpha(1), sigma);
    output = fullfile(datadir,[name '.mat']);
    R1rp = [];R2rp = [];R1wm = [];R2wm = [];
    if ~exist(output,'file')
        for ii = 1:length(V1Iter)
            fprintf('V1 %3.1f',V1Iter(ii));
            fprintf('.');
            Vinput = [V1Iter(ii), V2Iter(ii)];
            [R, G, I, rt, choice] = LcDsInhbt(Vinput, w, alpha, beta, sigma, ...
                Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
            R1rp(ii) = R(round((presentt+stimdur)/dt) - 1,1);
            R2rp(ii) = R(round((presentt+stimdur)/dt) - 1,2);
            R1wm(ii) = R(round((triggert)/dt) - 1,1);
            R2wm(ii) = R(round((triggert)/dt) - 1,2);
            fprintf('\n');
        end
        save(output,'R1rp','R2rp','R1wm','R2wm','V1Iter','V2Iter');
    else
        load(output);
    end
    if i == 5
        lgd(i)=plot(V1Iter(:)./(512), R1rp(:)./(R1rp(:)+R2rp(:)),'.-','Color',colorpalette{i},'MarkerSize',mksz/2);
        %plot(V1Iter(:)./(512), R1rp(:)./(R1rp(:)+R2rp(:)),'-','Color',colorpalette{i});
    else
        lgd(i)=plot(V1Iter(:)./(512), R1rp(:)./(R1rp(:)+R2rp(:)),'.-','Color',colorpalette{i},'MarkerSize',mksz/2);
        %plot(V1Iter(:)./(512), R1rp(:)./(R1rp(:)+R2rp(:)),'-','Color',colorpalette{i});
    end
end
legend([lgd(1),lgd(2),lgd(5),lgd(3),lgd(4)],{'\color[rgb]{.9373,.2784,.4353}w_{ij} << w_{ii}','\color[rgb]{1,.8196,.4}w_{ij} < w_{ii}',...
    '\color[rgb]{.0275,.2314,.298}w_{ij} = w_{ii}','\color[rgb]{.0235,.8392,.6275}w_{ij} > w_{ii}',...
    '\color[rgb]{.0667,.5412,.6980}w_{ij} >> w_{ii}'}, 'Box','off', 'Location','Best',...
    'FontSize', fontsize-7, 'FontName','Times New Roman','FontAngle','Italic',...
    'FontWeight','bold');
xlabel('Input ratio (V_1 vs. V_2)');ylabel('Coded ratio (R_1 vs. R_2)');
xticks([0,.25,.5,.75,1]); yticks([0,.25,.5,.75,1]);
savefig(h, filename, plotdir, fontsize, aspect3);