%% E-E only, self-excitation & input

%% Property 1. integration of noise
potentiation = [1, 2];
predur = 0;
dur = 100;
triggert = Inf;
stimdur = Inf;
presentt = dt;
Vprior = [1, 1]*scale;
Vinput = [1+c, 1-c]*scale;
R_conds = [];
for pi = 1:numel(potentiation)
    STDP_v = potentiation(pi);
    STDP_a = potentiation(pi);
    STDP_G = 1;
    rng(5);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInput_STDP(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    R_conds{pi} = R;
end
%
h = figure; 
filename = 'Output_E-E';
subplot(2,1,1); hold on;
h1 = histogram(R_conds{1}(2000:end,1),100);
h1.FaceColor = 'r'; % myred(6,:);
h1.EdgeColor = 'n';
h2 = histogram(R_conds{1}(2000:end,2),100);
h2.FaceColor = 'b'; % myblue(6,:);
h2.EdgeColor = 'n';
title('Baseline');
xlabel(' ');
ylabel(' ');
yticks([]);
xlim([15, 40]);
savefigs(h, filename, plotdir, fontsize, [2, 3]);

subplot(2,1,2); hold on;
h3 = histogram(R_conds{2}(2000:end,1),100);
h3.FaceColor = 'r'; % myred(6,:);
h3.EdgeColor = 'n';
h4 = histogram(R_conds{2}(2000:end,2),100);
h4.FaceColor = 'b'; % myblue(6,:);
h4.EdgeColor = 'n';
title('STDP');
xlabel('Firing rates (Hz)');
ylabel('Frequency');
yticks([]);
xlim([15, 40]);
savefigs(h, filename, plotdir, fontsize, [2, 3]);
%% Property 2. predicted choice behavior
potentiation = [1, 2];
task = 'RT_E-E';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
sgm = .01;
sgmInput = 1/3;
w = w0*ones(2);
a = a0*eye(2);
b = b0*eye(2);
thresh = 70;
stoprule = 1;
initialvals = [R0,R0; G0,G0; D0,D0];
filename = sprintf('LDDM_%s_STDPfrom%1.1fto%1.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',...
    task,min(potentiation), max(potentiation), a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
STDP_v = 1;
STDP_a = 1;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    [V1, ~] = meshgrid(scale*(1+cp), potentiation);
    [V2, STDP_G] = meshgrid(scale*(1-cp), potentiation);
    Vinput.V1 = V1;
    Vinput.V2 = V2;
    Vprior.V1 = ones(size(V1))*scale;
    Vprior.V2 = ones(size(V2))*scale;
    [rt, choice, ~] = LDDM_Rndinput_STDP_GPU(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    ACC = gather(mean(2-squeeze(choice),3,'omitnan'));
    meanRT = gather(mean(squeeze(rt),3,'omitnan'));
    save(output,'ACC','meanRT');
else
    load(output);
end

mygray = [.9:-.9/numel(potentiation):0];
h = figure;
filename = sprintf('RT_ACC_%s',name);
subplot(2,1,1); hold on;
for level = 1:numel(potentiation)
    plot(cp,ACC(level,:),'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz,'LineWidth',lwd);
end
set(gca,'XScale','log');
ylabel('% Correct');
lgd = legend(cellstr(string(potentiation)),...
    'Location','SouthEast', 'FontSize', fontsize-2, 'Box','off');
title(lgd,'iSTDP','FontSize',fontsize-8);
savefigs(h, filename, plotdir,fontsize, [2, 3]);
subplot(2,1,2); hold on;
for level = 1:numel(potentiation)
    plot(cp,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz,'LineWidth',lwd);
end
set(gca,'XScale','log');
ylabel('RT (s)');
savefigs(h, filename, plotdir,fontsize, [2, 3]);

%% Property 3. ACC change over noise
