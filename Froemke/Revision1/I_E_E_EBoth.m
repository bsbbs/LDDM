%% I-E and E-E both

%% Property 1. integration of noise
potentiation = [1, 2];
predur = 0;
dur = 100;
sgm = .01;
sgmInput = sgmInput_rprsnt;
c = c_rprsnt;
a = 0*eye(2);
b = b0*eye(2);
w = ones(2);
triggert = Inf;
stimdur = Inf;
presentt = dt;
initialvals = [1,1;2,2;0,0]*eqlb/3;
Vprior = [1, 1]*scale;
Vinput = [1+c, 1-c]*scale;
filename = sprintf('LDDM_timeCourse_I-E&E-E_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f',a0,b0,sgm,sgmInput);
output = fullfile(Simdir,[filename, '.mat']);
if ~exist(output, 'file')
    R_conds = [];
    for pi = 1:numel(potentiation)
        STDP_V = potentiation(pi);
        STDP_a = potentiation(pi);
        STDP_G = potentiation(pi);
        rng(5);
        [choice, rt, R, G, D, Vcourse] = LDDM_RndInput_STDP(Vprior, Vinput, STDP_V, STDP_a, STDP_G, w, a, b,...
            sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
        R_conds{pi} = R;
    end
    save(output, 'R_conds','potentiation');
else
    load(output);
end
% smoothed curve
pd1 = fitdist(R_conds{1}(2000:end,1),'kernel','Kernel','normal');
x = 10:.1:38;
y1 = pdf(pd1,x);
pd2 = fitdist(R_conds{1}(2000:end,2),'kernel','Kernel','normal');
y2 = pdf(pd2,x);
pd3 = fitdist(R_conds{2}(2000:end,1),'kernel','Kernel','normal');
y3 = pdf(pd3,x);
pd4 = fitdist(R_conds{2}(2000:end,2),'kernel','Kernel','normal');
y4 = pdf(pd4,x);
h = figure; hold on;
filename = 'Output_I-E&E-Edist';
plot(x,y1,'r--');
plot(x,y2,'b--');
plot(x,y3,'r-');
area(x,y3,'FaceColor','r','FaceAlpha',.5,'EdgeAlpha',.3);
plot(x,y4,'b-');
area(x,y4,'FaceColor','b','FaceAlpha',.5,'EdgeAlpha',.3);
xlim([10, 39]); 
% xlabel('Firing rates (Hz)');
ylabel('Density');
savefigs(h, filename, plotdir, fontsize, [2.9, 1.5]);
%% Property 2. predicted choice behavior
potentiation = linspace(1,4,5); %[1:.5:2.5];
task = 'RT_I-E&E-E';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
sgm = .01;
sgmInput = sgmInput_choice;
w = w0*ones(2);
a = 0*eye(2);
b = b0*eye(2);
thresh = 70;
stoprule = 1;
initialvals = [R0,R0; G0,G0; D0,D0];
filename = sprintf('LDDM_%s_STDPfrom%1.1fto%1.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',...
    task,min(potentiation), max(potentiation), a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
cp = [2 4 8 16 32 64 128 256 512]'/1000;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    [V1, STDP_V] = meshgrid(scale*(1+cp), potentiation);
    [V2, STDP_GR] = meshgrid(scale*(1-cp), potentiation);
    STDP_RG = STDP_V;
    Vinput.V1 = V1;
    Vinput.V2 = V2;
    Vprior.V1 = ones(size(V1))*scale;
    Vprior.V2 = ones(size(V2))*scale;
    [rt, choice, ~] = LDDM_Rndinput_STDP_GPUrv1(Vprior, Vinput, STDP_V, STDP_RG, STDP_GR, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    ACC = gather(mean(2-squeeze(choice),3,'omitnan'));
    meanRT = gather(mean(squeeze(rt),3,'omitnan'));
    save(output,'ACC','meanRT');
else
    load(output);
end
%
mygray = [.9:-.9/numel(potentiation):0];
h = figure;
filename = 'RT_ACC_I-E&E-E';
subplot(1,2,1); hold on;
for level = 1:numel(potentiation)
    plot(cp*100,ACC(level,:)*100,'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz/2,'LineWidth',lwd);
end
set(gca,'XScale','log');
ylabel('% Correct');
ylim([50, 100]);
xticks([1,10,100]);
xticklabels({'1','10','100'});
% lgd = legend(cellstr(string(potentiation)),...
%     'Location','SouthEast', 'FontSize', fontsize-4, 'Box','off');
% title(lgd,'LTP','FontSize',fontsize-4);
savefigs(h, filename, plotdir,fontsize, [4, 1.5]) ;
subplot(1,2,2); hold on;
for level = 1:numel(potentiation)
    plot(cp*100,meanRT(level,:),'.-','Color',mygray(level)*[1,1,1],'MarkerSize',mksz/2,'LineWidth',lwd);
end
set(gca,'XScale','log');
ylabel('RT (s)');
ylim([0, 4]);
xticks([1,10,100]);
xticklabels({'1','10','100'});
savefigs(h, filename, plotdir,fontsize, [4, 1.5]);

%% Property 3. ACC as a function of Input noise, two STDP levels
task = 'RT_I-E&E-E';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
sgm = .01;
w = w0*ones(2);
a = a0*eye(2);
b = b0*eye(2);
sims = 102400;
potentiation = [1, 4];
sgmInputvec = linspace(0,2,100);
filename = sprintf('LDDM_%s_STDP%.1f_%.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3fto%.3f_sims%i',task,min(potentiation),max(potentiation),a0,b0,sgm,min(sgmInputvec),max(sgmInputvec),sims);
output = fullfile(Simdir,[filename, '.mat']);
c = c_choice;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    for pi = 1:numel(potentiation)
        STDP_V = potentiation(pi);
        STDP_a = potentiation(pi);
        STDP_G = potentiation(pi);
        fprintf('potentiation %.3f\n',potentiation(pi));
        [V1, ~] = meshgrid(scale*(1+c), sgmInputvec);
        [V2, sgmInput] = meshgrid(scale*(1-c), sgmInputvec);
        Vinput.V1 = V1;
        Vinput.V2 = V2;
        Vprior.V1 = ones(size(V1))*scale;
        Vprior.V2 = ones(size(V2))*scale;
        [rt, choice, ~] = LDDM_Rndinput_STDP_GPU(Vprior, Vinput, STDP_V, STDP_a, STDP_G, w, a, b,...
            sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
        ACC(pi,:) = gather(mean(2-squeeze(choice),2,'omitnan'));
        meanRT(pi,:) = gather(mean(squeeze(rt),2,'omitnan'));
    end
    save(output,'ACC','meanRT');
else
    load(output);
end
%
h = figure; hold on;
plot(sgmInputvec, ACC(1,:)*100,'--','MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(2)*[1,1,1]);
plot(sgmInputvec, ACC(2,:)*100,'MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(5)*[1,1,1]);
%legend({'Baseline','+LTP'},'Box','off','Location','northeast');
%xlabel('Input noise (a.u.)');
ylabel('% Correct');
ylim([50, 100]);
savefigs(h, ['ACCoversgmInput_', filename], plotdir,fontsize, [2,1.5]);

%% largest amplitude and iSTDP
task = 'VR_I-E&E-E';
potentiation = [1:.5:2.5];
c = 3.2/100;
sgm = 0;
sgmInput = 0;
triggert = Inf;
stimdur = Inf;
dur = 4.5;
initialvals = [1,1;2,2;0,0]*eqlb/3;
Rstore = [];
h = figure;
filename = sprintf('R dynamic over STDP %s',task);
for level = 1:length(potentiation)
    Vinput = [1+c, 1-c]*scale;
    a = a0*eye(2);
    b = 0*eye(2);
    w = w0*ones(2);
    STDP_V = potentiation(level);
    STDP_a = potentiation(level);
    STDP_G = potentiation(level);
    [choice, rt, R, G, I, Vcourse] = LDDM_RndInput_STDP(Vprior, Vinput, STDP_V, STDP_a, STDP_G, w, a, b,...
    sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    Rstore(level) = max(R(:,1));
    subplot(1,4,level);
    plot(R);
    ylim([10,40]);
    xlabel('Time (a.u.)');
    ylabel('Firing rates (Hz)');
end
savefigs(h, [filename], plotdir,fontsize-5, [6,2]);

h = figure; hold on; 
filename = sprintf('MaxFR_%s',task);
plot(potentiation,Rstore,'k-','MarkerSize', mksz,'LineWidth',lwd);
xlabel('STDP');
ylabel('Max firing rates (Hz)');
savefigs(h, [filename], plotdir,fontsize, [2,1.5]);