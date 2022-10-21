%% Reward rate as a function of I-E vs. E-E connection weights ratio
%% Property 3. ACC as a function of Input noise, two STDP levels
task = 'RR';
predur = 0;
presentt = dt;
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sgm = .01;
sgmInput = sgmInput_choice;
w = w0*ones(2);
a = a0*eye(2);
b = b0*eye(2);
sims = 102400;
potentiation = 0:.1:4;

filename = sprintf('LDDM_%s_STDP%.1f_%.1f_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',task,min(potentiation),max(potentiation),a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
STDP_v = 1;
STDP_a = 1;
c = ones(size(potentiation))*c_choice;
ACC = [];
meanRT = [];
clear Vinput Vprior;
if ~exist(output,'file')
    
    STDP_G = potentiation;
    
    [V1, ~] = meshgrid(scale*(1+c), sgmInput);
    [V2, sgmInput] = meshgrid(scale*(1-c), sgmInput);
    Vinput.V1 = V1;
    Vinput.V2 = V2;
    Vprior.V1 = ones(size(V1))*scale;
    Vprior.V2 = ones(size(V2))*scale;
    [rt, choice, ~] = LDDM_Rndinput_STDP_GPU(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    ACC = gather(mean(2-squeeze(choice),2,'omitnan'));
    meanRT = gather(mean(squeeze(rt),2,'omitnan'));
    
    save(output,'ACC','meanRT');
else
    load(output);
end
%
h = figure; hold on;
plot(potentiation, ACC*100,'--','MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(2)*[1,1,1]);
xlabel('I-E/E-E ratio');
ylabel('% Correct');
ylim([50, 100]);
savefigs(h, ['ACCoverI-E_', filename], plotdir,fontsize, [2,1.5]);




RR = ACC./meanRT;
h = figure;
heatmap(RR)