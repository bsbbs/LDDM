%% Reward rate as a function of I-E vs. E-E connection weights ratio
%% Property 3. ACC as a function of Input noise, two STDP levels
% parameters initialization
dt = .001;
predur = 0;
presentt = dt;
w = w0*ones(2);
a0 = 0;
a = a0*eye(2);
b = b0*eye(2);
sgm = .01;
sgmInput = sgmInput_choice;
%%
task = 'RR';
triggert = presentt;
stimdur = Inf;
dur = 4.5;
sims = 102400;
potentiation = .01:.02:1.6;
eqlb = 32;
scale = (2*w0 - b0)*eqlb^2.*potentiation + (1-a0)*eqlb;
R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0).*potentiation))/2/(2*w0 - b0)./potentiation;
% R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0)))/2/(2*w0 - b0);
D0 = b0*R0;
G0 = (2*w0-b0)*R0;
initialvals = [R0(1),R0(1); G0(1),G0(1); D0(1),D0(1)];
filename = sprintf('LDDM_%s_STDP%.1f_%.1f_a%1.2f_b%1.2f_scalemin%3.0fmax%3.0f_sgm%1.1fsinpt%0.3f_sims%i',task,min(potentiation),max(potentiation),a0,b0,min(scale),max(scale),sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
STDP_v = 1;
STDP_a = 1;
c = c_choice;
ACC = [];
meanRT = [];
meanRTc = [];
meanRTw = [];
clear Vinput Vprior;
if ~exist(output,'file')
    STDP_G = potentiation';
    V1 = scale'*(1+c);
    V2 = scale'*(1-c);
    Vinput.V1 = V1;
    Vinput.V2 = V2;
    Vprior.V1 = scale';
    Vprior.V2 = scale';
    [rt, choice, ~] = LDDM_Rndinput_STDP_GPU(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
        sgm, sgmInput*scale', Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims);
    ACC = gather(mean(2-squeeze(choice),2,'omitnan'));
    meanRT = gather(mean(squeeze(rt),2,'omitnan'));
    for pi = 1:numel(potentiation)
        meanRTc(pi,1) = mean((rt(pi,1,choice(pi,1,:) == 1)));
        meanRTw(pi,1) = mean((rt(pi,1,choice(pi,1,:) == 2)));
    end
    save(output,'ACC','meanRT','meanRTc','meanRTw');
else
    load(output);
end

% 
h = figure;
subplot(2,2,1);
hold on;
mask = 1:50;
plot(potentiation(mask), ACC(mask)*100,'k-','MarkerSize',mksz,'LineWidth',lwd);
xlabel('I->E weight');
ylabel('% Correct');
ylim([50, 100]);
savefigs(h, ['ACCoverI-E_', filename], plotdir,fontsize, [2,1.5]);
subplot(2,2,2); hold on;
plot(potentiation(mask), meanRT(mask),'k-','MarkerSize',mksz,'LineWidth',lwd);
xlabel('I->E weight');
ylabel('RT (s)');
savefigs(h, ['ACCoverI-E_', filename], plotdir,fontsize, [2,1.5]);

subplot(2,2,3); hold on;
plot(meanRT(mask),  ACC(mask),'k-','MarkerSize',mksz,'LineWidth',lwd);
% for ri = 1:5
%     RR = ACC./(meanRTc.*ACC + meanRTw.*(1 - ACC) + ri);
%     bestT = max(RR)
%     time = -ri:.001:bestT;
%     plot(time,)
%     
% end
xlim([-.5,max(meanRT(mask))]);
ylim([0, 1]);
xlabel('RT (s)');
ylabel('Reward in total');
savefigs(h, ['ACCoverI-E_', filename], plotdir,fontsize, [2,1.5]);

subplot(2,2,4); hold on;
for ri = 0:5
    RR = ACC./(meanRT + ri); %(meanRTc.*ACC + meanRTw.*(1 - ACC) + ri);
    plot(potentiation(mask), RR(mask),'-','MarkerSize',mksz,'LineWidth',lwd,'Color',mygray(6-ri,:));
end
xlabel('I->E weight');
ylabel('Reward rate');
savefigs(h, ['ACCoverI-E_', filename], plotdir,fontsize, [5,5]);
%% Activities dynamics
% Reaction-time task
stimdur = Inf;
triggert = presentt;
dur = 1.8; % second
potentiation = 0:.02:1.6;
c = c_choice;
% scale = 20;
a0 = 0;
STDP_v = 1;
STDP_a = 1;
figure; hold on;
rtvec = [];
for pi = 1:numel(potentiation)
    STDP_G = potentiation(pi);
    eqlb = 32; % control the starting value as fixed
    scale = (2*w0 - b0)*STDP_G*eqlb^2 + (1-a0)*eqlb; % change scaling to fix starting value

    Vinput = [1+c, 1-c]*scale;

    R0 = ((a0-1)+sqrt((1-a0)^2 + 4*scale*(2*w0 - b0)*STDP_G))/2/(2*w0 - b0)/STDP_G; % the starting value as fixed
    D0 = b0*R0;
    G0 = (2*w0-b0)*R0;
    initialvals = [R0,R0; G0,G0; D0,D0];

    rng(5);
    [choice, rt, R, G, D, Vcourse] = LDDM_RndInput_STDP(Vprior, Vinput, STDP_v, STDP_a, STDP_G, w, a, b,...
        sgm, sgmInput*scale, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule);
    rtvec(pi) = rt;
    x = (1:numel(R(:,1)))*dt;
    plot(x, R(:,1), '-', 'LineWidth', lwd, 'Color', myred(6,:));
    plot(x, R(:,2), '-', 'LineWidth', lwd, 'Color', myblue(6,:));
    % plot(x, G(:,1), 'r--', 'LineWidth', lwd);
    % plot(x, G(:,2), 'b--', 'LineWidth', lwd);
    % plot(x, D(:,1), 'r-.', 'LineWidth', lwd);
    % plot(x, D(:,2), 'b-.', 'LineWidth', lwd);
end
plot(rtvec, thresh*ones(size(rtvec)),'k-','LineWidth',.5);
text(mean(rtvec), thresh*1.05,'Threshold','FontName','Arial');
xlim([-50*dt, dur]);
ylim([10,thresh]);
ylabel('Firing rates (Hz)');
xlabel('Time (s)');
savefigs(h, ['R_DynamicsI-E_', filename], plotdir, fontsize, [2.9,4]);