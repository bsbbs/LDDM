% Fit the mean firing rates from Louie, et al., 2011
cd('G:\My Drive\LDDM\Fit');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Fit');
addpath(genpath('./bads-master'));
% addpath(genpath('~/Documents/LDDM/'));
addpath(genpath('C:\Users\Bo\Documents\LDDM'));
out_dir = './Rslts/FitLouie2011/TrinaryRNM_TimeWindow1s';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
load('./Data/LIPMeanData_TrinaryConditions.mat');
mksz = 18;
lwd = 2.0;
%%
Vi = LIPMeanData_TrinaryConditions.Vi;
Vo = LIPMeanData_TrinaryConditions.Vo;
V1 = LIPMeanData_TrinaryConditions.V1;
V2 = LIPMeanData_TrinaryConditions.V2;
V3 = LIPMeanData_TrinaryConditions.V3;
FR = LIPMeanData_TrinaryConditions.FR;
Vilist = unique(Vi);
V1list = unique(V1);
%% Values in equilibirum
miu0 = 30;
sgm = 0;
I0 = .3255;
JNp = .2609;
JNn = .0497;
S = 1;
gamma = .64;
tauS = .1;
tauAMPA = .002;
dur = 4;
dt = .001;
presentt = 0;
stimdur = dur;
thresh = Inf;
initialvals = [2 ,2, 2; .1, .1, .1];
stoprule = 0;
x = [JNp,JNn,I0,S];
LB = [0, 0, 0, .6];
UB = [1, 1, 1, 1.7];
PLB = [.2, .1, .2, .9];
PUB = [.4, .2, .4, 1.1];
tic;
[RSS, R1s] = OLS(x,V1,V2,V3,FR);
toc
f = @(x)OLS(x, V1, V2, V3, FR);

count = 0;
summaries = [];
reports = [];
visulize = 1;
options = bads('defaults');     % Default options
options.Display = 'iter';
options.UncertaintyHandling = false;    % Function is deterministic
%% 
if exist(fullfile(out_dir,'FitRslt.mat'),'file')
    load(fullfile(out_dir,'FitRslt.mat'));
end
for i = 1:20
    fprintf('Fit # %i: ',i);
    if size(summaries,1) < i
        x0 = randomX(LB,UB);
        [x,fval,~,report] = bads(f,x0,LB,UB,PLB,PUB,options);
        
        % [X, tMin] = fmincon(f, X0, [],[],[],[], lowerBd, upperBd, [], opt);
        Rsquared = 1 - fval/var(FR)/(numel(FR)-1);
        summaries(end+1,:)  = [x fval Rsquared];
        reports{i} = report;
    else
        x = summaries(i,1:4);
        fval = summaries(i,5);
        Rsquared = 1 - fval/var(FR)/(numel(FR)-1);
        % Rsquared = summaries(i,end);
        % summaries(end+1,:)  = [x fval Rsquared];
    end
    %%
    if visulize
        h = figure(i);
        % subplot(1,2,1); 
        hold on;
        filename = sprintf('FitLouie_bads_%i_JNp%1.4f_JNn%1.4f_I0%1.4f_S%1.4f_fval%1.3f_Rsqd%1.4f',i,x,fval,Rsquared);
        [RSS, R1s, NU] = OLS(x,V1,V2,V3,FR);
        % Rsquared = 1 - sum((FR - R1s).^2)/(var(FR)*(numel(FR)-1))
        for ii = 1:numel(V1list)
            mask = V1 == V1list(ii);
            Vo = V2 + V3;
            pts(ii) = plot(Vo(mask),FR(mask),'.','MarkerSize',mksz-6);
            [ord, I] = sort(Vo(mask));
            R1 = R1s(mask);
            ln(ii) = plot(ord,R1(I),'-','Color',pts(ii).Color,'LineWidth',lwd);
        end
        ylabel('Rescaled activity');
        xlabel('V_{out}');
        legend([flip(pts) flip(ln)],flip({'V_{in} = 0','V_{in} = 50','V_{in} = 100','V_{in} = 200','','','',''}),'NumColumns', 2, 'Box','off','fontsize',7);
        savefigs(h, filename, plot_dir, 12, [3, 3]);
        saveas(h,fullfile(plot_dir,[filename, '.fig']));
        % the actual activity
        h = figure(i+40); hold on;
        for i = 1:numel(V1)
            plot(NU{i}(:,1),'-');
        end
        ylabel('FR (Hz)');
        xlabel('time (msecs)');
        savefigs(h, ['FR_' filename], plot_dir, 12, [2.5, 2.5]);
        saveas(h,fullfile(plot_dir,['FR_', filename, '.fig']));
    end
    %%
    fprintf('\n');
    pause(.3);
    if ~isnan(fval) && (count == 0 || (fval < RSS_min && sum(x ~= xBest)>0))
        xBest = x;
        RSS_min = fval;
        count = count + 1;
        %              sprintf(['Subject %d\t: Round-%d ' repmat('%8.3f\t', 1, length(XBest)) '\n'], jj, i, XBest)
    end
    save(fullfile(out_dir,'FitRslt.mat'),'summaries','xBest','RSS_min','reports');
end

%%

%%
function [RSS, R1s, NU] = OLS(x,V1,V2,V3,FR)
miu0 = 30;
sgm = 0;
JN = [x(1) -x(2) -x(2)
    -x(2)   x(1) -x(2)
    -x(2) -x(2) x(1)]; % nA
I0 = x(3);
S = x(4);
gamma = .64;
tauS = .1;
tauAMPA = .002;
dur = 4;
dt = .001;
presentt = 0;
stimdur = dur;
thresh = Inf;
initialvals = [2 ,2, 2; .1, .1, .1];
stoprule = 0;
R1s = [];
for i = 1:numel(V1)
    % fprintf('%i.',i)
    V = [V1(i) V2(i) V3(i)];
    cp = V/max(V1)*2;
    [~, ~, nu_wind, ~, ~, ~] = wong06(cp,miu0,sgm,I0,JN,...
        gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule);
    R1s(i,:) = mean(nu_wind(1:1000,1));
    NU{i} = nu_wind;
end
% fprintf('\n');
R1s = S*R1s/max(R1s)*max(FR);
RSS = sum((R1s - FR).^2);
end
%%
function X0 = randomX(lowerBd,upperBd)
X0 = lowerBd + (upperBd - lowerBd) .* rand(size(lowerBd));
end