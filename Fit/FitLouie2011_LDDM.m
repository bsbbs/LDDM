% Fit the mean firing rates from Louie, et al., 2011
% Glgdir = 'G:\My Drive\LDDM\Fit';
Glgdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Fit';
% Glgdir = '/Volumes/GoogleDrive/My Drive/LDDM/Fit/';
addpath(genpath(fullfile(Glgdir, 'bads-master_2019')));
% addpath(genpath('~/Documents/LDDM/utils'));
addpath(genpath('C:\Users\Bo\Documents\LDDM\utils'));
out_dir = fullfile(Glgdir,'Rslts/FitLouie2011/LDDM');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
load(fullfile(Glgdir,'./LouieData/LIPMeanData_TrinaryConditions.mat'));
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
S = 1;
BG = 0;
a = 1.3;
BR = 71;
x = [S,         BG-a,       BR];
LB = [0.5,  -70,         0];
UB = [2,       70,         200];
PLB = [.9,    -1,         60];
PUB = [1.2,     2,          80];
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

opt = optimset('Display', 'on', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-3, 'Algorithm', 'active-set');

% opt = optimoptions('fmincon');
% opt.Display = 'iter';
% opt.StepTolerance = .001;
%% 
if exist(fullfile(out_dir,'FitRslt.mat'),'file')
    load(fullfile(out_dir,'FitRslt.mat'));
end
for i = 1:20
    fprintf('Fit # %i: ',i);
    if size(summaries,1) < i
        x0 = randomX(PLB,PUB);
        [x,fval,~,report] = bads(f,x0,LB,UB,PLB,PUB,options);
%         [x, fval, ~, report] = fmincon(f, x0, [],[],[],[], LB, UB, [], opt); % 
        Rsquared = 1 - fval/var(FR)/(numel(FR)-1);
        summaries(end+1,:)  = [x fval Rsquared];
        reports{i} = report;
    else
        x = summaries(i,1:end-2);
        fval = summaries(i,end-1);
        Rsquared = summaries(i,end);
    end
    
    if ~isnan(fval) && (count == 0 || (fval < RSS_min && sum(x ~= xBest)>0))
        xBest = x;
        S = xBest(1);
        BG_a = xBest(2);
        BR = xBest(3);
        RSS_min = fval;
        count = count + 1;
        sprintf(['Iteration %d: ' repmat('%8.3f\t', 1, length(xBest)) '\n'], i, xBest)
    end
    save(fullfile(out_dir,'FitRslt.mat'),'summaries','xBest', 'S', 'BG_a','BR','RSS_min','reports');
end
%% parameters recovery
% over G0_a
BG_a_rng = [-40, 40];
N = 1000;
output = fullfile(out_dir,'RcvryBG_aRslt.mat');
if exist(output,'file')
    load(output);
else
    summaries = nan(N, 8);
    reports = cell(N,1);
    parfor i = 1:N
        BG_a = rand(1,1)*range(BG_a_rng) +BG_a_rng(1);
        fprintf('Fit BG_a generated = %2.1f: ',BG_a);
        x_gnrt = [xBest(1), BG_a, xBest(3)];
        [~, gnrt_activity] = OLS(x_gnrt,V1,V2,V3,FR);
        f = @(x)OLS(x, V1, V2, V3, gnrt_activity);
        x0 = xBest;
        [x,fval,~,report] = bads(f,x0,LB,UB,PLB,PUB,options);
        Rsquared = 1 - fval/var(gnrt_activity)/(numel(gnrt_activity)-1);
        summaries(i,:)  = [x_gnrt, x, fval, Rsquared];
        reports{i} = report;
    end
    save(output,'summaries','reports');
end

h = figure;
filename = 'PrmtrRcvry_Louie2011';
subplot(1,2,1);
hold on;
plot(BG_a_rng,BG_a_rng,'--','LineWidth',.5);
<<<<<<< HEAD
plot(summaries(:,2), summaries(:,5),'o','MarkerSize',2);
=======
plot(summaries(:,2), summaries(:,5),'o','MarkerSize',5);
>>>>>>> 35d07b0b2577f7c984c16cb486aa12b3ddabb46d
xlabel('Generated B_G-\alpha','FontAngle','italic','FontName','Times');
ylabel('Recovered B_G-\alpha','FontAngle','italic','FontName','Times');
xticks(-40:20:40);
yticks(-40:20:40);
xlim([-41,41]);
ylim([-41,41]);
mysavefig(h, filename, plot_dir, 12, [7, 3]);

% over B
BR_rng = [0, 140];
N = 1000;
output = fullfile(out_dir,'RcvryBR_Rslt.mat');
if exist(output,'file')
    load(output);
else
    summaries = nan(N, 8);
    reports = cell(N,1);
    parfor i = 1:N
        BR = rand(1,1)*range(BR_rng) + min(BR_rng);
        fprintf('Fit B generated = %2.1f: ',BR);
        x_gnrt = [xBest(1), xBest(2), BR];
        [~, gnrt_activity] = OLS(x_gnrt,V1,V2,V3,FR);
        f = @(x)OLS(x, V1, V2, V3, gnrt_activity);
        x0 = xBest;
        [x,fval,~,report] = bads(f,x0,LB,UB,PLB,PUB,options);
        Rsquared = 1 - fval/var(gnrt_activity)/(numel(gnrt_activity)-1);
        summaries(i,:)  = [x_gnrt, x, fval,  Rsquared];
        reports{i} = report;
    end
    save(output,'summaries','reports');
end
subplot(1,2,2);
hold on;
plot([0,140],[0,140],'--','LineWidth',.5);
<<<<<<< HEAD
plot(summaries(:,3), summaries(:,6),'o','MarkerSize',2);
=======
plot(summaries(:,3), summaries(:,6),'o','MarkerSize',5);
>>>>>>> 35d07b0b2577f7c984c16cb486aa12b3ddabb46d
xlabel('Generated B_R','FontAngle','italic','FontName','Times');
ylabel('Recovered B_R','FontAngle','italic','FontName','Times');
xticks(0:40:140);
yticks(0:40:140);
xlim([-1,141]);
ylim([-1,141]);
mysavefig(h, filename, plot_dir, 12, [7, 3]);
%%
if visulize
    h = figure;  hold on;
    filename = sprintf('FitLouie_bads_%i_S%0.3e_BG-a%1.3f_BR%2.1f_fval%1.3e_Rsqd%1.6f',i,S,BG_a,BR,fval,Rsquared);
    [RSS, R1s] = OLS(xBest,V1,V2,V3,FR);
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
    legend(flip(ln),flip({'V_{in} = 0','V_{in} = 50','V_{in} = 100','V_{in} = 200'}),'Box','off','fontsize',7);
    legend([flip(pts) flip(ln)],flip({'V_{in} = 0','V_{in} = 50','V_{in} = 100','V_{in} = 200','','','',''}),'NumColumns', 2, 'Box','off','fontsize',7);
    mysavefig(h, filename, plot_dir, 12, [3, 3]);
    saveas(h,fullfile(plot_dir,[filename, '.fig']));
end
%%
function [RSS, R1s] = OLS(x,V1,V2,V3,FR)
S = x(1);
BG_a = x(2);
BR = x(3);
R1s = [];
for i = 1:numel(V1)
    % fprintf('%i.',i)
    V = [V1(i) V2(i) V3(i)]+BR;
    tmp = SolveEqlb(V,BG_a);
    if numel(tmp) > 1
        warning('Real positive solution more than 1. Recorded only the first value.');
    end
    R1s(i,:) = tmp(1);
end
% fprintf('\n');
R1s = S*R1s/max(R1s)*max(FR);
RSS = sum((R1s - FR).^2);
end
%%
function R1 = SolveEqlb(V,BG_a)
w = 1;
v = 1;
b = 0;
syms R1 R2 R3
if V(1) > 0 && V(2) > 0 && V(3) > 0
    eqns = [(V(1)/R1 - (w - b)*R1 - (1+BG_a) - v*R3)/v == R2, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
        (V(2)/R2 - (w - b)*R2 - (1+BG_a) - v*R1)/v == R3, ... % dR2/dt = 0: R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v
        (V(3)/R3 - (w - b)*R3 - (1+BG_a) - v*R2)/v == R1]; % dR3/dt = 0: R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v;
    vars = [R1 R2 R3];
    [AnswR1, ~, ~] = solve(eqns, vars);
elseif V(1) > 0 && (V(2) <= 0 || V(3) <= 0)
    if V(2) <= 0 && V(3) > 0
        warning('V2 <= 0');
        eqns = [(V(1)/R1 - (w - b)*R1 - (1+BG_a) - v*R3)/v == 0, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
            (V(3)/R3 - (w - b)*R3 - (1+BG_a))/v == R1]; % dR3/dt = 0: R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v;
        vars = [R1 R3];
        [AnswR1, ~] = solve(eqns, vars);
    elseif V(3) <= 0 && V(2) > 0
        warning('V3 <= 0');
        eqns = [(V(1)/R1 - (w - b)*R1 - (1+BG_a))/v == 0, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
            (V(2)/R2 - (w - b)*R2 - (1+BG_a) - v*R1)/v == 0]; % dR2/dt = 0: R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v
        vars = [R1 R2];
        [AnswR1, ~] = solve(eqns, vars);
    elseif V(2) <= 0 && V(3) <= 0
        warning('V2, V3 <= 0');
        eqns = [(V(1)/R1 - (w - b)*R1 - (1+BG_a))/v == 0]; % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
        vars = [R1];
        [AnswR1] = solve(eqns, vars);
    end
elseif V(1) <= 0
    warning('V1 <= 0');
    AnswR1 = 0;
end
PositiveRealAnsw = double(AnswR1) >=0 & double(imag(AnswR1)) == 0;
R1 = [];
for i = [find(PositiveRealAnsw)]'
    R1(end+1) = double(AnswR1(i));
end
end
%%
function X0 = randomX(lowerBd,upperBd)
X0 = lowerBd + (upperBd - lowerBd) .* rand(size(lowerBd));
end