% Fit the mean firing rates from Louie, et al., 2011
cd('G:\My Drive\LDDM\Fit');
% cd('/Volumes/GoogleDrive/My Drive/LDDM/Fit');
addpath(genpath('./bads-master'));
addpath(genpath('../utils'));
out_dir = './Rslts/FitLouie2011/TrinaryI0GAMMA';
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
a = 50;
GAMMA = 1;
I0 = 0.1;
S = 1;
x = [a,GAMMA,I0,S];
LB = [0, 1, 0.1, .6];
UB = [80, 500, 500, 1.7];
PLB = [5, 100, 80, .9];
PUB = [40, 300, 120, 1.1];
tic;
[RSS, R1s] = OLS(x,V1,V2,V3,FR);
toc
f = @(x)OLS(x, V1, V2, V3, FR);

% opt = optimset('Display', 'on', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-3, 'Algorithm', 'active-set');
count = 0;
summaries = [];
reports = [];
visulize = 1;
options = bads('defaults');     % Default options
options.Display = 'iter';
options.UncertaintyHandling = false;    % Function is deterministic
% options.UncertaintyHandling = true;    % Function is deterministic
for i = 10:20
    fprintf('Fit # %i: ',i);
    if size(summaries,1) < i
        x0 = randomX(LB,UB);
        [x,fval,~,report] = bads(f,x0,LB,UB,PLB,PUB,options);
        
        % [X, tMin] = fmincon(f, X0, [],[],[],[], lowerBd, upperBd, [], opt);
        summaries(end+1,:)  = [x fval];
        reports{i} = report;
    else
        x = summaries(i,1:end-1);
        fval = summaries(i,end);
    end
    %%
    if visulize
        h = figure(i);  hold on;
        filename = sprintf('FitLouie_bads_%i_a%1.2f_GAMMA%3.1f_I0%3.1f_S%1.2f_fval%1.3f',i,x,fval);
        [RSS, R1s] = OLS(x,V1,V2,V3,FR);
        for ii = 1:numel(V1list)
            mask = V1 == V1list(ii);
            Vo = V2 + V3;
            pts(ii) = plot(Vo(mask),FR(mask),'.','MarkerSize',mksz-6);
            [ord, I] = sort(Vo(mask));
            R1 = R1s(mask);
            ln(ii) = plot(ord,R1(I),'-','Color',pts(ii).Color,'LineWidth',lwd);
        end
        ylabel('rescaled FR');
        xlabel('V_{out}');
        legend(flip(ln),flip({'V_{in} = 0','V_{in} = 50','V_{in} = 100','V_{in} = 200'}),'Box','off','fontsize',7);
        savefigs(h, filename, plot_dir, 12, [2.5, 2.5]);
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
end
save(fullfile(out_dir,'FitRslt.mat'),'summaries','xBest','RSS_min','reports');
%%

%%
function [RSS, R1s] = OLS(x,V1,V2,V3,FR)
a = x(1);
GAMMA = x(2);
I0 = x(3);
S = x(4);
R1s = [];
for i = 1:numel(V1)
    % fprintf('%i.',i)
    V = [V1(i) V2(i) V3(i)]+I0;
    tmp = SolveEqlb(a,GAMMA,V);
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
function R1 = SolveEqlb(a,GAMMA,V)
w = 1;
v = 1;
b = 0;
syms R1 R2 R3
if V(1) > 0 && V(2) > 0 && V(3) > 0
    eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a) - v*R3)/v == R2, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
        (V(2)/R2 - (w - b)*R2 - (GAMMA-a) - v*R1)/v == R3, ... % dR2/dt = 0: R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v
        (V(3)/R3 - (w - b)*R3 - (GAMMA-a) - v*R2)/v == R1]; % dR3/dt = 0: R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v;
    vars = [R1 R2 R3];
    [AnswR1, ~, ~] = solve(eqns, vars);
elseif V(1) > 0 && (V(2) <= 0 || V(3) <= 0)
    if V(2) <= 0 && V(3) > 0
        eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a) - v*R3)/v == 0, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
            (V(3)/R3 - (w - b)*R3 - (GAMMA-a))/v == R1]; % dR3/dt = 0: R1 = (V(3)./R3 - (w - b)*R3 - (1-a) - v*R2)/v;
        vars = [R1 R3];
        [AnswR1, ~] = solve(eqns, vars);
    elseif V(3) <= 0 && V(2) > 0
        eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a))/v == 0, ... % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
            (V(2)/R2 - (w - b)*R2 - (GAMMA-a) - v*R1)/v == 0]; % dR2/dt = 0: R3 = (V(2)./R2 - (w - b)*R2 - (1-a) - v*R1)/v
        vars = [R1 R2];
        [AnswR1, ~] = solve(eqns, vars);
    elseif V(2) <= 0 && V(3) <= 0
        eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a))/v == 0]; % dR1/dt = 0: R2 = (V(1)./R1 - (w - b)*R1 - (1-a) - v*R3)/v
        vars = [R1];
        [AnswR1] = solve(eqns, vars);
    end
elseif V(1) <= 0
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