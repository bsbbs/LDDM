% Fit the mean firing rates from Louie, et al., 2011
cd('/Volumes/GoogleDrive/My Drive/LDDM/Fit');
addpath(genpath('./bads-master'));
out_dir = './Rslts/FitLouie2011';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
plot_dir = fullfile(out_dir,'graphics');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end
load('./Data/LIPMeanData.mat');
mksz = 18;
lwd = 2.0;
%%
Vi = LIPMeanData.Vi;
Vo = LIPMeanData.Vo;
FR = LIPMeanData.FR;
Vilist = unique(Vi);

%% Values in equilibirum
a = 50;
GAMMA = 1;
I0 = 0.1;
x = [a,GAMMA,I0];
LB = [0, 1, 0.1];
UB = [80, 500, 500];
PLB = [5, 100, 80];
PUB = [40, 300, 120];
tic;
[RSS, R1s] = OLS(x,Vi,Vo,FR);
toc
f = @(x)OLS(x, Vi, Vo, FR);

% opt = optimset('Display', 'on', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-3, 'Algorithm', 'active-set');
count = 0;
summaries = [];
visulize = 1;
options = bads('defaults');     % Default options
options.Display = 'iter';
options.UncertaintyHandling = false;    % Function is deterministic
% options.UncertaintyHandling = true;    % Function is deterministic
for i = 1:20
    fprintf('Fit # %i: ',i);
    x0 = randomX(LB,UB);
    [x,fval,~,output] = bads(f,x0,LB,UB,PLB,PUB,options);
    
    % [X, tMin] = fmincon(f, X0, [],[],[],[], lowerBd, upperBd, [], opt);
    summaries(end+1,:)  = [x fval];
    %%
    if visulize
        h = figure(i);  hold on;
        filename = sprintf('FitLouie_bads_%i_a%1.2f_GAMMA%3.1f_I0%3.1f_fval%1.3f',i,x,fval);
        [RSS, R1s] = OLS(x,Vi,Vo,FR);
        for ii = 1:numel(Vilist)
            mask = Vi == Vilist(ii);
            ln = plot(Vo(mask),FR(mask),'.','MarkerSize',mksz-6);
            [ord, I] = sort(Vo(mask));
            R1 = R1s(mask);
            plot(ord,R1(I),'-','Color',ln.Color,'LineWidth',lwd);
        end
        ylabel('rescaled FR');
        xlabel('V_{out}');
        savefigs(h, filename, plot_dir, 12, [2.5, 2.5]);
    end
    %%
    fprintf('\n');
    pause(.3);
    if ~isnan(tMin) && (count == 0 || (fval < ll_min && sum(x ~= xBest)>0))
        xBest = x;
        ll_min = fval;
        count = count + 1;
        %              sprintf(['Subject %d\t: Round-%d ' repmat('%8.3f\t', 1, length(XBest)) '\n'], jj, i, XBest)
    end
end

%%

%%
function [RSS, R1s] = OLS(x,Vi,Vo,FR)
a = x(1);
GAMMA = x(2);
I0 = x(3);
if numel(Vi) ~= numel(Vo)
    error('Input Vi dimension does not match Vo');
end
R1s = [];
for i = 1:numel(Vi)
    % fprintf('%i.',i)
    V = [Vi(i) Vo(i)]+I0;
    tmp = SolveEqlb(a,GAMMA,V);
    if numel(tmp) > 1
        warning('Real positive solution more than 1. Recorded only the first value.');
    end
    R1s(i,:) = tmp(1);
end
% fprintf('\n');
R1s = max(FR)*R1s/max(R1s);
RSS = sum((R1s - FR).^2);
end
%%
function R1 = SolveEqlb(a,GAMMA,V)
w = 1;
v = 1;
b = 0;
syms R1 R2
if V(1) > 0 && V(2) > 0
    eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a))/v == R2, ... % dR1/dt = 0
        (V(2)/R2 - (w - b)*R2 - (GAMMA-a))/v == R1];% dR2/dt = 0
    vars = [R1 R2];
    [AnswR1, ~] = solve(eqns, vars);
elseif V(1) > 0 && V(2) <= 0
    eqns = [(V(1)/R1 - (w - b)*R1 - (GAMMA-a))/v == 0]; % dR1/dt = 0
    vars = [R1];
    AnswR1 = solve(eqns, vars);
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