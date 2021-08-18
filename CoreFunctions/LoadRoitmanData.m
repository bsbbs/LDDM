function data = LoadRoitmanData(datadir)
addpath(datadir);
ColumnNames608
load T1RT.mat;
x(:,R_RT) = x(:,R_RT)/1000;
cohlist = unique(x(:,R_COH));
quantilemat = [];
proportionmat = [];
Expected = [];
numbins = 30;
histmat = [];
rtrange = [];
histbinsmat = [];
bincentermat = [];
Valuesmat = [];
% define parameters for method QMLE, 
% reference: Heathcote & Australia, and Mewhort, 2002.
qntls = .1:.1:.9; %[.1,.3,.5,.7,.9]; % quantile probability
expd_qntls = [0, qntls, 1];
P = expd_qntls(2:end) - expd_qntls(1:end-1); % proportion in each bin
h = figure;
for vi = 1:length(cohlist)
    Lcoh = x(:,R_COH)==cohlist(vi);
    RT_all = x(Lcoh,R_RT);
    On(vi) = numel(RT_all); % total number of observation
    rtmax = max(RT_all);
    rtmin = min(RT_all);
    histbins = linspace(rtmin,rtmax,numbins+1);
    if vi == 1
        quantilemat(1:5,vi) = quantile(RT_all,[.1,.3,.5,.7,.9]);
        quantilemat(6:10,vi) = quantile(RT_all,[.1,.3,.5,.7,.9]);
        proportionmat(vi) = .5;
        total_trials = length(RT_all);
        histg_all = histogram(RT_all, 'BinEdges', histbins, 'Visible',0);
        tmp1 = histg_all.Values/total_trials/2;
        tmp2 = tmp1;
        Values1 = histg_all.Values/2;
        Values2 = Values1;
        ON(:,1,vi) = numel(RT_all)*P/2; % numbers of observation in each bin for correct trial
        ON(:,2,vi) = ON(:,1,vi); % numbers of observation in each bin for error trial
        q(:,1,vi) = quantile(RT_all,qntls); % RT value on quantiles, correct trial
        q(:,2,vi) = q(:,1,vi); % RT value on quantiles, error trial
        OP(:,1,vi) = P/2; % the distributed proportion in each bin, correct trial
        OP(:,2,vi) = OP(:,1,vi); % the distributed proportion in each bin, error trial
    else
        Corr = x(:,R_DIR) == x(:,R_TRG);
        Wro = x(:,R_DIR) ~= x(:,R_TRG);
        RT_corr = x(Lcoh & Corr,R_RT);
        RT_wro = x(Lcoh & Wro, R_RT);
        quantilemat(1:5,vi) = quantile(RT_corr,[.1,.3,.5,.7,.9]);
        quantilemat(6:10,vi) = quantile(RT_wro,[.1,.3,.5,.7,.9]);
        total_trials = (length(RT_corr)+length(RT_wro));
        proportionmat(vi) = length(RT_corr)/total_trials;
        histg_corr = histogram(RT_corr, 'BinEdges', histbins, 'Visible',0);
        tmp1 = histg_corr.Values/total_trials;
        Values1 = histg_corr.Values;
        histg_wro = histogram(RT_wro, 'BinEdges', histbins, 'Visible',0);
        tmp2 = histg_wro.Values/total_trials;
        Values2 = histg_wro.Values;
        ON(:,1,vi) = numel(RT_corr)*P; % numbers of observation in each bin for correct trial
        ON(:,2,vi) = numel(RT_wro)*P; % numbers of observation in each bin for error trial
        q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
        q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
        OP(:,1,vi) = numel(RT_corr)/On(vi)*P; % the distributed proportion in each bin, correct trial
        OP(:,2,vi) = numel(RT_wro)/On(vi)*P; % the distributed proportion in each bin, error trial
    end
    Expected(:,vi) = [[.1,.2,.2,.2,.2,.1]'.*proportionmat(vi); [.1,.2,.2,.2,.2,.1]'.*(1-proportionmat(vi))];
    Valuesmat = [Valuesmat; [Values1 Values2]];
    histmat = [histmat; [tmp1 tmp2]];
    rtrange = [rtrange; [rtmin rtmax]]; 
    histbinsmat = [histbinsmat; [histbins]];
    bincenter = (histbins(1:end-1) + histbins(2:end))/2;
    bincentermat = [bincentermat; [bincenter]];
end
close(h);
quantilemat(isnan(quantilemat)) = 0;
data.quantilemat = quantilemat;
data.proportionmat = proportionmat;
data.Expected = Expected;
data.numbins = numbins;
data.histmat = histmat;
data.Valuesmat = Valuesmat;
data.rtrange = rtrange;
data.histbins = histbinsmat;
data.bincenter = bincentermat;
data.q = q;
data.On = On;
data.ON = ON;
data.OP = OP;
data.qntls = qntls;