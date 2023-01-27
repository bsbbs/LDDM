function data = Load_SimData(rtmat, choicemat)
[r, c] = size(rtmat);
if r < c
    rtmat = rtmat';
    choicemat = choicemat';
end
[r, c] = size(choicemat);
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
for vi = 1:c
    RT_all = rtmat(:,vi);
    On(vi) = sum(choicemat(:,vi) == 1 | choicemat(:,vi) == 2); % total number of observation, omit nan
    rtmax = max(RT_all);
    rtmin = min(RT_all);
    histbins = linspace(rtmin,rtmax,numbins+1);
    if vi == 1
        tmp = quantile(RT_all,[.1,.3,.5,.7,.9]);
        quantilemat(1:5,vi) = tmp;
        quantilemat(6:10,vi) = tmp;
        proportionmat(vi) = .5;
        histg_all = histogram(RT_all, 'BinEdges', histbins, 'Visible',0);
        tmp1 = histg_all.Values/On(vi)/2;
        tmp2 = tmp1;
        Values1 = histg_all.Values/2;
        Values2 = Values1;
        ON(:,1,vi) = On(vi)/2*P; % numbers of observation in each bin for correct trial
        ON(:,2,vi) = On(vi)/2*P; % numbers of observation in each bin for error trial
        OP(:,1,vi) = P/2; % the distributed proportion in each bin, correct trial
        OP(:,2,vi) = P/2; % the distributed proportion in each bin, error trial
        q(:,1,vi) = quantile(RT_all,qntls); % RT value on quantiles, correct trial
        q(:,2,vi) = q(:,1,vi); % RT value on quantiles, error trial
        
    else
        Corr = choicemat(:,vi) == 1;
        Wro = choicemat(:,vi) == 2;
        RT_corr = rtmat(Corr,vi);
        RT_wro = rtmat(Wro,vi);
        quantilemat(1:5,vi) = quantile(RT_corr,[.1,.3,.5,.7,.9]);
        quantilemat(6:10,vi) = quantile(RT_wro,[.1,.3,.5,.7,.9]);
        proportionmat(vi) = numel(RT_corr)/On(vi);
        histg_corr = histogram(RT_corr, 'BinEdges', histbins, 'Visible',0);
        tmp1 = histg_corr.Values/On(vi);
        Values1 = histg_corr.Values;
        histg_wro = histogram(RT_wro, 'BinEdges', histbins, 'Visible',0);
        tmp2 = histg_wro.Values/On(vi);
        Values2 = histg_wro.Values;
        ON(:,1,vi) = numel(RT_corr)*P; % numbers of observation in each bin for correct trial
        ON(:,2,vi) = numel(RT_wro)*P; % numbers of observation in each bin for error trial
        OP(:,1,vi) = proportionmat(vi)*P; % the distributed proportion in each bin, correct trial
        OP(:,2,vi) = (1-proportionmat(vi))*P; % the distributed proportion in each bin, error trial
        q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
        q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
        
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
end