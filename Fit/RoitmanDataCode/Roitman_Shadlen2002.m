%% Fitting the data from Roitman & Shadlen, 2002

%% read in easier format for the RT task
clear all;
close all;
directory = '/Users/bs3667/Documents/Projects/RecurrentModel/Analysis/RoitmanDataCode';
cd(directory);
load('T1RT.mat');
ColumnNames608
%% speed and accuracy
close all;
cohlist = unique(x(:,R_COH));
maxrt = max(x(:,R_RT));
minrt = min(x(:,R_RT));
segrt = maxrt - minrt;
bins = 30;
BinEdge = [minrt:segrt/bins:maxrt];
bank1 = [];
bank2 = [];
for i = 1:length(cohlist)
    Lcoh = x(:,R_COH)==cohlist(i);
    if i == 1
        Dir1 = x(:,R_TRG) == 1;
        Dir2 = x(:,R_TRG) == 2;
        RT.corr = x(Lcoh & Dir1,R_RT);
        RT.wro = x(Lcoh & Dir2, R_RT);
    else
        Corr = x(:,R_DIR) == x(:,R_TRG);
        Wro = x(:,R_DIR) ~= x(:,R_TRG);
        RT.corr = x(Lcoh & Corr,R_RT);
        RT.wro = x(Lcoh & Wro, R_RT);
    end
    hg = histogram(RT.corr,'BinEdges',BinEdge);
    bank1(i,:) = hg.Values;
    if ~isempty(RT.wro)
        hg = histogram(RT.wro,'BinEdges',BinEdge);
        bank2(i,:) = hg.Values;
    else
        bank2(i,:) = zeros(1,bins);
    end
end
BinMiddle = hg.BinEdges(1:end-1) + hg.BinWidth/2;

h=figure;
for ii = 1:length(cohlist)
    subplot(2,3,ii);
    bar(BinMiddle,bank1(ii,:));
    hold on;
    bar(BinMiddle,-bank2(ii,:));
    ylim([-120,230]);
    xlabel('time bin/ms');
    ylabel('frequency');
    title(sprintf('coherence %i',cohlist(ii)));
end
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 8 5];
saveas(h,sprintf('./RoitmanDataCode/RTDistrb_Coherence'),'epsc2');

%% quantile probability plot
h = figure;
hold on;
Rprop = [];
Rpropwro = [];
Rpropcorr = []
ywro = [];
ycorr = [];
for i = 1:length(cohlist)
    Lcoh = x(:,R_COH)==cohlist(i);
    if i == 1
        Dir1 = x(:,R_TRG) == 2;
        Dir2 = x(:,R_TRG) == 1;
        RT.corr = x(Lcoh & Dir1,R_RT);
        RT.wro = x(Lcoh & Dir2, R_RT);
    else
        Corr = x(:,R_DIR) == x(:,R_TRG);
        Wro = x(:,R_DIR) ~= x(:,R_TRG);
        RT.corr = x(Lcoh & Corr,R_RT);
        RT.wro = x(Lcoh & Wro, R_RT);
    end
    if length(RT.wro) > 0
        ywro(:,i) = quantile(RT.wro,[.1 .3 .5 .7 .9]);
    else
        ywro(:,i) = ones(1,5)*nan;
    end
    ycorr(:,i) = quantile(RT.corr,[.1 .3 .5 .7 .9]);
    Rprop = [length(RT.wro)/(length(RT.corr)+length(RT.wro)), Rprop, length(RT.corr)/(length(RT.corr)+length(RT.wro))];
    %Rpropwro = [Rpropwro, length(RT.wro)/(length(RT.corr)+length(RT.wro))];
    %Rpropcorr = [Rpropcorr, length(RT.corr)/(length(RT.corr)+length(RT.wro))];
    plot([ones(1,5)*length(RT.wro)/(length(RT.corr)+length(RT.wro)) ones(1,5)*length(RT.corr)/(length(RT.corr)+length(RT.wro))], [ywro(:,i)' ycorr(:,i)'],'.k', 'MarkerSize', 20);
end

%Rprop = [Rpropwro, Rpropcorr];
for qi = 1:5
    mask = [2:6, 7:12];
    plot(Rprop(mask),[ywro(qi,5:-1:1),ycorr(qi,:)],'-k');
end
ylabel('RT quantile (ms)');
xlabel('Response proportion');
%% save into txt files
dlmwrite('/Users/bs3667/Documents/Glimcher Lab/My work/RecurrentModel/ModelFitting/x.txt', x, 'delimiter','\t');
