%% Speed-accuracy tradeoff buddget line
c = [.032, .064, .128, .256, .512]';
c1 = [1 - flip(c); 1; 1 + c];
c2 = ones(size(c1));
c3 = 0:.5:2;
filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_sgm%2.1f_sim%i',task,length(c1),length(unique(c2)),length(unique(c3)),sgm,sims);
simrslt = fullfile(outdir,[filename, '.mat']);
load(simrslt);
mycol = colormap(jet(41));
h = figure; hold on;
for i = 1:length(c3)
    chset1 = squeeze(choice(i,10,:) == 1);
    chset2 = squeeze(choice(i,10,:) == 2);
    qntls = quantile(rt(i,10,chset1 | chset2),9);
    qntls = [min(rt(i,10,chset1 | chset2)) qntls max(rt(i,10,chset1 | chset2))];
    meanrt = [];
    acc = [];
    for q = 1:(length(qntls)-1)
        meanrt(q) = mean(rt(i,10,(chset1 | chset2) & squeeze(rt(i, 10, :) <= qntls(q+1) & rt(i, 10, :) >= qntls(q))));
        acc(q) = mean(2-choice(i,10,(chset1 | chset2) & squeeze(rt(i, 10, :) <= qntls(q+1) & rt(i, 10, :) >= qntls(q))));
    end
    plot(meanrt, acc, '.', 'MarkerSize',mksz*3, 'color',mycol(1+(i-1)*10,:));
end
h = figure; hold on;
for i = 1:length(c3)
    meanrt = [];
    acc = [];
    for j = 1:length(c1)
        chset1 = squeeze(choice(i,j,:) == 1);
        chset2 = squeeze(choice(i,j,:) == 2);
        meanrt(j) = mean(rt(i,j,(chset1 | chset2) ));
        acc(j) = mean(2-choice(i,j,(chset1 | chset2)));
    end
    plot(meanrt, acc, '.', 'MarkerSize',mksz*3, 'color',mycol(1+(i-1)*10,:));
end
c = [.256]';
c1 = [1 + c];
c2 = ones(size(c1));
c3 = 0:.04:2;
[cp.cp1, ~] = meshgrid(c1, c3);
[cp.cp2, cp.cp3] = meshgrid(c2, c3);
filename = sprintf('LDDM_%s_%ic1_%ic2_%ic3_sgm%2.1f_sim%i',task,length(c1),length(unique(c2)),length(unique(c3)),sgm,sims);
simrslt = fullfile(outdir,[filename, '.mat']);
load(simrslt);
mycol = colormap(jet(length(c3)));
for i = 1:length(c3)
    meanrt = [];
    acc = [];
    for j = 1
        chset1 = squeeze(choice(i,j,:) == 1);
        chset2 = squeeze(choice(i,j,:) == 2);
        meanrt(j) = mean(rt(i,j,(chset1 | chset2) ));
        acc(j) = mean(2-choice(i,j,(chset1 | chset2)));
    end
    plot(meanrt, acc, '.', 'MarkerSize',mksz*3, 'color',mycol(i,:));
end
xlabel('RT (s)');
ylabel({'Conditional choice probability';'Opt. 1 vs. Opt. 2'});
%% LDDM, 7 params fitting results
dat_dir = '/gpfs/data/glimcherlab/BoShen/RecurrentModel/FitDynamic/Rslts/FitBhvr7ParamsIV_QMLE_SvrGPU/graphics';
params = [	0	1.433631	25.35945	3251.289056	0.185325	0.224459	0.323132	16539.138186];
name = sprintf('a%2.2f_b%1.2f_sgm%2.1f_scale%4.1f_tau%1.2f_%1.2f_%1.2f_nLL%4.0f',params);
load(fullfile(dat_dir,sprintf('PlotData_%s.mat',name)));
h = figure; hold on;
for i = 2:6
    chset1 = squeeze(choicemat(:,i) == 1);
    chset2 = squeeze(choicemat(:,i) == 2);
    qntls = quantile(rtmat(chset1 | chset2,i),9);
    qntls = [min(rtmat(chset1 | chset2,i)) qntls max(rtmat(chset1 | chset2,i))];
    meanrt = [];
    acc = [];
    for q = 1:(length(qntls)-1)
        meanrt(q) = mean(rtmat((chset1 | chset2) & squeeze(rtmat(:,i) <= qntls(q+1) & rtmat(:,i) >= qntls(q)), i));
        acc(q) = mean(2-choicemat((chset1 | chset2) & squeeze(rtmat(:,i) <= qntls(q+1) & rtmat(:,i) >= qntls(q)), i));
    end
    plot(meanrt, acc, '.', 'MarkerSize',mksz*3, 'color',mycol(1+(i-2)*10,:));
end
%% Wong06, 7 params, fitting results
dat_dir = '/gpfs/data/glimcherlab/BoShen/RecurrentModel/FitDynamic/Rslts/WW06FitBhvr7ParamsII_QMLE_GPU(OLD)/graphics';
params = [0.428032	0.1131	0.350868	0.022471	105.691412	0.042461	0.014988	16877.754482];
name = sprintf('JNp%2.1f_JNn%1.2f_I0%1.2f_noise%1.2f_miu0%2.2f_nLL%4.1f',params);
load(fullfile(dat_dir,sprintf('PlotData_%s.mat',name)));
h = figure; hold on;
for i = 2:6
    chset1 = squeeze(choicemat(:,i) == 1);
    chset2 = squeeze(choicemat(:,i) == 2);
    qntls = quantile(rtmat(chset1 | chset2,i),9);
    qntls = [min(rtmat(chset1 | chset2,i)) qntls max(rtmat(chset1 | chset2,i))];
    meanrt = [];
    acc = [];
    for q = 1:(length(qntls)-1)
        meanrt(q) = mean(rtmat((chset1 | chset2) & squeeze(rtmat(:,i) <= qntls(q+1) & rtmat(:,i) >= qntls(q)), i));
        acc(q) = mean(2-choicemat((chset1 | chset2) & squeeze(rtmat(:,i) <= qntls(q+1) & rtmat(:,i) >= qntls(q)), i));
    end
    plot(meanrt, acc, '.', 'MarkerSize',mksz*3, 'color',mycol(1+(i-2)*10,:));
end
