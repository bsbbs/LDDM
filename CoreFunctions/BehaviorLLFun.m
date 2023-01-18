function [ll, Chi2, N] = BehaviorLLFun(q, On ,ON, rtmat, choicemat)
% QMLE, quantile maximum likelihood estimation
% reference: Heathcote & Australia, and Mewhort, 2002.
% Chi-square, reference: Ratcliff & McKoon, 2007.
h = figure;
for vi = 1:6
    En(vi) = numel(rtmat(:,vi)); % the total numbers of trials in each condition
    RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    if ~all(isnan(q(:,1,vi)))
        tmp = histogram(RT_corr, [0; q(:,1,vi); Inf], 'Visible',0);
        EN(:,1,vi) = tmp.Values; % the number of correct trials in each condition
    else
        EN(:,1,vi) = NaN(numel(q(:,1,vi))+1,1);
    end
    if ~all(isnan(q(:,2,vi)))
        tmp = histogram(RT_wro, [0; q(:,2,vi); Inf], 'Visible',0);
        EN(:,2,vi) = tmp.Values; % the number of error trials in each condition
    else
        EN(:,2,vi) =  NaN(numel(q(:,2,vi))+1,1);
    end
    f(:,:,vi) = log((EN(:,:,vi)/En(vi))); % the logrithmed proportion of each condition, with correct and error trials seperately counted
    f(f(:,1,vi) == -Inf,1,vi) = log(.1/En(vi)); % set floor value of f at each point, to prevent -Inf
    f(f(:,2,vi) == -Inf,2,vi) = log(.1/En(vi)); % set floor value of f at each point, to prevent -Inf
    ON_adj(:,:,vi) = ON(:,:,vi)*En(vi)./On(vi); % the observed number of trials under each condition, scaled to match the total number of trials used in estimation
end
close(h);
ll = sum(ON(:).*f(:),'omitnan');
Chi2vec = (EN - ON_adj).^2./EN; % not working because of the EN = 0 problem
Chi2 = sum(Chi2vec(:),'omitnan');
N = sum(~isnan(ON(:).*f(:)));
end