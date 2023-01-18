function [ll, Chi2, N] = DynamicLLFun(m_mr1ccut,m_mr2ccut,m_mr1cDcut,m_mr2cDcut,sm_mr1c,sm_mr2c,sm_mr1cD,sm_mr2cD, dot_axcut, dot_gap, sac_axcut, sac_gap)
%% least squared error fitting for neural dynamics
Chi2 = 0;
N = 0;
for i = 1:4
    switch i
        case 1
            fit = sm_mr1c(dot_axcut>=dot_gap,:);
            dat = m_mr1ccut(dot_axcut>=dot_gap,:);
        case 2
            fit = sm_mr2c(dot_axcut>=dot_gap,:);
            dat = m_mr2ccut(dot_axcut>=dot_gap,:);
        case 3
            fit = sm_mr1cD(sac_axcut<=-sac_gap,:);
            dat = m_mr1cDcut(sac_axcut<=-sac_gap,:);
        case 4
            fit = sm_mr2cD(sac_axcut<=-sac_gap,:);
            dat = m_mr2cDcut(sac_axcut<=-sac_gap,:);
    end
    fit(isnan(fit)) = 0; % replace the missing values as zero, which will drive a large deviance from the data values
    sse = (fit - dat).^2;
    Nsse = sum(~isnan(sse(:)));
    Chi2 = Chi2 + sum(sse(:), 'omitnan');
    N = N + Nsse;
end
sigma2 = Chi2/N; % the normalization parameter in gaussian distribution, which connect the root mean squared error and lieklihood
ll = -.5*Chi2/sigma2 - .5*N*log(2*pi*sigma2);
% equivalently, ll1 can be written as
% ll1 = -N/2*log(Chi2/N) - N/2*log(2*pi) - N/2;
% maximizing likelihood is equivalent to minimizing RMSE
end