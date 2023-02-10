function [rt, choice, argmaxR, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = LCA_Dynmc_Trim_GPU(RhoMat, k, beta, sgm, T0, thresh, dur, dt, sims, dot_ax, sac_ax)
%% define parameters
sgmArray = gpuArray(sgm);
tau = gpuArray(.1);
dtArray = gpuArray(dt);
dt_tau = dt/tau;
k11 = gpuArray(k(1,1));
k12 = gpuArray(k(1,2));
k21 = gpuArray(k(2,1));
k22 = gpuArray(k(2,2));
beta11 = gpuArray(beta(1,1));
beta12 = gpuArray(beta(1,2));
beta21 = gpuArray(beta(2,1));
beta22 = gpuArray(beta(2,2));
threshArray = gpuArray(thresh);
if isstruct(RhoMat)
    name = fieldnames(RhoMat);
    Rho1mat = RhoMat.(name{1});
    Rho2mat = RhoMat.(name{2});
else
    Rho1mat = RhoMat(:,1);
    Rho2mat = RhoMat(:,2);
end
Rho1inputArray = gpuArray(repmat(Rho1mat,1,1,sims));
Rho2inputArray = gpuArray(repmat(Rho2mat,1,1,sims));
sizeVinput = size(Rho1mat);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);

time_spc = 100 -30; % ms, to exclude activity within 100 msecs of eye movement initiation in calculating mrc
time_spcD = 200 - 90; % ms, to exclude activity within 200 msecs of motion onset in calculating mrcD

%% initialize variables
rt = gpuArray.Inf(sizeComput);
choice = gpuArray.zeros(sizeComput);
x1Out = gpuArray.nan(sizeComput); % records the values at decision, or the end of simulation if choice was not made
x2Out = gpuArray.nan(sizeComput);
x1 = gpuArray.zeros(sizeComput);
x2 = gpuArray.zeros(sizeComput);

mr1c = gpuArray.nan([numel(dot_ax), sizeComput]); % mean rates sorted at beginning of the dot motion task
mr2c = gpuArray.nan([numel(dot_ax), sizeComput]);
mr1cD = gpuArray.nan([-min(sac_ax)+sum(sac_ax >= 0), sizeComput]); 
mr2cD = gpuArray.nan([-min(sac_ax)+sum(sac_ax >= 0), sizeComput]); % mean rates sorted at the time of decision,
% buffering sample rate before decision is at 1ms, after decision is based on sac_ax


%% simulation: premotion stage (ti < 0), dot motion task begin at ti = 0
% ti = 0 sorted at the the beginning of evidence accumulation, T0 capture the
% non-decision delay in the beginning and the execution time at the end
Continue = gpuArray(ones(sizeComput)); % to mark the trials that choices haven't made yet
tpafterward = gpuArray(zeros(sizeComput)); % time stamp intermediate variable after decision
for ti = 1:(dur/dt)
    % sample mrc according to dot_ax
    if any(ti == dot_ax)
        rec_axi = find(ti == (dot_ax));
        stllgo = (rt == Inf);
        mr1c(rec_axi,stllgo) = x1(stllgo);  % replace the running trials' values, keep others as nan
        mr2c(rec_axi,stllgo) = x2(stllgo);
    end

    % when all of the channel hit the decision boundary or timed out, stop simulation
    if NComput == 0 && all(tpafterward(:) > sum(sac_ax>0))
        break;
    end

    
    % update xs
    x1 = x1 + (Rho1inputArray.*Continue - k11*x1 - beta12*x2)*dt_tau + gpuArray.randn(sizeComput)*sgmArray*sqrt(dt_tau);
    x2 = x2 + (Rho2inputArray.*Continue - k22*x2 - beta21*x1)*dt_tau + gpuArray.randn(sizeComput)*sgmArray*sqrt(dt_tau);
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = x1 >= 0;
    x1 = x1 .* inside;
    inside = x2 >= 0;
    x2 = x2 .* inside;
    % threshold detecting
    inside = (x1 >= threshArray) + (x2 >= threshArray);
    flip = (inside > 0) & (choice == 0);
    NComput = NComput - sum(flip(:));
    rt(flip) = ti;
    
    choice = choice + ((x2 > x1) - (x1 > x2) + 3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    x1Out(flip) = x1(flip); % update the values at choice, keep others as nan
    x2Out(flip) = x2(flip);
    % exclude data in mrc 100 ms before decision (defined in time_spc)
    loc = find(flip);
    for excl = 1:numel(loc)
        excldt = rt(loc(excl)) - time_spc;
        mr1c(dot_ax > excldt,loc(excl)) = NaN;
        mr2c(dot_ax > excldt,loc(excl)) = NaN;
    end
    % sample mrcD according to sac_ax
    if ti > time_spcD % wait for recording until 200ms after stimuli on
        % so automatically exclude data in m_mrcD within 200 ms of onset_of_stimuli (defined in time_spcD)
        % for chains are still going, update the values saved before sac
        mr1cD(1:(-min(sac_ax)-1), Continue) = mr1cD(2:-min(sac_ax),Continue); % push the values into the queue
        mr1cD(-min(sac_ax),Continue) = x1(Continue); % update R1 from saved buffer
        mr2cD(1:(-min(sac_ax)-1), Continue) = mr2cD(2:-min(sac_ax),Continue); % push into the queue
        mr2cD(-min(sac_ax),Continue) = x2(Continue);
        % for chains just made decision
        mr1cD(-min(sac_ax)+1,flip) = x1(flip);
        mr2cD(-min(sac_ax)+1,flip) = x2(flip);
        tpafterward(flip) = 1; % mark the time stamp at decision as 0, after decision, push the time stamp one step forward
        % for chains already stopped, sample according to sac_ax untill max(sac_ax)
        smpl = (ti == rt + sac_ax(sum(sac_ax<=0)+min(tpafterward,sum(sac_ax>0)))) & (tpafterward <= sum(sac_ax>0));
        tplist = unique(tpafterward(smpl));
        for si = 1:numel(tplist) % loop over different time stamps
            updatecells = smpl & (tpafterward == tplist(si));
            mr1cD(-min(sac_ax)+1+tplist(si), updatecells) = x1(updatecells);
            mr2cD(-min(sac_ax)+1+tplist(si), updatecells) = x2(updatecells);
        end
        tpafterward = tpafterward + smpl; % push forward the time step after recorded
    end
end
% for those trials that choices were not made
x1Out(choice == 0) = x1(choice == 0);
x2Out(choice == 0) = x2(choice == 0);
%% calculate
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2,  0 choice is not made
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
inside = x1Out < x2Out;
inside = inside + (x1Out == x2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
argmaxR = gather(inside) + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
%% post calculus, only preserve the correct trials
m_mr1c = gpuArray.nan([numel(dot_ax), sizeComput(1), sizeComput(2)]); % average over all trials of the same condition
m_mr2c = gpuArray.nan([numel(dot_ax), sizeComput(1), sizeComput(2)]);
m_mr1cD = gpuArray.nan([numel(sac_ax), sizeComput(1), sizeComput(2)]);
m_mr2cD = gpuArray.nan([numel(sac_ax), sizeComput(1), sizeComput(2)]);
for i = 1:size(rt,1)
    for j = 1:size(rt,2)
        choose1 = choice(i,j,:) == 1; % mr1cD(memolen,i,j,:) > mr2cD(memolen,i,j,:); % argmaxR(i,j,:) == 1;
        m_mr1c(:,i,j) = mean(mr1c(:,i,j,choose1),4,'omitnan');
        m_mr2c(:,i,j) = mean(mr2c(:,i,j,choose1),4,'omitnan');
        m_mr1cD(:,i,j) = mean(mr1cD([-min(sac_ax) + sac_ax(sac_ax < 0)' + 1, -min(sac_ax)+(1:sum(sac_ax>=0))],i,j,choose1),4,'omitnan');
        m_mr2cD(:,i,j) = mean(mr2cD([-min(sac_ax) + sac_ax(sac_ax < 0)' + 1, -min(sac_ax)+(1:sum(sac_ax>=0))],i,j,choose1),4,'omitnan');
        % only look at the data more than half numbers of trials
        mRT = median(rt(i,j,choose1)); % median rt for trials choosing R1,  noticing rt = Inf in non-choice trials
        m_mr1c(dot_ax >= mRT + time_spc,i,j) = NaN;
        m_mr2c(dot_ax >= mRT + time_spc,i,j) = NaN;
        m_mr1cD(sac_ax <= -mRT + time_spcD,i,j) = NaN;
        m_mr2cD(sac_ax <= -mRT + time_spcD,i,j) = NaN;
    end
end
m_mr1c = gather(squeeze(m_mr1c));
m_mr2c = gather(squeeze(m_mr2c));
m_mr1cD = gather(squeeze(m_mr1cD));
m_mr2cD = gather(squeeze(m_mr2cD));
rt = gather(rt)*dt;
rt(rt == Inf) = nan;
choice = gather(choice);