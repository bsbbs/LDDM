function [rt, choice, argmaxR, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = LDDM_Full_Dynmc_gap_GPU(Vinput, w, a, b,...
    sgm, Tau, dur, dt, thresh, initialvals, G0, stoprule, sims, dot_ax, sac_ax, dot_gap, sac_gap)
%%%%%%%%%%%%%%%%%%
%% GPU calculation, only binary choice is allowed. N is limited as 2.
% Created by Bo Shen, NYU, 2019
% Vinput, w, a (alpha), b (beta), Tau, dt are all in the same meaning as
% the names of parameters in the paper.
% - Vinput and Vprior: allows two formats.
%   Format A: input values as a MxN array. N is the number of choice items,
%       M is the number of V1-V2 pairs
%       - example      Vinput =  [320, 192
%                                 288, 224
%                                 270, 240];
%    Format B: input values as a structure with the first field defining a V1
%       matrix and the second field defining a V2 matrix. V1 and V2 matrix
%       have to be the same dimension. V1 and V2 values at the same position
%       of the matrix will be paired as inputs.
%       - example      [V1 V2] =  meshgrid([320,288,270],[192,224,240])
%                      Vinput.V1 = V1; Vinput.V2 = V2;
% - w: connection weight from R neurons to G neurons. w must be a NxN matrix,
%   describing connections between each pair of R and G.
%   - example      [w11, w12;   [1, 1;
%                 w21, w22] =  1,1];
% - a: self-excitation parameter, describing the weight of
%   self-excitation for R neurons. alpha must be a NxN diagonal matrix.
%   Non-zero values only on the diagnal means that each R neuron has its own
%   excitation loop. Lateral or cross excitation does not exist.
%   - example      [a11, a12;   [15, 15;
%                 a21, a22] =  15, 15];
% - b: the weight of local disinhibition,
%   - example      [b11, b12;   [2, 2;
%                 b21, b22] =  2, 2];
% - sgm: the magnitude of background noise as a scalar, determines the
%   variance of Ornstein-Uhlenbeck process
% - Tau: a 1x3 array contains the time constants for R, G, and I neurons.
%   - example      Tau =  [.2,.2,.2];
% - dur: total duration of simulations, from the beginning of simulation to a
%   force ending, unit in second.
% - dt: the time step for numerical integration in unit of second, for example 0.001 (1 msec).
% - presentt: the time onset of stimuli presentation from the start of
%   simulation as time 0. The unit is second.
% - triggert: the time onset of turning on the local disinhibition from
%   the start of simulation as time 0. The weight of local disinhibition
%   turn from zero to positive values set in beta. The unit is
%   second.
% - thresh: threshold for R neural firing rates. Choice is assumed to be made
%   when any of the R neural firing rate hit the threshold.
% - initialvals: initial values for the neural firing rates of R, G, and I.
%   must be a 3xN matrix.
%   - example      initialvals = [R1_0, R2_0
%                                 G1_0, G2_0]
%                                 I1_0, I2_0]
% - stimdur: defines the duration of stimuli presenting, starting from presentt.
%   After withdrawal of stimuli, all input values turn into zeros.
% - stoprule:  a control parameter to tell the program if need
%   to stop simulating when one of the R neurons' firing rates hit threshold
% 	1 for to stop, 0 for to continue simulating until total duration set in dur.
% - sims: numbers of simulations for each pair of input
% - mrc: mean rate wrt dot task onset
% - mrcD: mean rate wrt decision
% - the time line is sorted at the beginning of the dot motion task. ti = 0
% indicated the task enters motion stage from premotion stage
%%%%%%%%%%%%%%%%%%%
tauN = 0.002; % time constant for Ornstein-Uhlenbeck process of noise
%% define parameters
sgmArray = gpuArray(sgm);
tauN = gpuArray(tauN);
dtArray = gpuArray(dt);
w11 = gpuArray(w(1,1));
w12 = gpuArray(w(1,2));
w21 = gpuArray(w(2,1));
w22 = gpuArray(w(2,2));
Tau1 = gpuArray(Tau(1));
Tau2 = gpuArray(Tau(2));
Tau3 = gpuArray(Tau(3));
alpha11 = gpuArray(a(1,1));
alpha12 = gpuArray(a(1,2));
alpha21 = gpuArray(a(2,1));
alpha22 = gpuArray(a(2,2));
beta11 = gpuArray(b(1,1));
beta12 = gpuArray(b(1,2));
beta21 = gpuArray(b(2,1));
beta22 = gpuArray(b(2,2));
threshArray = gpuArray(thresh);
G0 = gpuArray(G0);
if isstruct(Vinput)
    name = fieldnames(Vinput);
    V1mat = Vinput.(name{1});
    V2mat = Vinput.(name{2});
else
    V1mat = Vinput(:,1);
    V2mat = Vinput(:,2);
end
sizeVinput = size(V1mat);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);
V1Array = gpuArray(repmat(V1mat,1,1,sims));
V2Array = gpuArray(repmat(V2mat,1,1,sims));
posttask_steps = gpuArray(round(dur/dt));
dot_ax = gpuArray(dot_ax);
sac_ax = gpuArray(sac_ax);
sac_gap = gpuArray(sac_gap);
dot_gap = gpuArray(dot_gap);
sac_axp = sac_ax + sac_gap;
time_spc = gpuArray(100); % ms, to exclude activity within 100 msecs of eye movement initiation in calculating mrc
% since there is a delay (sac_gap) between the dynamic hit boundry and
% saccade behavior, the gap is counted.
time_spcD = gpuArray(200); % ms, to exclude activity within 200 msecs of motion onset in calculating mrcD
% since there is a delay (dot_gap) between the stimuli onset and the
% starting of WTA dynamics, the gap is counted.
%% stablizing noise for 200 ms
InoiseR1 = gpuArray(zeros(sizeComput));
InoiseG1 = gpuArray(zeros(sizeComput));
InoiseD1 = gpuArray(zeros(sizeComput));
InoiseR2 = gpuArray(zeros(sizeComput));
InoiseG2 = gpuArray(zeros(sizeComput));
InoiseD2 = gpuArray(zeros(sizeComput));
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseD1 = InoiseD1 + (-InoiseD1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseD2 = InoiseD2 + (-InoiseD2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
X = gpuArray(ones(sizeComput));
R1 = (X*initialvals(1,1)) + InoiseR1;
R2 = (X*initialvals(1,2)) + InoiseR2;
G1 = (X*initialvals(2,1)) + InoiseG1;
G2 = (X*initialvals(2,2)) + InoiseG2;
D1 = (X*initialvals(3,1)) + InoiseD1;
D2 = (X*initialvals(3,2)) + InoiseD2;

%% initialize variables
mr1c = gpuArray.nan([numel(dot_ax), sizeComput]); % mean rates sorted at beginning of the dot motion task
mr2c = gpuArray.nan([numel(dot_ax), sizeComput]);
mr1cD = gpuArray.nan([-min(sac_axp)+1+sum(sac_axp > 0), sizeComput]); 
mr2cD = gpuArray.nan([-min(sac_axp)+1+sum(sac_axp > 0), sizeComput]); % mean rates sorted at the time of decision,
% buffering sample rate before decision is at 1ms, after decision is based on sac_ax
rt = gpuArray.Inf(sizeComput);
choice = gpuArray(zeros(sizeComput));

%% simulation: premotion stage (ti < 0), dot motion task begin at ti = 0
% ti = 0 sorted at the the beginning of the dot motion task, while the 
% beginning of simulation will be later, determined by dot_gap
R1Out = gpuArray.nan(sizeComput); % records the values at decision, or the end of simulation if choice was not made
R2Out = gpuArray.nan(sizeComput);
Continue = gpuArray(ones(sizeComput)); % to mark the trials that choices haven't made yet
tpafterward = gpuArray(zeros(sizeComput)); % time stamp intermediate variable after decision
for ti = dot_gap:max(posttask_steps(:))
    % sample mrc according to dot_ax
    if any(ti == dot_ax)
        rec_axi = find(ti == (dot_ax));
        stllgo = (rt == Inf);
        mr1c(rec_axi,stllgo) = R1(stllgo);  % replace the running trials' values, keep others as nan
        mr2c(rec_axi,stllgo) = R2(stllgo);
    end

    % when all of the channel hit the decision boundary or timed out, stop simulation
    if stoprule == 1
        if NComput == 0 && all(tpafterward(:) > sum(sac_axp>0))
            break;
        end
    end

    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 - D1 + G0)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2  - D2 + G0)/Tau2*dtArray + InoiseG2;
    D1 = D1 + (-D1 + beta11*R1.*Continue + beta12*R2.*Continue)/Tau3*dtArray + InoiseD1;
    D2 = D2 + (-D2 + beta21*R1.*Continue + beta22*R2.*Continue)/Tau3*dtArray + InoiseD2;
    R1 = R1 + (-R1 + (V1Array + alpha11*R1+alpha12*R2).*Continue./(1+G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (V2Array + alpha21*R1+alpha22*R2).*Continue./(1+G2old))/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseD1 = InoiseD1 + (-InoiseD1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseD2 = InoiseD2 + (-InoiseD2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = D1 >= 0;
    D1 = D1 .* inside;
    inside = D2 >= 0;
    D2 = D2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
    % threshold detecting
    inside = (R1 >= threshArray) + (R2 >= threshArray);
    flip = (inside > 0) & (rt == Inf);
    NComput = NComput - sum(flip(:));
    rt(flip) = ti;
        
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    R1Out(flip) = R1(flip); % update the values at choice, keep others as nan
    R2Out(flip) = R2(flip);
    % exclude data in mrc 100 ms before decision (defined in time_spc)
    % count the gap from the boundary hitting to saccade
    loc = find(flip);
    for excl = 1:numel(loc)
        excldt = rt(loc(excl)) - time_spc + sac_gap;
        mr1c(dot_ax > excldt,loc(excl)) = NaN;
        mr2c(dot_ax > excldt,loc(excl)) = NaN;
    end
    % sample mrcD according to sac_ax
    if ti > time_spcD % wait for recording until 200ms after stimuli on
        % so automatically exclude data in m_mrcD within 200 ms of onset_of_stimuli (defined in time_spcD)
        % for chains are still going, update the values saved before sac
        mr1cD(1:(-min(sac_axp)-1), Continue) = mr1cD(2:-min(sac_axp),Continue); % push the values into the queue
        mr1cD(-min(sac_axp),Continue) = R1(Continue); % update R1 from saved buffer
        mr2cD(1:(-min(sac_axp)-1), Continue) = mr2cD(2:-min(sac_axp),Continue); % push into the queue
        mr2cD(-min(sac_axp),Continue) = R2(Continue);
        % for chains just made decision
        mr1cD(-min(sac_axp)+1,flip) = R1(flip);
        mr2cD(-min(sac_axp)+1,flip) = R2(flip);
        tpafterward(flip) = 1; % mark the time stamp at decision as 0, after decision, push the time stamp one step forward
        % for chains already stopped, sample according to sac_ax untill max(sac_ax)
        smpl = (ti == rt + sac_axp(sum(sac_axp<=0)+min(tpafterward,sum(sac_axp>0)))) & (tpafterward <= sum(sac_axp>0));
        tplist = unique(tpafterward(smpl));
        for si = 1:numel(tplist) % loop over different time stamps
            updatecells = smpl & (tpafterward == tplist(si));
            mr1cD(-min(sac_axp)+1+tplist(si), updatecells) = R1(updatecells);
            mr2cD(-min(sac_axp)+1+tplist(si), updatecells) = R2(updatecells);
        end
        tpafterward = tpafterward + smpl; % push forward the time step after recorded
    end
end
% for those trials that choices were not made
R1Out(choice == 0) = R1(choice == 0);
R2Out(choice == 0) = R2(choice == 0);
%% calculate
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2,  0 choice is not made
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
inside = R1Out < R2Out;
inside = inside + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
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
        m_mr1cD(:,i,j) = mean(mr1cD([(-min(sac_axp) + 1) + sac_axp(sac_axp <= 0)', -min(sac_axp)+1+(1:sum(sac_axp>0))],i,j,choose1),4,'omitnan');
        m_mr2cD(:,i,j) = mean(mr2cD([(-min(sac_axp) + 1) + sac_axp(sac_axp <= 0)', -min(sac_axp)+1+(1:sum(sac_axp>0))],i,j,choose1),4,'omitnan');
        % only look at the data more than half numbers of trials
        % mRT = median(rt(i,j,choose1)); % median rt for trials choosing R1,  noticing rt = Inf in non-choice trials
        % m_mr1c(dot_ax >= mRT + time_spc,i,j) = NaN;
        % m_mr2c(dot_ax >= mRT + time_spc,i,j) = NaN;
        % m_mr1cD(sac_ax <= -mRT + time_spcD,i,j) = NaN;
        % m_mr2cD(sac_ax <= -mRT + time_spcD,i,j) = NaN;
    end
end
m_mr1c = gather(squeeze(m_mr1c));
m_mr2c = gather(squeeze(m_mr2c));
m_mr1cD = gather(squeeze(m_mr1cD));
m_mr2cD = gather(squeeze(m_mr2cD));
rt = gather(rt)*dt + sac_gap*dt;
rt(rt == Inf) = nan;
choice = gather(choice);