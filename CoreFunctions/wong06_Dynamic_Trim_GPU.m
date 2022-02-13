function [rt, choice, argmaxR, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = wong06_Dynamic_Trim_GPU(cp,miu0,sgm,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims, dot_ax, sac_ax)
%%%%%%%%%%%% 
% adapted version of the appendix function in Wong & Wang, 2006
% Created by Bo Shen, NYU, 2019
%% GPU calculation, only binary choice is allowed. N is limited as 2.
%%%%%%%%%%%%%%%%%%%%%
%%% CAUTION: because of lacking gpu memory, firing rates are not exported as time series
%%%%%%%%%%%%%%%%%%%%%
% - Vinput: allows two formats.
%   Format A: input values as a MxN array. N is the number of choice items,
%       M is the number of c1-c2 pairs
%       - example      Vinput =  [320, 192
%                                 288, 224
%                                 270, 240];
%    Format B: input values as a structure with fields defining c1, c2, c3 matrices.
%       The input matrices have to be the same dimension of 2. Values at the same position
%       of the matrices will be used as a set of inputs.
%       - example      [V1 V2] =  meshgrid([320,288,270],[192,224,240])
%                      Vinput.V1 = V1; Vinput.V2 = V2;
%%%%%%%%%%%%%%%
%% set parameters
% parameter in H
a = 270; %VnC^-1
b = 108; %Hz
d = .154; %s
% time unit
unit = 1; % set as 1 if the unit of time related variables is second.
% change it to .001 if the unit of time related varibales is msec.
JAext = 5.2 * 10^-4; %nA/Hz

%% define parameters
I0Array = gpuArray(I0);
JN11 = (JN(1,1));
JN12 = (JN(1,2));
JN21 = (JN(2,1));
JN22 = (JN(2,2));
if isstruct(cp)
    name = fieldnames(cp);
    cp1 = cp.(name{1});
    cp2 = cp.(name{2});
else
    cp1 = cp(:,1);
    cp2 = cp(:,2);
end
cp1Array = gpuArray(repmat(cp1,1,1,sims));
cp2Array = gpuArray(repmat(cp2,1,1,sims));
sizeVinput = size(cp1);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);

presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
% sliding window
time_wind = .050/dt/unit;  % Temporal window size for averaging, 50 msec
% not used, slide_wind = .005/dt/unit;  % Sliding step for window, 5 msec
% time line
posttask_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
onset_of_trigger = onset_of_stimuli;
if isequal(size(onset_of_trigger), size(V1mat))
    onset_of_trigger = gpuArray(repmat(onset_of_trigger,1,1,sims));
end
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
time_spc = 100; % ms, to exclude activity within 100 msecs of eye movement initiation in calculating mrc
time_spcD = 200; % ms, to exclude activity within 200 msecs of motion onset in calculating mrcD
%% stablizing noise for 200 ms
Inoise1 = gpuArray(zeros(sizeComput));
Inoise2 = gpuArray(zeros(sizeComput));
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    % update noise
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
end
%% get initial values at ti = 0
mr1c = gpuArray.nan([numel(dot_ax), sizeComput]); % mean rates sorted at beginning of the dot motion task
mr2c = gpuArray.nan([numel(dot_ax), sizeComput]);
mr1cD = gpuArray.nan([-min(sac_ax)+sum(sac_ax >= 0), sizeComput]); 
mr2cD = gpuArray.nan([-min(sac_ax)+sum(sac_ax >= 0), sizeComput]); % mean rates sorted at the time of decision,
% buffering sample rate before decision is at 1ms, after decision is based on sac_ax

% drop in buffer to sliding window
H1buffer = gpuArray(repmat(initialvals(1,1),[sizeComput,time_wind]));
H2buffer = gpuArray(repmat(initialvals(1,2),[sizeComput,time_wind]));
H1_wind = mean(H1buffer,4); % smoothed firing rates by sliding windows
H2_wind = mean(H2buffer,4);
H1new = gpuArray(ones(sizeComput)*initialvals(1,1));
H2new = gpuArray(ones(sizeComput)*initialvals(1,2));
S1 = gpuArray(ones(sizeComput)*initialvals(2,1));
S2 = gpuArray(ones(sizeComput)*initialvals(2,2));
%% simulation begin
R1Out = gpuArray.nan(sizeComput); % records the values at decision, or the end of simulation if choice was not made
R2Out = gpuArray.nan(sizeComput);
Continue = gpuArray(ones(sizeComput)); % to mark the trials that choices haven't made yet
rt = gpuArray.Inf(sizeComput);
choice = gpuArray(zeros(sizeComput));
tpafterward = gpuArray(zeros(sizeComput)); % time stamp intermediate variable after decision
for ti = 0:max(posttask_steps(:))
    % sample mrc according to dot_ax
    if any(ti == dot_ax)
        rec_axi = find(ti == dot_ax);
        stllgo = (choice == 0);
        mr1c(rec_axi,stllgo) = H1_wind(stllgo);  % replace the running trials' values, keep others as nan
        mr2c(rec_axi,stllgo) = H2_wind(stllgo);
    end

    % when all of the channel hit the decision boundary or timed out, stop simulation
    if stoprule == 1
        if NComput == 0 && all(tpafterward(:) > sum(sac_ax>0))
            break;
        end
    end
    
    % update the input values
    if stoprule == 1
        if NComput == 0 && all(tpafterward(:) > sum(sac_ax>0))
            break;
        end
    end
    
    if ti >= onset_of_stimuli && ti < offset_of_stimuli
        V1 = cp1Array;
        V2 = cp2Array;
    else
        V1 = 0;
        V2 = 0; 
    end
    % update nerual firing rates H
    H1 = H1new;
    H2 = H2new;
    I1 = JAext*miu0*V1;
    I2 = JAext*miu0*V2;
    x1 = JN11*S1 + JN12*S2 + I0Array + I1 + Inoise1;
    x2 = JN21*S1 + JN22*S2 + I0Array + I2 + Inoise2;
    H1new = (a*x1 - b)./(1 - exp(-d*(a*x1 - b)));
    H2new = (a*x2 - b)./(1 - exp(-d*(a*x2 - b)));
    % update synaptic activities S
    S1 = S1 + (-S1/tauS + (1-S1).*H1*unit*gamma)*dt;
    S2 = S2 + (-S2/tauS + (1-S2).*H2*unit*gamma)*dt;
    % update noise
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    % To ensure firing rates are always positive (noise may cause negative)
    inside = H1new >= 0;
    H1new = H1new .* inside;
    inside = H2new >= 0;
    H2new = H2new .* inside;
    loci = mod(ti,time_wind) + (mod(ti,time_wind) == 0)*time_wind;
    H1buffer(:,:,:,loci) = H1;
    H2buffer(:,:,:,loci) = H2;
    H1_wind = mean(H1buffer,4); % smoothed firing rates by sliding windows
    H2_wind = mean(H2buffer,4);
    % threshold detecting
    inside = (H1_wind >= thresh) + (H2_wind >= thresh);
    flip = (inside > 0) & (rt == Inf);
    NComput = NComput - sum(flip(:));
    if numel(unique(onset_of_trigger)) == 1
        rt(flip) = ti-onset_of_trigger(1);
    elseif all(size(onset_of_trigger) == sizeComput)
        rt(flip) = ti - onset_of_trigger(flip);
    end
    choice = choice + ((H2_wind > H1_wind) - (H1_wind > H2_wind) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    R1Out(flip) = H1_wind(flip); % update the values at choice, keep others as nan
    R2Out(flip) = H2_wind(flip);
    % exclude data in mrc 100 ms before decision (defined in time_spc)
    loc = find(flip);
    for excl = 1:numel(loc)
        if numel(unique(onset_of_trigger)) == 1
            excldt = rt(loc(excl)) + onset_of_trigger(1) - time_spc;
        elseif all(size(onset_of_trigger) == sizeComput)
            excldt = rt(loc(excl)) + onset_of_trigger(loc(excl)) - time_spc;
        end
        mr1c(dot_ax > excldt,loc(excl)) = NaN;
        mr2c(dot_ax > excldt,loc(excl)) = NaN;
    end
    % sample mrcD according to sac_ax
    if ti > onset_of_stimuli + time_spcD % wait for recording until 200ms after stimuli on
        % so automatically exclude data in m_mrcD within 200 ms of onset_of_stimuli (defined in time_spcD)
        % for chains are still going, update the values saved before sac
        mr1cD(1:(-min(sac_ax)-1), Continue) = mr1cD(2:-min(sac_ax),Continue); % push the values into the queue
        mr1cD(-min(sac_ax),Continue) = R1(Continue); % update R1 from saved buffer
        mr2cD(1:(-min(sac_ax)-1), Continue) = mr2cD(2:-min(sac_ax),Continue); % push into the queue
        mr2cD(-min(sac_ax),Continue) = R2(Continue);
        % for chains just made decision
        mr1cD(-min(sac_ax)+1,flip) = H1_wind(flip);
        mr2cD(-min(sac_ax)+1,flip) = H2_wind(flip);
        tpafterward(flip) = 1; % mark the time stamp at decision as 0, after decision, push the time stamp one step forward
        % for chains already stopped, sample according to sac_ax untill max(sac_ax)
        smpl = (ti == onset_of_trigger + rt + sac_ax(find(sac_ax == 0)+min(tpafterward,sum(sac_ax>0)))) & (tpafterward <= sum(sac_ax>0));
        tplist = unique(tpafterward(smpl));
        for si = 1:numel(tplist) % loop over different time stamps
            updatecells = smpl & (tpafterward == tplist(si));
            mr1cD(-min(sac_ax)+1+tplist(si), updatecells) = H1_wind(updatecells);
            mr2cD(-min(sac_ax)+1+tplist(si), updatecells) = H2_wind(updatecells);
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
%% post calculus
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
        mRT = median(rt(i,j,choose1)) + round(mean(triggert - presentt,'all')/dt); % median rt for trials choosing R1,  noticing rt = Inf in non-choice trials
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