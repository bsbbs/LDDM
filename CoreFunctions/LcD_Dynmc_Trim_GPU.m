function [rt, choice, argmaxR, m_mr1c, m_mr2c, m_mr1cD, m_mr2cD] = LcD_Dynmc_Trim_GPU(Vprior, Vinput, w, a, b,...
    sgm, Tau, predur, dur, dt, tgap, triggert, thresh, initialvals, stimdur, stoprule, sims, dot_ax, sac_ax)
%%%%%%%%%%%%%%%%%%
%% GPU calculation, only binary choice is allowed. N is limited as 2.
% Created by Bo Shen, NYU, 2019
% Vinput, w, a (alpha), b (beta), Tau, dt are all in the same meaning as
% the names of parameters in the paper.
% - Vinput: allows two formats.
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
%%%%%%%%%%%%%%%%%%%
tauN = 0.002; % time constant for Ornstein-Uhlenbeck process of noise
%% preparation
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
if isstruct(Vinput)
    name = fieldnames(Vinput);
    V1mat = Vinput.(name{1});
    V2mat = Vinput.(name{2});
    name = fieldnames(Vprior);
    V1prmat = Vprior.(name{1});
    V2prmat = Vprior.(name{2});
else
    V1mat = Vinput(:,1);
    V2mat = Vinput(:,2);
    V1prmat = Vprior(:,1);
    V2prmat = Vprior(:,2);
end
sizeVinput = size(V1mat);
sizeComput = [sizeVinput, sims];
X = ones(sizeComput);
NComput = prod(sizeComput);
premotion_steps = round(predur/dt);
presentt = 0;
onset_of_stimuli = gpuArray(round(tgap/dt));
onset_of_trigger = gpuArray(round(triggert/dt));
total_time_steps = gpuArray(round(dur/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = 0 + stim_duration;
% dot_ax = [-100:20:1000];
% sac_ax = [-1000:20:300];
step = dot_ax(2) - dot_ax(1); % ms, time step for sampling
time_spc = 100; % ms, time space leaving for 100ms after decision in mrc
time_spcD = 200; % ms, time space leaving for 200ms after stimui on in mrcD
before = -min(dot_ax); % ms, for mrc
after = max(sac_ax); % ms, for mrcD
bufflen = time_spc/step; % for mrc
bufflenD = time_spcD/step; % for mrcD
mrclen = numel(dot_ax); % for mrc
memolen = abs(min(sac_ax)/step); % for mrcD
mmrllen = 1+abs(max(sac_ax)/step); % for mrcD
%% stablizing noise
InoiseR1=gpuArray(X*0);
InoiseG1=gpuArray(X*0);
InoiseI1=gpuArray(X*0);
InoiseR2=gpuArray(X*0);
InoiseG2=gpuArray(X*0);
InoiseI2=gpuArray(X*0);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
R1 = gpuArray(X*initialvals(1,1)) + InoiseR1;
R2 = gpuArray(X*initialvals(1,2)) + InoiseR2;
G1 = gpuArray(X*initialvals(2,1)) + InoiseG1;
G2 = gpuArray(X*initialvals(2,2)) + InoiseG2;
I1 = gpuArray(X*initialvals(3,1)) + InoiseI1;
I2 = gpuArray(X*initialvals(3,2)) + InoiseI2;
%% simulation
R1slide = gpuArray.nan([step,sizeComput]);
R2slide = gpuArray.nan([step,sizeComput]);
R1buff = gpuArray.nan([bufflen,sizeComput]);
R2buff = gpuArray.nan([bufflen,sizeComput]);
mr1c = gpuArray.nan([mrclen, sizeComput]);
mr2c = gpuArray.nan([mrclen, sizeComput]);
mr1cD = gpuArray.nan([memolen+mmrllen, sizeComput]);
mr2cD = gpuArray.nan([memolen+mmrllen, sizeComput]);
m_mr1c = gpuArray.nan([mrclen, sizeComput(1), sizeComput(2)]);
m_mr2c = gpuArray.nan([mrclen, sizeComput(1), sizeComput(2)]);
m_mr1cD = gpuArray.nan([memolen+mmrllen, sizeComput(1), sizeComput(2)]);
m_mr2cD = gpuArray.nan([memolen+mmrllen, sizeComput(1), sizeComput(2)]);
rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
%% premotion stage,
V1prArray = gpuArray(repmat(V1prmat,1,1,sims));
V2prArray = gpuArray(repmat(V2prmat,1,1,sims));
for ti = -premotion_steps:0
    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 - I1)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2  - I2)/Tau2*dtArray + InoiseG2;
    I1 = I1 + (-I1)/Tau3*dtArray + InoiseI1;
    I2 = I2 + (-I2)/Tau3*dtArray + InoiseI2;
    R1 = R1 + (-R1 + (V1prArray + alpha11*R1+alpha12*R2)./(1+G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (V2prArray + alpha21*R1+alpha22*R2)./(1+G2old))/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = I1 >= 0;
    I1 = I1 .* inside;
    inside = I2 >= 0;
    I2 = I2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
    % slide time bin, 20 ms
    stepi = mod(ti, step)+1;
    R1slide(stepi,:,:,:) = R1;
    R2slide(stepi,:,:,:) = R2;
    % current time step for mrc
    % dot_axi = ceil((ti - (triggert/dt - before))/step);
    dot_axi = ceil((ti - (presentt - before))/step); 
    if dot_axi > 0 && dot_axi <= numel(dot_ax)+bufflen && stepi == 1 % refresh records for mrc
        buffi = mod(dot_axi, bufflen) + 1;
        rec_axi = dot_axi - bufflen; % it was 100ms or 5 steps before the current time stamp
        if rec_axi > 0
            stllgo = rt == 0;
            mr1c(rec_axi,stllgo) = R1buff(buffi,stllgo);  % count trails that rt still hasn't made 100ms ago
            mr2c(rec_axi,stllgo) = R2buff(buffi,stllgo);
        end
        R1buff(buffi,:,:,:) = mean(R1slide,1); % update R buffer
        R2buff(buffi,:,:,:) = mean(R2slide,1);
    end
%     % current time step for mrc
%     dot_axi = ceil((ti - (triggert/dt - before))/step); 
%     if dot_axi > 0 && dot_axi <= numel(dot_ax)+bufflen && stepi == 1 % refresh records for mrc
%         buffi = mod(dot_axi, bufflen) + 1;
%         rec_axi = dot_axi - bufflen; % it was 100ms or 5 steps before the current time stamp
%         if rec_axi > 0
%             stllgo = squeeze(rt == 0);
%             dim = numel(size(stllgo));
%             m_mr1c(rec_axi,:,:) = mean(squeeze(R1buff(buffi,:,:,:)).*stllgo,dim)./mean(stllgo,dim); % count trails that rt still hasn't made
%             m_mr2c(rec_axi,:,:) = mean(squeeze(R2buff(buffi,:,:,:)).*stllgo,dim)./mean(stllgo,dim);
%         end
%         R1buff(buffi,:,:,:) = mean(R1slide,1); % update R buffer
%         R2buff(buffi,:,:,:) = mean(R2slide,1);
%     end
end
%% ndt gap, V = 0, beta = 0
for ti = 1:(onset_of_stimuli-1)
    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 - I1)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2 - I2)/Tau2*dtArray + InoiseG2;
    I1 = I1 + (-I1)/Tau3*dtArray + InoiseI1;
    I2 = I2 + (-I2)/Tau3*dtArray + InoiseI2;
    R1 = R1 + (-R1 + (alpha11*R1+alpha12*R2)./(1+G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (alpha21*R1+alpha22*R2)./(1+G2old))/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = I1 >= 0;
    I1 = I1 .* inside;
    inside = I2 >= 0;
    I2 = I2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
    % slide time bin, 20 ms
    stepi = mod(ti, step)+1;
    R1slide(stepi,:,:,:) = R1;
    R2slide(stepi,:,:,:) = R2;
    % current time step for mrc
    % dot_axi = ceil((ti - (triggert/dt - before))/step);
    dot_axi = ceil((ti - (presentt - before))/step); 
    if dot_axi > 0 && dot_axi <= numel(dot_ax)+bufflen && stepi == 1 % refresh records for mrc
        buffi = mod(dot_axi, bufflen) + 1;
        rec_axi = dot_axi - bufflen; % it was 100ms or 5 steps before the current time stamp
        if rec_axi > 0
            stllgo = rt == 0;
            mr1c(rec_axi,stllgo) = R1buff(buffi,stllgo);  % count trails that rt still hasn't made 100ms ago
            mr2c(rec_axi,stllgo) = R2buff(buffi,stllgo);
        end
        R1buff(buffi,:,:,:) = mean(R1slide,1); % update R buffer
        R2buff(buffi,:,:,:) = mean(R2slide,1);
    end
end
if isempty(ti)
    ti0 = 0;
else
    ti0 = ti;
end
%% 1st stage, after presenting, before decision, V turned on, working memory (alpha) up, beta still = 0
V1Array = gpuArray(repmat(V1mat,1,1,sims));
V2Array = gpuArray(repmat(V2mat,1,1,sims));
if numel(stimdur) == 1 % to cover the situation when onset_of_stimuli = offset_of_stimuli
    if stimdur < dt
        V1Array = gpuArray(zeros(sizeComput));
        V2Array = gpuArray(zeros(sizeComput));
    end
elseif prod(size(stimdur) == size(V1mat))
    stimdurmat = repmat(stimdur,1,1,sims);
    V1Array(stimdurmat < dt) = 0;
    V2Array(stimdurmat < dt) = 0;
    offset_of_stimulimat = repmat(offset_of_stimuli,1,1,sims);
else 
    error('the size of stimdur is different from the size of Vinput, plz check');
end
R1Out = gpuArray(zeros(sizeComput));
R2Out = gpuArray(zeros(sizeComput));
Continue = gpuArray(ones(sizeComput));
BetasUp = gpuArray(zeros(sizeComput));
for ti = (ti0 + 1):max(total_time_steps(:)) %(onset_of_trigger -1)
    if stoprule == 1
        if NComput == 0 && ti > round((triggert + max(rt(:)))/dt + after + step)
            break;
        end
    end
    if numel(offset_of_stimuli) == 1
        if ti == offset_of_stimuli
            V1Array = gpuArray(zeros(sizeComput));
            V2Array = gpuArray(zeros(sizeComput));
        end
    elseif prod(size(offset_of_stimuli) == size(V1mat))
        flip = offset_of_stimulimat == ti;
        V1Array(flip) = 0;
        V2Array(flip) = 0;
        % R1Rprsnt(flip) = R1(flip);
        % R2Rprsnt(flip) = R2(flip);
    end
    if numel(onset_of_trigger) == 1
        if ti >= onset_of_trigger
            BetasUp = gpuArray(ones(sizeComput));
        end
    elseif prod(size(onset_of_trigger) == size(V1mat))
        flip = TI >= onset_of_trigger;
        BetasUp(flip) = 1;
    end
    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 - I1)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2  - I2)/Tau2*dtArray + InoiseG2;
    I1 = I1 + (-I1 + beta11*R1.*Continue.*BetasUp + beta12*R2.*Continue.*BetasUp)/Tau3*dtArray + InoiseI1;
    I2 = I2 + (-I2 + beta21*R1.*Continue.*BetasUp + beta22*R2.*Continue.*BetasUp)/Tau3*dtArray + InoiseI2;
    R1 = R1 + (-R1 + (V1Array + alpha11*R1+alpha12*R2).*Continue./(1+G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (V2Array + alpha21*R1+alpha22*R2).*Continue./(1+G2old))/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = I1 >= 0;
    I1 = I1 .* inside;
    inside = I2 >= 0;
    I2 = I2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
    % threshold detecting
    inside = (R1 >= threshArray) + (R2 >= threshArray);
    flip = (inside > 0) .* (rt == 0) .* (ti > onset_of_trigger);
    NComput = NComput - sum(flip(:));
    rt = rt + gpuArray(ti-onset_of_trigger).*flip*dtArray;
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    R1Out = R1Out + R1.*flip;
    R2Out = R2Out + R2.*flip;
    % slide time bin, 20 ms
    stepi = mod(ti, step)+1;
    R1slide(stepi,:,:,:) = R1;
    R2slide(stepi,:,:,:) = R2;
    % current time step for mrc
    % dot_axi = ceil((ti - (triggert/dt - before))/step);
    dot_axi = ceil((ti - (presentt - before))/step); 
    if dot_axi > 0 && dot_axi <= numel(dot_ax)+bufflen && stepi == 1 % refresh records for mrc
        buffi = mod(dot_axi, bufflen) + 1;
        rec_axi = dot_axi - bufflen; % it was 100ms or 5 steps before the current time stamp
        if rec_axi > 0
            stllgo = rt == 0;
            mr1c(rec_axi,stllgo) = R1buff(buffi,stllgo);  % count trails that rt still hasn't made 100ms ago
            mr2c(rec_axi,stllgo) = R2buff(buffi,stllgo);
        end
        R1buff(buffi,:,:,:) = mean(R1slide,1); % update R buffer
        R2buff(buffi,:,:,:) = mean(R2slide,1);
    end
    % current time step for mrcD
    if dot_axi > 0 && stepi == 1
        GtherR1 = mean(R1slide,1); % average R1 from the sliding window
        GtherR2 = mean(R2slide,1); % average R2 from the sliding window
        % for chains are still going, update the values saved before sac
        delay = ti - triggert/dt - bufflenD*step; % wait for 200ms or 10 steps after stimuli on
        if delay > 0
            stllgo = rt == 0;
            mr1cD(1:memolen-1, stllgo) = mr1cD(2:memolen,stllgo); % push into the queue
            mr1cD(memolen,stllgo) = GtherR1(stllgo); % update R1 from saved buffer
            mr2cD(1:memolen-1, stllgo) = mr2cD(2:memolen,stllgo); % push into the queue
            mr2cD(memolen,stllgo) = GtherR2(stllgo);
        end
        % for chains already stopped, continue saving untill max(sac_ax)
        mmrial = ceil((ti - (triggert+rt)/dt + .5)/step).*(rt > 0);
        callist = gather(unique(mmrial(mmrial>0 & mmrial<=mmrllen)));
        for fi = 1:length(callist)
            cal = callist(fi);
            mr1cD(memolen+cal,mmrial == cal) = GtherR1(mmrial == cal);
            mr2cD(memolen+cal,mmrial == cal) = GtherR2(mmrial == cal);
        end
    end
end
R1Out = R1Out + R1.*(R1Out == 0);
R2Out = R2Out + R2.*(R2Out == 0);
%% calculate
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2,  0 choice is not made
rt(rt==0) = NaN;
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
inside = R1Out < R2Out;
inside = inside + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
argmaxR = inside + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
% argmaxR = 2 - mean(inside,3);
rt = gather(rt);
choice = gather(choice);
%% post calculus
for i = 1:size(rt,1)
    for j = 1:size(rt,2)
        choose1 = mr1cD(memolen,i,j,:) > mr2cD(memolen,i,j,:); % argmaxR(i,j,:) == 1;
        m_mr1c(:,i,j) = mean(mr1c(:,i,j,choose1),4,'omitnan');
        m_mr2c(:,i,j) = mean(mr2c(:,i,j,choose1),4,'omitnan');
        m_mr1cD(:,i,j) = mean(mr1cD(:,i,j,choose1),4,'omitnan');
        m_mr2cD(:,i,j) = mean(mr2cD(:,i,j,choose1),4,'omitnan');
        mRT = median(rt(i,j,choose1),'omitnan'); % median rt for trials choosing R1
        % trim closer to time_spc ot time_spcD
        if isnan(mRT)
            mRT = dur;
        end
        m_mr1c(round((before + (triggert+mRT)/dt-time_spc)/step):end,i,j) = NaN;
        m_mr2c(round((before + (triggert+mRT)/dt-time_spc)/step):end,i,j) = NaN;
        m_mr1cD(1:(memolen - round((mRT/dt - time_spcD)/step)),i,j) = NaN;
        m_mr2cD(1:(memolen - round((mRT/dt - time_spcD)/step)),i,j) = NaN;
    end
end
m_mr1c = gather(squeeze(m_mr1c));
m_mr2c = gather(squeeze(m_mr2c));
m_mr1cD = gather(squeeze(m_mr1cD));
m_mr2cD = gather(squeeze(m_mr2cD));