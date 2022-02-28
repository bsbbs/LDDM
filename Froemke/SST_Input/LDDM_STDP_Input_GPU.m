function [rt, choice, argmaxR] = LDDM_STDP_Input_GPU(Vprior, Vinput, iSTDP, w, a, b,...
    sgm, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims)
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
V1inputArray = gpuArray(repmat(V1mat,1,1,sims));
V2inputArray = gpuArray(repmat(V2mat,1,1,sims));
V1prArray = gpuArray(repmat(V1prmat,1,1,sims));
V2prArray = gpuArray(repmat(V2prmat,1,1,sims));
sizeVinput = size(V1mat);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);
if all(size(iSTDP)==sizeVinput)
    iSTDP = gpuArray(repmat(iSTDP,1,1,sims));
end
pretask_steps = round(predur/dt);
onset_of_stimuli = gpuArray(round(presentt/dt));
if isequal(size(onset_of_stimuli), size(V1mat))
    onset_of_stimuli = gpuArray(repmat(onset_of_stimuli,1,1,sims));
end
onset_of_trigger = gpuArray(round(triggert/dt));
if isequal(size(onset_of_trigger), size(V1mat))
    onset_of_trigger = gpuArray(repmat(onset_of_trigger,1,1,sims));
end
posttask_steps = gpuArray(round(dur/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
if isequal(size(offset_of_stimuli), size(V1mat))
    offset_of_stimuli = gpuArray(repmat(offset_of_stimuli,1,1,sims));
end
%% stablizing noise for 200 ms
InoiseR1 = gpuArray.zeros(sizeComput);
InoiseG1 = gpuArray.zeros(sizeComput);
InoiseR2 = gpuArray.zeros(sizeComput);
InoiseG2 = gpuArray.zeros(sizeComput);
InoiseI1 = gpuArray.zeros(sizeComput);
InoiseI2 = gpuArray.zeros(sizeComput);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
X = gpuArray(ones(sizeComput));
R1 = (X*initialvals(1,1)) + InoiseR1;
R2 = (X*initialvals(1,2)) + InoiseR2;
G1 = (X*initialvals(2,1)) + InoiseG1;
G2 = (X*initialvals(2,2)) + InoiseG2;
I1 = (X*initialvals(3,1)) + InoiseI1;
I2 = (X*initialvals(3,2)) + InoiseI2;

%% initialize variables
rt = gpuArray.zeros(sizeComput);
choice = gpuArray.zeros(sizeComput);
R1Out = gpuArray.nan(sizeComput); % records the values at decision, or the end of simulation if choice was not made
R2Out = gpuArray.nan(sizeComput);
Continue = gpuArray(ones(sizeComput)); % to mark the trials that choices haven't made yet
BetasUp = gpuArray.zeros(sizeComput); % change from 0 to beta on the time of trigger
%% simulation: premotion stage (ti < 0), dot motion task begin at ti = 0
% ti = 0 sorted at the the beginning of the dot motion task, while the 
% onset of stimuli (onset_of_stimuli) can be later than 0, which will cause
% initial dip
for ti = -pretask_steps:max(posttask_steps(:))
    if stoprule == 1
        if NComput == 0
            break;
        end
    end
    % update the values
    if numel(unique(onset_of_stimuli)) == 1
        if ti >= onset_of_stimuli(1)
            V1Array = V1inputArray;
            V2Array = V2inputArray;
        else
            V1Array = 0;
            V2Array = 0;
        end
    elseif isequal(size(onset_of_stimuli), sizeComput)
        flip = ti >= onset_of_stimuli;
        V1Array = V1inputArray.*flip;
        V2Array = V2inputArray.*flip;
    end

    if numel(unique(offset_of_stimuli)) == 1
        if ti >= offset_of_stimuli(1)
            V1Array = 0;
            V2Array = 0;
        end
    elseif isequal(size(offset_of_stimuli), sizeComput)
        flip = ti >= offset_of_stimuli;
        V1Array(flip) = 0;
        V2Array(flip) = 0;
    end

    if numel(unique(onset_of_trigger)) == 1
        if ti >= onset_of_trigger(1)
            BetasUp = gpuArray(ones(sizeComput));
        end
    elseif isequal(size(onset_of_trigger), sizeComput)
        flip = ti >= onset_of_trigger;
        BetasUp(flip) = 1;
    end
    % update R, G, I
    G1old = G1;
    G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 - I1)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2  - I2)/Tau2*dtArray + InoiseG2;
    I1 = I1 + (-I1 + beta11*R1.*Continue.*BetasUp + beta12*R2.*Continue.*BetasUp)/Tau3*dtArray + InoiseI1;
    I2 = I2 + (-I2 + beta21*R1.*Continue.*BetasUp + beta22*R2.*Continue.*BetasUp)/Tau3*dtArray + InoiseI2;
    R1 = R1 + (-R1 + (V1Array + V1prArray*(ti<0) + alpha11*R1+alpha12*R2).*Continue./(1+iSTDP.*G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (V2Array + V2prArray*(ti<0) + alpha21*R1+alpha22*R2).*Continue./(1+iSTDP.*G2old))/Tau1*dtArray + InoiseR2;
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
    flip = (inside > 0) & (choice == 0) & (ti > onset_of_trigger);
    NComput = NComput - sum(flip(:));
    rt = rt + gpuArray(ti-onset_of_trigger).*flip*dtArray;
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    R1Out(flip) = R1(flip); % update the values at choice, keep others as nan
    R2Out(flip) = R2(flip);
end
% for those trials that choices were not made
R1Out(choice == 0) = R1(choice == 0);
R2Out(choice == 0) = R2(choice == 0);
%% calculate
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2, 0 choice is not made
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
rt(rt==0) = NaN;
inside = R1Out < R2Out;
inside = inside + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
argmaxR = inside + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
rt = gather(rt);
choice = gather(choice);
argmaxR = gather(argmaxR);