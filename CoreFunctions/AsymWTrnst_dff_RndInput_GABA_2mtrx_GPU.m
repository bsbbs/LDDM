function [rt, choice, argmaxR, dR] = AsymWTrnst_dff_RndInput_GABA_2mtrx_GPU(Vinput, Gaba, w, wasym, a, sgm, sgmInput, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims)
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

Tau1 = gpuArray(Tau(1));
Tau2 = gpuArray(Tau(2));
threshArray = gpuArray(thresh);

V1mat = Vinput(1);
V2mat = Vinput(2);
n = size(a,1);
m = size(wasym,2);
sizeComput = [n, m, sims];
alpha11 = gpuArray(repmat(a,1,1,sims));
alpha22 = gpuArray(repmat(a,1,1,sims));
v = w;
X = ones(sizeComput);
NComput = prod(sizeComput);
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
onset_of_trigger = gpuArray(round(triggert/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%% stablizing noise
InoiseR1=gpuArray(X*0);
InoiseG1=gpuArray(X*0);
InoiseR2=gpuArray(X*0);
InoiseG2=gpuArray(X*0);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
R1 = gpuArray(X*initialvals(1,1)) + InoiseR1;
R2 = gpuArray(X*initialvals(1,2)) + InoiseR2;
G1 = gpuArray(X*initialvals(2,1)) + InoiseG1;
G2 = gpuArray(X*initialvals(2,2)) + InoiseG2;
%% simulation
%% 0st stage, before presenting, V = 0, alpha = 0
if onset_of_trigger == 0
    v = gpuArray(repmat(w*wasym,1,1,sims));
end
for ti = 1:(onset_of_stimuli-1)
    if ti == onset_of_trigger
        v = gpuArray(repmat(w*wasym,1,1,sims));
    end
    % update R, G, I
    G1 = G1 + (-G1 + w*R1 + v.*R2)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + v.*R1 + w*R2)/Tau2*dtArray + InoiseG2;
    R1 = R1 + (-R1)/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2)/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
end
if isempty(ti)
    ti0 = 0;
else
    ti0 = ti;
end

%% 1st stage, after presenting, before decision, V turned on, working memory (alpha) up, beta still = 0
V1Input = gpuArray(repmat(V1mat,1,1,sims));
V2Input = gpuArray(repmat(V2mat,1,1,sims));
Scale = (V1Input+V2Input);
V1Array = V1Input;
V2Array = V2Input;
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
V1Array = (V1Input + gpuArray.randn(size(V1Array)).*Scale.*sgmInput).*(V1Array ~= 0);
V2Array = (V2Input + gpuArray.randn(size(V2Array)).*Scale.*sgmInput).*(V2Array ~= 0);
R1Out = gpuArray(zeros(sizeComput));
R2Out = gpuArray(zeros(sizeComput));
rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
Continue = gpuArray(ones(sizeComput));
for ti = (ti0 + 1):max(total_time_steps(:)) %(onset_of_trigger -1)
    if stoprule == 1
        if NComput == 0
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
    if ti == onset_of_trigger
        v = gpuArray(repmat(w*wasym,1,1,sims));
    end
    if (mod(ti*dt, .05) == 0)
        V1Array = (V1Input + gpuArray.randn(size(V1Array)).*Scale.*sgmInput).*(V1Array ~= 0);
        V2Array = (V2Input + gpuArray.randn(size(V2Array)).*Scale.*sgmInput).*(V2Array ~= 0);
        % V1Array = V1Input.*(V1Array ~= 0) + gpuArray(randn(size(V1Array))*sgmInput*256/30).*(V1Array ~= 0);
        % V2Array = V2Input.*(V2Array ~= 0) + gpuArray(randn(size(V2Array))*sgmInput*256/30).*(V2Array ~= 0);
    end
    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w*R1 + v.*R2)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + v.*R1 + w*R2)/Tau2*dtArray + InoiseG2;
    R1 = R1 + (-R1 + (V1Array + alpha11.*R1)./(1+Gaba*G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + (V2Array + alpha22.*R2)./(1+Gaba*G2old))/Tau1*dtArray + InoiseR2;
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = R1 >= 0;
    R1 = R1 .* inside;
    inside = R2 >= 0;
    R2 = R2 .* inside;
    % threshold detecting
    inside = abs(R1 - R2) >= threshArray; %(R1 >= threshArray) + (R2 >= threshArray);
    flip = inside .* (rt == 0) .* (ti > onset_of_trigger);
    NComput = NComput - sum(flip(:));
    rt = rt + gpuArray(ti-onset_of_trigger).*flip*dtArray;
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    Continue = choice == 0;
    R1Out = R1Out + R1.*flip;
    R2Out = R2Out + R2.*flip;
end
R1Out = R1Out + R1.*(R1Out == 0);
R2Out = R2Out + R2.*(R2Out == 0);
%% Calculations for output
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2, 0 choice is not made
rt(rt==0) = NaN;
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
%inside = R1Out < R2Out;
%inside = inside + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
%argmaxR = inside + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
% argmaxR = 2 - mean(inside,3);
argmaxR = NaN;
dR = abs(R2Out - R1Out);