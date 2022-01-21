function [choice, rt] = LDDM_GPU3(cp, scale, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule, sims)
%%%%%%%%%%%%%%%%%%
%% GPU calculation, only trinary choice is allowed. N is limited as 3.
% Created by Bo Shen, NYU, 2021
% Vinput, w, a (alpha), b (beta), Tau, dt are all in the same meaning as
% the names of parameters in the paper.
% - cp: input coherence vectors/matrices, allows two formats.
%   Format A: input values as a MxN array. N is the number of choice items,
%       M is the number of input sets
%   - example
%   cp =  [0.6250    0.3750 0.3750
%         0.5625    0.4375  0.3750
%         0.5273    0.4688  0.3750];
%    Format B: input coherences as a structure with fields defining c1, c2, c3 matrices.
%       The input matrices have to be the same dimension of 2. Values at the same position
%       of the matrices will be used as a set of inputs.
%    - example
%    c = [0, .032, .064, .128, .256, .512]';
%    c1 = [1 - flip(c); 1 + c];
%    c2 = ones(size(c1));
%    c3 = c;
%    [cp.cp1, ~] = meshgrid(c1, c3);
%    [cp.cp2, cp.cp3] = meshgrid(c2, c3);
% - scale: input scaling parameter.
%    Vinput.V1 = cp1*scale; Vinput.V2 = cp2*scale; Vinput.V3 = cp3*scale;
% - w: connection weight from R neurons to G neurons. w must be a NxN matrix,
%   describing connections between each pair of R and G.
%   - example    [w11, w12, w13    [1, 1, 1
%                 w21, w22, w23 =  1, 1, 1
%                 w31, w32, w33]    1, 1, 1]
% - a: self-excitation parameter, describing the weight of
%   self-excitation for R neurons. alpha must be a NxN diagonal matrix.
%   Non-zero values only on the diagnal means that each R neuron has its own
%   excitation loop. Lateral or cross excitation does not exist.
%   - example      [a11, a12, a13;   [15, 0, 0
%                 a21, a22, a23 =  0, 15, 0
%                   a31, a32, a33]  0, 0, 15]
% - b: the weight of local disinhibition,
%   - example      [b11, b12, b13   [2, 0, 0
%                 b21, b22, b23 =  0, 2, 0
%                   b31, b32, b33]  0, 0, 2]
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
%   - example      initialvals = [R1_0, R2_0, R3_0
%                                 G1_0, G2_0, G3_0]
%                                 I1_0, I2_0, I3_0]
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
w13 = gpuArray(w(1,3));
w21 = gpuArray(w(2,1));
w22 = gpuArray(w(2,2));
w23 = gpuArray(w(2,3));
w31 = gpuArray(w(3,1));
w32 = gpuArray(w(3,2));
w33 = gpuArray(w(3,3));
TauR = gpuArray(Tau(1));
TauG = gpuArray(Tau(2));
TauI = gpuArray(Tau(3));
alpha11 = gpuArray(a(1,1));
alpha12 = gpuArray(a(1,2));
alpha13 = gpuArray(a(1,3));
alpha21 = gpuArray(a(2,1));
alpha22 = gpuArray(a(2,2));
alpha23 = gpuArray(a(2,3));
alpha31 = gpuArray(a(3,1));
alpha32 = gpuArray(a(3,2));
alpha33 = gpuArray(a(3,3));
beta11 = gpuArray(b(1,1));
beta12 = gpuArray(b(1,2));
beta13 = gpuArray(b(1,3));
beta21 = gpuArray(b(2,1));
beta22 = gpuArray(b(2,2));
beta23 = gpuArray(b(2,3));
beta31 = gpuArray(b(3,1));
beta32 = gpuArray(b(3,2));
beta33 = gpuArray(b(3,3));
threshArray = gpuArray(thresh);
if isstruct(cp)
    name = fieldnames(cp);
    V1mat = cp.(name{1}).*scale;
    V2mat = cp.(name{2}).*scale;
    V3mat = cp.(name{3}).*scale;
else
    
    V1mat = cp(:,1).*scale;
    V2mat = cp(:,2).*scale;
    V3mat = cp(:,3).*scale;
end
sizeVinput = size(V1mat);
sizeComput = [sizeVinput, sims];
X = ones(sizeComput);
NComput = prod(sizeComput);
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
onset_of_trigger = gpuArray(round(triggert/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%%
trigismat = 0;
if prod(size(onset_of_trigger) == size(V1mat))
    onset_of_trigger = repmat(onset_of_trigger,1,1,sims);
    trigismat = 1;
end
%% stablizing noise
InoiseV1=gpuArray(X*0);
InoiseV2=gpuArray(X*0);
InoiseV3=gpuArray(X*0);
InoiseR1=gpuArray(X*0);
InoiseG1=gpuArray(X*0);
InoiseI1=gpuArray(X*0);
InoiseR2=gpuArray(X*0);
InoiseG2=gpuArray(X*0);
InoiseI2=gpuArray(X*0);
InoiseR3=gpuArray(X*0);
InoiseG3=gpuArray(X*0);
InoiseI3=gpuArray(X*0);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    InoiseV1 = InoiseV1 + (-InoiseV1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseV2 = InoiseV2 + (-InoiseV2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseV3 = InoiseV3 + (-InoiseV3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR3 = InoiseR3 + (-InoiseR3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG3 = InoiseG3 + (-InoiseG3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI3 = InoiseI3 + (-InoiseI3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
R1new = gpuArray(X*initialvals(1,1).*scale) + InoiseR1;
R2new = gpuArray(X*initialvals(1,2).*scale) + InoiseR2;
R3new = gpuArray(X*initialvals(1,3).*scale) + InoiseR3;
G1 = gpuArray(X*initialvals(2,1).*scale) + InoiseG1;
G2 = gpuArray(X*initialvals(2,2).*scale) + InoiseG2;
G3 = gpuArray(X*initialvals(2,3).*scale) + InoiseG3;
I1 = gpuArray(X*initialvals(3,1).*scale) + InoiseI1;
I2 = gpuArray(X*initialvals(3,2).*scale) + InoiseI2;
I3 = gpuArray(X*initialvals(3,3).*scale) + InoiseI3;
%% simulation
rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
Continue = gpuArray(ones(sizeComput));
BetasUp = gpuArray(zeros(sizeComput));
for ti = 0:max(total_time_steps(:))
    if ti >= onset_of_stimuli
        V1Array = gpuArray(repmat(V1mat,1,1,sims));
        V2Array = gpuArray(repmat(V2mat,1,1,sims));
        V3Array = gpuArray(repmat(V3mat,1,1,sims));
    end
    if numel(offset_of_stimuli) == 1
        if ti >= offset_of_stimuli
            V1Array = gpuArray(zeros(sizeComput));
            V2Array = gpuArray(zeros(sizeComput));
            V3Array = gpuArray(zeros(sizeComput));
        end
    elseif prod(size(offset_of_stimuli) == size(V1mat))
        offset_of_stimulimat = repmat(offset_of_stimuli,1,1,sims);
        withdraw = offset_of_stimulimat >= ti;
        V1Array(withdraw) = 0;
        V2Array(withdraw) = 0;
        V3Array(withdraw) = 0;
    end
    if ~trigismat % numel(onset_of_trigger) == 1
        if ti >= onset_of_trigger
            BetasUp = gpuArray(ones(sizeComput));
        end
    elseif trigismat %prod(size(onset_of_trigger) == size(V1mat))
        BetasUp(ti >= onset_of_trigger) = 1;
    end
%     V1Array = V1Array + InoiseV1;
%     V2Array = V2Array + InoiseV2;
%     V3Array = V3Array + InoiseV3;
    % update R, G, I
    R1 = R1new;
    R2 = R2new;
    R3 = R3new;
    R1new = R1 + (-R1 + Continue.*(V1Array + alpha11*R1+alpha12*R2+alpha13*R3)./(1+G1))/TauR*dtArray + InoiseR1;
    R2new = R2 + (-R2 + Continue.*(V2Array + alpha21*R1+alpha22*R2+alpha23*R3)./(1+G2))/TauR*dtArray + InoiseR2;
    R3new = R3 + (-R3 + Continue.*(V3Array + alpha31*R1+alpha32*R2+alpha33*R3)./(1+G3))/TauR*dtArray + InoiseR3;
    G1 = G1 + (-G1 + w11*R1 + w12*R2 + w13*R3  - I1)/TauG*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2 + w23*R3  - I2)/TauG*dtArray + InoiseG2;
    G3 = G3 + (-G3 + w31*R1 + w32*R2 + w33*R3  - I3)/TauG*dtArray + InoiseG3;
    I1 = I1 + (-I1 + (beta11*R1 + beta12*R2 + beta13*R3).*Continue.*BetasUp)/TauI*dtArray + InoiseI1;
    I2 = I2 + (-I2 + (beta21*R1 + beta22*R2 + beta23*R3).*Continue.*BetasUp)/TauI*dtArray + InoiseI2;
    I3 = I3 + (-I3 + (beta31*R1 + beta32*R2 + beta33*R3).*Continue.*BetasUp)/TauI*dtArray + InoiseI3;
    % update noise
    InoiseV1 = InoiseV1 + (-InoiseV1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseV2 = InoiseV2 + (-InoiseV2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseV3 = InoiseV3 + (-InoiseV3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR3 = InoiseR3 + (-InoiseR3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG3 = InoiseG3 + (-InoiseG3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI1 = InoiseI1 + (-InoiseI1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI2 = InoiseI2 + (-InoiseI2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseI3 = InoiseI3 + (-InoiseI3 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = R1new >= 0;
    R1new = R1new .* inside;
    inside = R2new >= 0;
    R2new = R2new .* inside;
    inside = R3new >= 0;
    R3new = R3new .* inside;
    inside = G1 >= 0;
    G1 = G1 .* inside;
    inside = G2 >= 0;
    G2 = G2 .* inside;
    inside = G3 >= 0;
    G3 = G3 .* inside;
    inside = I1 >= 0;
    I1 = I1 .* inside;
    inside = I2 >= 0;
    I2 = I2 .* inside;
    inside = I3 >= 0;
    I3 = I3 .* inside;
    % threshold detecting
    inside = (R1 >= threshArray) + (R2 >= threshArray) + (R3 >= threshArray);
    flip = (inside > 0) .* (rt == 0) .* (ti > onset_of_trigger);
    NComput = NComput - sum(flip(:));
    rt = rt + (ti-onset_of_trigger)*dtArray.*flip;
    choice = choice + ((R1 > R2).*(R1 > R3) + 2*(R2 > R1).*(R2 > R3) + 3*(R3 > R1).*(R3 > R2)) .* flip;
    % 1 for choosing R1; 2 for choosing R2; 3 for choosing R3, 0 for having high equal values, choice is not made
    Continue = rt == 0;
    % after every channel hit the decision threshold, stop simulation
    if stoprule == 1
        if NComput == 0
            break;
        end
    end
end
choice(rt == 0) = NaN; % time out trials
rt(rt==0) = NaN;
choice = gather(choice);
rt = gather(rt);