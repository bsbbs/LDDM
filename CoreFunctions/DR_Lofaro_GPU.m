function [argmaxR] = DR_Lofaro_GPU(Vinput, w, sgm, Tau, presentt, dur, dt, initialvals, stimdur, sims)
%%%%%%%%%%%%%%%%%%%
%% preparation
sgmArray = gpuArray(sgm);
tauN = gpuArray(0.002);% time constant for Ornstein-Uhlenbeck process of noise
dtArray = gpuArray(dt);
w11 = gpuArray(w(1,1));
w12 = gpuArray(w(1,2));
w21 = gpuArray(w(2,1));
w22 = gpuArray(w(2,2));
Tau1 = gpuArray(Tau(1));
Tau2 = gpuArray(Tau(2));
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
V1Array = gpuArray(zeros(sizeComput));
V2Array = gpuArray(zeros(sizeComput));
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%% simulation
X = ones(sizeComput);
R1 = gpuArray(X*initialvals(1,1));
R2 = gpuArray(X*initialvals(1,2));
G1 = gpuArray(X*initialvals(2,1));
G2 = gpuArray(X*initialvals(2,2));
%% stablizing noise
InoiseR1=gpuArray(X*0);
InoiseG1=gpuArray(X*0);
InoiseR2=gpuArray(X*0);
InoiseG2=gpuArray(X*0);
stablizetime = round(1/dt);
for kk = 1:stablizetime
    % update noise
    InoiseR1 = InoiseR1 + (-InoiseR1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseR2 = InoiseR2 + (-InoiseR2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG1 = InoiseG1 + (-InoiseG1 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
    InoiseG2 = InoiseG2 + (-InoiseG2 + gpuArray.randn(sizeComput)*sqrt(dtArray)*sgmArray)/tauN*dtArray;
end
R1 = R1 + InoiseR1;
R2 = R2 + InoiseR2;
G1 = G1 + InoiseG1;
G2 = G2 + InoiseG2;
%% 0st stage, before presenting, V = 0, alpha = 0
for ti = 1:(onset_of_stimuli-1)
    % update R, G, I
    G1 = G1 + (-G1 + w11*R1 + w12*R2)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2)/Tau2*dtArray + InoiseG2;
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
for ti = (ti0 + 1):max(total_time_steps(:))
    if numel(offset_of_stimuli) == 1
        if ti == offset_of_stimuli
            V1Array = gpuArray(zeros(sizeComput));
            V2Array = gpuArray(zeros(sizeComput));
        end
    elseif prod(size(offset_of_stimuli) == size(V1mat))
        flip = offset_of_stimulimat == ti;
        V1Array(flip) = 0;
        V2Array(flip) = 0;
        R1Out(flip) = R1(flip);
        R2Out(flip) = R2(flip);
    end
    % update R, G, I
    G1old = G1; G2old = G2;
    G1 = G1 + (-G1 + w11*R1 + w12*R2)/Tau2*dtArray + InoiseG1;
    G2 = G2 + (-G2 + w21*R1 + w22*R2)/Tau2*dtArray + InoiseG2;
    R1 = R1 + (-R1 + V1Array./(1+G1old))/Tau1*dtArray + InoiseR1;
    R2 = R2 + (-R2 + V2Array./(1+G2old))/Tau1*dtArray + InoiseR2;
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
choice = R1Out < R2Out;
choice = choice + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
choice = choice + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
argmaxR = 2 - mean(choice,3);
