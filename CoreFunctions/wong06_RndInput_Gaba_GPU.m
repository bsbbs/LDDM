function [rt, choice, argmaxR, dR] = wong06_RndInput_Gaba_GPU(Vinput,Gaba,miu0,sgm,sgmInput,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims)
%%%%%%%%%%%%
% adapted version of the appendix function in Wong & Wang, 2006
% Created by Bo Shen, NYU, 2019
%% GPU calculation, only binary choice is allowed. N is limited as 2.
%%%%%%%%%%%%%%%%%%%%%
%%% CAUTION: because of lacking gpu memory, firing rates are not smoothed by
%%% sliding windows in this GPU version
%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%
%% set parameters
% parameter in H
a = 270; %VnC^-1
b = 108; %Hz
d = .154; %s
% gamma = .641;
unit = 1; % set as 1 if the unit of time related variables is second.
% change it to .001 if the unit of time related varibales is msec.

% connectivity strength
% JN = [.2609, -.0497;
%     -.0497, .2609]; %nA
JAext = 5.2 * 10^-4; %nA/Hz

%% preparation
%GabaArray = gpuArray(Gaba);
miu0Array = gpuArray(miu0);
I0Array = gpuArray(I0);
sigmaArray = gpuArray(sgm);
tauSArray = gpuArray(tauS/unit);
tauAMPAArray = gpuArray(tauAMPA/unit);
dtArray = gpuArray(dt/unit);
JN11 = (JN(1,1));
JN12 = (JN(1,2));
JN21 = (JN(2,1));
JN22 = (JN(2,2));
JAextArray = gpuArray(JAext);
threshArray = gpuArray(thresh);
if isstruct(Vinput)
    name = fieldnames(Vinput);
    V1mat = Vinput.(name{1});
    V2mat = Vinput.(name{2});
else
    V1mat = Vinput(:,1);
    V2mat = Vinput(:,2);
end
cp1Array = gpuArray(repmat(V1mat,1,1,sims)/256);
cp2Array = gpuArray(repmat(V2mat,1,1,sims)/256);
sizeVinput = size(V1mat);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);

presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
% sliding window
time_wind = .050/dt/unit;  % Temporal window size for averaging, 50 msec
slide_wind = .005/dt/unit;  % Sliding step for window, 5 msec
% time line
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%% get initial values
H1i = initialvals(1,1);
H2i = initialvals(1,2);
S1i = initialvals(2,1);
S2i = initialvals(2,2);
syms x;
if H1i < 10^-5
    x1 = 0;
else
    s1 = solve((a*x - b)./(1 - exp(-d*(a*x - b))) == H1i, x);
    x1 = (double(vpa(s1)));
end
if H2i < 10^-5
    x2 = 0;
else
    s2 = solve((a*x - b)./(1 - exp(-d*(a*x - b))) == H2i, x);
    x2 = (double(vpa(s2)));
end
%% stablizing noise
Inoise1 = gpuArray(zeros(sizeComput));
Inoise2 = gpuArray(zeros(sizeComput));
stablizetime = round(.2/unit/dt);
for kk = 1:stablizetime
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    H1 = real((a*(x1 + Inoise1) - b)./(1 - exp(-d*(a*(x1 + Inoise1) - b))));
    H2 = real((a*(x2 + Inoise2) - b)./(1 - exp(-d*(a*(x2 + Inoise2) - b))));
    inside = H1 >= 0;
    H1 = H1 .* inside;
    inside = H2 >= 0;
    H2 = H2 .* inside;
    % S1 = S1i + (1-S1i).*(H1-H1i)*unit*gamma*dt;
    % S2 = S2i + (1-S2i).*(H2-H2i)*unit*gamma*dt;
    loci = mod(kk,time_wind) + (mod(kk,time_wind) == 0)*time_wind;
    H1buffer(:,:,:,loci) = H1;
    H2buffer(:,:,:,loci) = H2;
end
H1 = real(mean(H1buffer,4));
H2 = real(mean(H2buffer,4));
S1 = H1*gamma*tauS./(H1*gamma*tauS+1);
S2 = H2*gamma*tauS./(H2*gamma*tauS+1);
%% simulation begin
%% 0st stage, before presenting, V = 0
for ti = 1:(onset_of_stimuli-1)
    % update nerual firing rates and synaptic activities
    x1 = JN11*S1 + Gaba*JN12*S2 + I0 + Inoise1;
    x2 = Gaba*JN21*S1 + JN22*S2 + I0 + Inoise2;
    H1old = H1;
    H2old = H2;
    H1 = real((a*x1 - b)./(1 - exp(-d*(a*x1 - b))));
    H2 = real((a*x2 - b)./(1 - exp(-d*(a*x2 - b))));
    S1 = S1 + (-S1/tauS + (1-S1).*H1old*unit*gamma)*dt;
    S2 = S2 + (-S2/tauS + (1-S2).*H2old*unit*gamma)*dt;
    % update noise
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    % To ensure firing rates are always positive (noise may cause negative)
    inside = H1 >= 0;
    H1 = H1 .* inside;
    inside = H2 >= 0;
    H2 = H2 .* inside;
    loci = mod(kk+ti,time_wind) + (mod(kk+ti,time_wind) == 0)*time_wind;
    H1buffer(:,:,:,loci) = H1;
    H2buffer(:,:,:,loci) = H2;
end
if isempty(ti)
    ti0 = 0;
else
    ti0 = ti;
end

%% 1st stage, stimuli turned on.
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
I1 = JAext*miu0*(cp1Array + gpuArray(randn(size(cp1Array))*sgmInput));
I2 = JAext*miu0*(cp2Array + gpuArray(randn(size(cp2Array))*sgmInput));

rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
R1 = mean(H1buffer,4);
R2 = mean(H2buffer,4);
R1Out = gpuArray(zeros(sizeComput));
R2Out = gpuArray(zeros(sizeComput));
for ti = (ti0 + 1):(max(total_time_steps(:)))
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
    if stoprule == 1
        if NComput == 0
            break;
        end
    end
    % update nerual firing rates and synaptic activities
    if mod(ti*dt,.05) == 0
        I1 = JAext*miu0*(cp1Array + gpuArray(randn(size(cp1Array))*sgmInput));
        I2 = JAext*miu0*(cp2Array + gpuArray(randn(size(cp2Array))*sgmInput));
    end
    x1 = JN11*S1 + Gaba*JN12*S2 + I0 + I1 + Inoise1;
    x2 = Gaba*JN21*S1 + JN22*S2 + I0 + I2 + Inoise2;
    H1old = H1;
    H2old = H2;
    H1 = real((a*x1 - b)./(1 - exp(-d*(a*x1 - b))));
    H2 = real((a*x2 - b)./(1 - exp(-d*(a*x2 - b))));
    S1 = S1 + (-S1/tauS + (1-S1).*H1old*unit*gamma)*dt;
    S2 = S2 + (-S2/tauS + (1-S2).*H2old*unit*gamma)*dt;
    % update noise
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    % To ensure firing rates are always positive (noise may cause negative)
    inside = H1 >= 0;
    H1 = H1 .* inside;
    inside = H2 >= 0;
    H2 = H2 .* inside;
    loci = mod(kk+ti,time_wind) + (mod(kk+ti,time_wind) == 0)*time_wind;
    H1buffer(:,:,:,loci) = H1;
    H2buffer(:,:,:,loci) = H2;
    % threshold detecting
    R1 = mean(H1buffer,4);
    R2 = mean(H2buffer,4);
    inside = (R1 >= thresh) + (R2 >= thresh);
    flip = (inside > 0) .* (rt == 0);
    NComput = NComput - sum(flip(:));
    rt = rt + (ti-ti0)*dt*flip;
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    R1Out = R1Out + R1.*flip;
    R2Out = R2Out + R2.*flip;
end
R1Out = R1Out + R1.*(R1Out == 0);
R2Out = R2Out + R2.*(R2Out == 0);
%% Calculations for output
choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2, 0 choice is not made
rt(rt==0) = NaN;
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
inside = R1Out < R2Out;
inside = inside + (R1Out == R2Out)/2;% R1 < R2 recorded as 1, R1 > R2 recorded as 0, R1 == R2 recorded as .5;
argmaxR = inside + 1; % keep consistant as other models, 1 for choose R1, 2 for choose R2, and 1.5 for equal
% argmaxR = 2 - mean(inside,3);
dR = abs(R2Out - R1Out);