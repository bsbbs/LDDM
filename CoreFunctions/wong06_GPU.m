function [choice, rt] = wong06_GPU(cp,miu0,sgm,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims)
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

%% preparation
I0Array = gpuArray(I0);
JN11 = (JN(1,1));
JN12 = (JN(1,2));
JN21 = (JN(2,1));
JN22 = (JN(2,2));
if isstruct(cp)
    name = fieldnames(cp);
    c1 = cp.(name{1});
    c2 = cp.(name{2});
else
    c1 = cp(:,1);
    c2 = cp(:,2);
end
cp1Array = gpuArray(repmat(c1,1,1,sims)/256);
cp2Array = gpuArray(repmat(c2,1,1,sims)/256);
sizeVinput = size(c1);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);

presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
% sliding window
time_wind = .050/dt/unit;  % Temporal window size for averaging, 50 msec
% not used, slide_wind = .005/dt/unit;  % Sliding step for window, 5 msec
% time line
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%% get initial values
H1new = gpuArray(ones(sizeComput)*initialvals(1,1));
H2new = gpuArray(ones(sizeComput)*initialvals(1,2));
S1 = gpuArray(ones(sizeComput)*initialvals(2,1));
S2 = gpuArray(ones(sizeComput)*initialvals(2,2));
% noise
Inoise1 = gpuArray.randn(sizeComput).*sgm;
Inoise2 = gpuArray.randn(sizeComput).*sgm;
%% simulation begin
rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
c1Array = gpuArray(zeros(sizeComput)); % input values set as zeros before onset of stimuli
c2Array = gpuArray(zeros(sizeComput));
for ti = 1:(max(total_time_steps(:)))
    if ti >= onset_of_stimuli
        cp1Array = gpuArray(repmat(cp1,1,1,sims));
        cp2Array = gpuArray(repmat(cp2,1,1,sims));
    end
    if numel(offset_of_stimuli) == 1
        if ti >= offset_of_stimuli
            c1Array = gpuArray(zeros(sizeComput));
            c2Array = gpuArray(zeros(sizeComput));
        end
    elseif prod(size(offset_of_stimuli) == size(cp1))
        offset_of_stimulimat = repmat(offset_of_stimuli,1,1,sims);
        withdraw = offset_of_stimulimat >= ti;
        c1Array(withdraw) = 0;
        c2Array(withdraw) = 0;
    end
    % update nerual firing rates H
    H1 = H1new;
    H2 = H2new;
    I1 = JAext*miu0*cp1Array;
    I2 = JAext*miu0*cp2Array;
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
    inside = (R1 >= thresh) + (R2 >= thresh);
    flip = (inside > 0) .* (rt == 0);
    NComput = NComput - sum(flip(:));
    rt = rt + (ti-onset_of_stimuli)*dt*flip;
    choice = choice + ((R2 > R1) - (R1 > R2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    % when all of the channel hit the decision boundary, stop simulation
    if stoprule == 1
        if NComput == 0
            break;
        end
    end
end
choice((choice==0)&(rt==0)) = NaN; % 1 for choosing H1; 2 for choosing H2; 3 for choosing H3; 0 for having equal high values; NaN choice is not made
rt(rt==0) = NaN;
rt = gather(rt);
choice = gather(choice);

choice = choice/2; %1 choose R1, 2 choose R2, 1.5 R1 = R2, 0 choice is not made
rt(rt==0) = NaN;
choice(choice == 0) = NaN; %1 choose R1, 2 choose R2, 1.5 R1 = R2, NaN choice is not made
rt = gather(rt);
choice = gather(choice);