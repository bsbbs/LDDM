function [choice, rt] = wong06_GPU3(cp,miu0,sgm,I0,JN,...
    gamma, tauS, tauAMPA, dur, dt, presentt, stimdur, thresh, initialvals, stoprule, sims)
%%%%%%%%%%%%
% adapted version of the appendix function in Wong & Wang, 2006
% Created by Bo Shen, NYU, 2019
% GPU calculation for trinary choice. The number of inputs is limited as 3.
%%%%%%%%%%%%%%%%%%%%%
%%% CAUTION: because of lacking gpu memory, firing rates are not exported as time series
%%%%%%%%%%%%%%%%%%%%%
% - cp as input coherence, follow the format listed below
%   Format A: input coherences as a MxN array. N is the number of choice items,
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
%%%%%%%%%%%%%%%
%% set parameters
% parameter in H
a = 270; %VnC^-1
b = 108; %Hz
d = .154; %s
% time unit parameter
unit = 1; % set as 1 if the unit of time related variables is second.
% change it to .001 if the unit of time related varibales is msec.
% connectivity strength paramter
JAext = 5.2 * 10^-4; %nA/Hz

%% preparation
I0Array = gpuArray(I0);
JN11 = (JN(1,1));
JN12 = (JN(1,2));
JN13 = (JN(1,3));
JN21 = (JN(2,1));
JN22 = (JN(2,2));
JN23 = (JN(2,3));
JN31 = (JN(3,1));
JN32 = (JN(3,2));
JN33 = (JN(3,3));
if isstruct(cp) % Format B
    name = fieldnames(cp);
    cp1 = cp.(name{1});
    cp2 = cp.(name{2});
    cp3 = cp.(name{2});
else % Formart A
    cp1 = cp(:,1);
    cp2 = cp(:,2);
    cp3 = cp(:,3);
end
size_input = size(cp1);
sizeComput = [size_input, sims];
NComput = prod(sizeComput);
% sliding window
time_wind = .050/dt/unit;  % Temporal window size for averaging, 50 msec
% slide_wind = .005/dt/unit;  % Sliding window for output, 5 msec, discarded to use dt instead
% time line
presentt = presentt/unit;
dur = dur/unit;
stimdur = stimdur/unit;
total_time_steps = gpuArray(round(dur/dt));
onset_of_stimuli = gpuArray(round(presentt/dt));
stim_duration = gpuArray(round(stimdur/dt));
offset_of_stimuli = onset_of_stimuli + stim_duration;
%% get initial values
H1new = gpuArray(ones(sizeComput)*initialvals(1,1));
H2new = gpuArray(ones(sizeComput)*initialvals(1,2));
H3new = gpuArray(ones(sizeComput)*initialvals(1,3));
S1 = gpuArray(ones(sizeComput)*initialvals(2,1));
S2 = gpuArray(ones(sizeComput)*initialvals(2,2));
S3 = gpuArray(ones(sizeComput)*initialvals(2,3));
% noise
Inoise1 = gpuArray.randn(sizeComput).*sgm;
Inoise2 = gpuArray.randn(sizeComput).*sgm;
Inoise3 = gpuArray.randn(sizeComput).*sgm;
%% simulation begin
rt = gpuArray(zeros(sizeComput));
choice = gpuArray(zeros(sizeComput));
c1Array = gpuArray(zeros(sizeComput)); % input values set as zeros before onset of stimuli
c2Array = gpuArray(zeros(sizeComput));
c3Array = gpuArray(zeros(sizeComput));
for ti = 1:(max(total_time_steps(:)))
    if ti >= onset_of_stimuli
        cp1Array = gpuArray(repmat(cp1,1,1,sims));
        cp2Array = gpuArray(repmat(cp2,1,1,sims));
        cp3Array = gpuArray(repmat(cp3,1,1,sims));
    end
    if numel(offset_of_stimuli) == 1
        if ti >= offset_of_stimuli
            c1Array = gpuArray(zeros(sizeComput));
            c2Array = gpuArray(zeros(sizeComput));
            c3Array = gpuArray(zeros(sizeComput));
        end
    elseif prod(size(offset_of_stimuli) == size(cp1))
        offset_of_stimulimat = repmat(offset_of_stimuli,1,1,sims);
        withdraw = offset_of_stimulimat >= ti;
        c1Array(withdraw) = 0;
        c2Array(withdraw) = 0;
        c3Array(withdraw) = 0;
    end
    % update nerual firing rates H
    H1 = H1new;
    H2 = H2new;
    H3 = H3new;
    I1 = JAext*miu0*cp1Array;
    I2 = JAext*miu0*cp2Array;
    I3 = JAext*miu0*cp3Array;
    x1 = JN11*S1 + JN12*S2 + JN13*S3 + I0Array + I1 + Inoise1;
    x2 = JN21*S1 + JN22*S2 + JN23*S3 + I0Array + I2 + Inoise2;
    x3 = JN31*S1 + JN32*S2 + JN33*S3 + I0Array + I3 + Inoise3;
    H1new = (a*x1 - b)./(1 - exp(-d*(a*x1 - b)));
    H2new = (a*x2 - b)./(1 - exp(-d*(a*x2 - b)));
    H3new = (a*x3 - b)./(1 - exp(-d*(a*x3 - b)));
    % update synaptic activities S
    S1 = S1 + (-S1/tauS + (1-S1).*H1*unit*gamma)*dt;
    S2 = S2 + (-S2/tauS + (1-S2).*H2*unit*gamma)*dt;
    S3 = S3 + (-S3/tauS + (1-S3).*H3*unit*gamma)*dt;
    % update noise
    Inoise1 = Inoise1 +(-Inoise1/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise2 = Inoise2 +(-Inoise2/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    Inoise3 = Inoise3 +(-Inoise3/tauAMPA*dt + gpuArray.randn(sizeComput).*sqrt(dt/tauAMPA*sgm^2));
    % To ensure firing rates are always positive (noise may cause negative)
    inside = H1new >= 0;
    H1new = H1new .* inside;
    inside = H2new >= 0;
    H2new = H2new .* inside;
    inside = H3new >= 0;
    H3new = H3new .* inside;
    loci = mod(ti,time_wind) + (mod(ti,time_wind) == 0)*time_wind;
    H1buffer(:,:,:,loci) = H1;
    H2buffer(:,:,:,loci) = H2;
    H3buffer(:,:,:,loci) = H3;
    H1_wind = mean(H1buffer,4); % smoothed firing rates by sliding windows
    H2_wind = mean(H2buffer,4);
    H3_wind = mean(H3buffer,4);
    % threshold detecting
    inside = (H1_wind >= thresh) + (H2_wind >= thresh) + (H3_wind >= thresh); 
    flip = (inside > 0) .* (rt == 0);
    NComput = NComput - sum(flip(:));
    rt = rt + (ti-onset_of_stimuli)*dt*flip;
    choice = choice + ( (H1_wind > H2_wind).*(H1_wind > H3_wind) + 2*(H2_wind > H1_wind).*(H2_wind > H3_wind) + 3*(H3_wind > H1_wind).*(H3_wind > H2_wind) ) .* flip;
    % 1 for choosing H1; 2 for choosing H2; 3 for choosing H3; 0 for having equal high values, choice is not made
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