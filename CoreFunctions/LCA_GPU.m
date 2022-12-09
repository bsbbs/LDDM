function [rt, choice, argmaxR] = LCA_GPU(RhoMat, k, beta, sgm, T0, thresh, dur, dt, sims)
sgmArray = gpuArray(sgm);
tau = gpuArray(.1);
dtArray = gpuArray(dt);
dt_tau = dt/tau;
k11 = gpuArray(k(1,1));
k12 = gpuArray(k(1,2));
k21 = gpuArray(k(2,1));
k22 = gpuArray(k(2,2));
beta11 = gpuArray(beta(1,1));
beta12 = gpuArray(beta(1,2));
beta21 = gpuArray(beta(2,1));
beta22 = gpuArray(beta(2,2));
threshArray = gpuArray(thresh);
if isstruct(RhoMat)
    name = fieldnames(RhoMat);
    Rho1mat = RhoMat.(name{1});
    Rho2mat = RhoMat.(name{2});
else
    Rho1mat = RhoMat(:,1);
    Rho2mat = RhoMat(:,2);
end
Rho1inputArray = gpuArray(repmat(Rho1mat,1,1,sims));
Rho2inputArray = gpuArray(repmat(Rho2mat,1,1,sims));
sizeVinput = size(Rho1mat);
sizeComput = [sizeVinput, sims];
NComput = prod(sizeComput);

%% initialize variables
rt = gpuArray.zeros(sizeComput);
choice = gpuArray.zeros(sizeComput);
x1Out = gpuArray.nan(sizeComput); % records the values at decision, or the end of simulation if choice was not made
x2Out = gpuArray.nan(sizeComput);
x1 = gpuArray.zeros(sizeComput);
x2 = gpuArray.zeros(sizeComput);
% Continue = gpuArray(ones(sizeComput)); % to mark the trials that choices haven't made yet
%% simulation: premotion stage (ti < 0), dot motion task begin at ti = 0
% ti = 0 sorted at the the beginning of evidence accumulation, T0 capture the
% non-decision delay in the beginning and the execution time at the end
for ti = 1:(dur/dt)
    if NComput == 0
        break;
    end
    % update x
    x1 = x1 + (Rho1inputArray - k11*x1 - beta12*x2)*dt_tau + gpuArray.randn(sizeComput)*sgmArray*sqrt(dt_tau);
    x2 = x2 + (Rho1inputArray - k22*x2 - beta21*x1)*dt_tau + gpuArray.randn(sizeComput)*sgmArray*sqrt(dt_tau);
    
    % setting lower boundary, forcing neural firing rates to be non-negative
    inside = x1 >= 0;
    x1 = x1 .* inside;
    inside = x2 >= 0;
    x2 = x2 .* inside;
    
    % threshold detecting
    inside = (x1 >= threshArray) + (x2 >= threshArray);
    flip = (inside > 0) & (choice == 0);
    NComput = NComput - sum(flip(:));
    rt = rt + ti.*flip*dtArray;
    choice = choice + ((x2 > x1) - (x1 > x2) +3) .* flip; % 2 choose R1, 4 choose R2, 3 R1 = R2, 0 choice is not made
    % Continue = choice == 0;
    x1Out(flip) = x1(flip); % update the values at choice, keep others as nan
    x2Out(flip) = x2(flip);
end
% for those trials that choices were not made
x1Out(choice == 0) = x1(choice == 0);
x2Out(choice == 0) = x2(choice == 0);
%% calculate
choice = choice/2; %1 choose x1, 2 choose x2, 1.5 x1 = x2, 0 choice is not made
choice(choice == 0) = NaN; %1 choose x1, 2 choose x2, 1.5 x1 = x2, NaN choice is not made
rt(rt==0) = NaN;
rt = rt + T0;
inside = x1Out < x2Out;
inside = inside + (x1Out == x2Out)/2;% x1 < x2 recorded as 1, x1 > x2 recorded as 0, x1 == x2 recorded as .5;
argmaxR = inside + 1; % keep consistant as other models, 1 for choose x1, 2 for choose x2, and 1.5 for equal
rt = gather(rt);
choice = gather(choice);
argmaxR = gather(argmaxR);