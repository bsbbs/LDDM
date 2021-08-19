function [choice, rt, R, G, I] = LDDM(Vinput, w, a, b, sgm, Tau, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule)
%%%%%%%%%%%%%%%%%%
% The core function of local disinhibition decision model (LDDM)
% Created by Bo Shen, NYU, 2019
% Vinput, w, a(alpha), b(beta), Tau, dt are all in the same meaning as
% the names of parameters in the paper.
% - Vinput: the input value as a 1xN array. N is the number of choice items,
%   allowed value from 1 to any positive integers
%   - example      Vinput =  [316, 196];
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
%%%%%%%%%%%%%%%%%%%
tauN =0.002; % time constant for Ornstein-Uhlenbeck process of noise
%% preparation
total_time_steps = round(dur/dt);
onset_of_stimuli = round(presentt/dt);
onset_of_trigger = round(triggert/dt);
stim_duration = round(stimdur/dt);
sizeVinput = size(Vinput);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end
rt = NaN;
choice = NaN;
%% simulation begin
G(1,:) = initialvals(2,:);
R(1,:) = initialvals(1,:);
I(1,:) = initialvals(3,:);
InoiseG = randn(sizeVinput)*sgm;
InoiseR = randn(sizeVinput)*sgm;
InoiseI = randn(sizeVinput)*sgm;
for ti = 1:total_time_steps
    % update R, G, I
    dR = (-R(ti,:)' + (Vinput'*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration) + a*(ti >= onset_of_stimuli)*R(ti,:)')./(1+G(ti,:)'))/Tau(1)*dt;
    dG = (-G(ti,:)' + w * R(ti,:)' - I(ti,:)')/Tau(2)*dt;
    dI = (-I(ti,:)' + b*(ti >= onset_of_trigger)*R(ti,:)')/Tau(3)*dt;
    R(ti+1,:) = R(ti,:) + dR' + InoiseR(ti,:);
    G(ti+1,:) = G(ti,:) + dG' + InoiseG(ti,:);
    I(ti+1,:) = I(ti,:) + dI' + InoiseI(ti,:);
    
    % update noise
    dInoiseG = (-InoiseG(ti,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    dInoiseI = (-InoiseI(ti,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    dInoiseR = (-InoiseR(ti,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    InoiseG(ti+1,:) = InoiseG(ti,:) + dInoiseG;
    InoiseI(ti+1,:) = InoiseI(ti,:) + dInoiseI;
    InoiseR(ti+1,:) = InoiseR(ti,:) + dInoiseR;
    % setting lower boundary, forcing neural firing rates to be non-negative
    G(ti+1,G(ti+1,:) < 0) = 0;
    I(ti+1,I(ti+1,:) < 0) = 0;
    R(ti+1,R(ti+1,:) < 0) = 0;
    % threshold detecting
    if ti > onset_of_trigger && isnan(rt)
        if max(R(ti,:)) >= thresh
            rt = (ti-onset_of_trigger)*dt;
            choice = find(R(ti,:) == max(R(ti,:)));
            if stoprule == 1
                break;
            else
                a = zeros(sizeVinput(2)); % reset 'a' after decision is made
                b = zeros(sizeVinput(2)); % reset 'b' after decision is made
                Vinput = zeros(sizeVinput); % reset 'Vinput' after decision is made
            end
        end
    end
end
