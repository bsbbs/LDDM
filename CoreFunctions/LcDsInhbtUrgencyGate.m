function [R, G, I, rt, choice] = LcDsInhbtUrgencyGate(Vinput, w, a, b1, b2, sgm, Tau, dur, dt, presentt, triggert1, triggert2, thresh, initialvals, stimdur, stoprule)
%%%%%%%%%%%%%%%%%%
% The core function of local disinhibitory model
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
onset_of_trigger1 = round(triggert1/dt);
onset_of_trigger2 = round(triggert2/dt);
stim_duration = round(stimdur/dt);
sizeVinput = size(Vinput);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end
rt = NaN;
choice = NaN;
%% simulation begin
G(1,:) = initialvals(2,:);
R(1,:) = initialvals(1,:);
I(1,:) = initialvals(3,:);
InoiseG = zeros(sizeVinput);
InoiseR = zeros(sizeVinput);
InoiseI = zeros(sizeVinput);
for ti = 2:total_time_steps
    % update R, G, I
    dG = (-G(ti-1,:)' + w * R(ti-1,:)' - I(ti-1,:)')/Tau(2)*dt;
    dI = (-I(ti-1,:)' + ((ti >= onset_of_trigger1)*b1*min((ti - onset_of_trigger1)/(onset_of_trigger2 - onset_of_trigger1),1)+ (b2-b1)*(ti >= onset_of_trigger2))*R(ti-1,:)')/Tau(3)*dt;
    dR = (-R(ti-1,:)' + (Vinput'*(ti >= onset_of_stimuli & ti < onset_of_stimuli+stim_duration) + a*(ti >= onset_of_stimuli)*R(ti-1,:)')./(1+G(ti-1,:)'))/Tau(1)*dt;
    G(ti,:) = G(ti-1,:) + dG' + InoiseG(ti-1,:);
    I(ti,:) = I(ti-1,:) + dI' + InoiseI(ti-1,:);
    R(ti,:) = R(ti-1,:) + dR' + InoiseR(ti-1,:);
    % update noise
    dInoiseG = (-InoiseG(ti-1,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    dInoiseI = (-InoiseI(ti-1,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    dInoiseR = (-InoiseR(ti-1,:) + randn(sizeVinput).*sqrt(dt).*sgm)/tauN*dt;
    InoiseG(ti,:) = InoiseG(ti-1,:) + dInoiseG;
    InoiseI(ti,:) = InoiseI(ti-1,:) + dInoiseI;
    InoiseR(ti,:) = InoiseR(ti-1,:) + dInoiseR;
    % setting lower boundary, forcing neural firing rates to be non-negative
    G(ti,G(ti,:) < 0) = 0;
    I(ti,I(ti,:) < 0) = 0;
    R(ti,R(ti,:) < 0) = 0;
    % threshold detecting
    if ti > onset_of_stimuli && isnan(rt)
        if max(R(ti,:)) >= thresh
            rt = (ti-onset_of_stimuli)*dt;
            choice = find(R(ti,:) == max(R(ti,:)));
            if stoprule == 1
                break;
            else
                a = zeros(sizeVinput(2)); % reset 'a' after decision is made
                b2 = zeros(sizeVinput(2)); % reset 'b' after decision is made
                Vinput = zeros(sizeVinput); % reset 'Vinput' after decision is made
            end
        end
    end
end
