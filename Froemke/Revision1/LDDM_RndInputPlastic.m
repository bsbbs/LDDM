function [R, G, D, Vcourse, winputp, wrgp, wgrp] = LDDM_RndInputPlastic(Vprior, Vinput, BR, BG, winput0, wrg0, wgr0, a, b,...
    sgmR, sgmG, sgmInput, Tau, predur, dur, dt, presentt, triggert, thresh, initialvals, stimdur, stoprule)
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
tauN = 0.002; % time constant for Ornstein-Uhlenbeck process of noise
eta = .00001; % learning rate of synaptic weight on Oja's rule
targetrate = 32; % target firing rate of inhibitory synaptic plasticity
%% define parameters
winput = winput0;
wrg = wrg0;
wgr = wgr0;

pretask_steps = round(predur/dt);
onset_of_stimuli = round(presentt/dt); % align to the beginning of task as t = 0.
stim_duration = round(stimdur/dt);
offset_of_stimuli = onset_of_stimuli + stim_duration;
onset_of_trigger = round(triggert/dt);
posttask_steps = round(dur/dt);
sizeVinput = size(Vinput);
if sizeVinput(1) > 1 error('Error: the size of Vinput has to be 1xN'); end
rt = NaN;
choice = NaN;
smplrate = 1000; % per seccond
smplsz = round((predur+dur)/smplrate);
winputp = nan(smplsz,sizeVinput(2));
wrgp = nan(smplsz,sizeVinput(2), sizeVinput(2));
wgrp = nan(smplsz,sizeVinput(2));
%% stablizing noise for 200 ms
InoiseG = zeros(sizeVinput);
InoiseR = zeros(sizeVinput);
stablizetime = round(.2/dt);
for kk = 1:stablizetime
    % update noise
    InoiseG = InoiseG + (-InoiseG + randn(sizeVinput).*sqrt(dt).*sgmG)/tauN*dt;
    InoiseR = InoiseR + (-InoiseR + randn(sizeVinput).*sqrt(dt).*sgmR)/tauN*dt;
end
G = initialvals(2,:) + InoiseG;
R = initialvals(1,:) + InoiseR;
D = initialvals(3,:);
%% simulation begin
Vnoise = randn(sizeVinput)*sgmInput;
t_stamp = pretask_steps + 1;
for ti = (-pretask_steps):posttask_steps % align the beginning of the task as ti = 0
    % input values
    if ti > -pretask_steps && ti < 0
        V = Vprior;
    elseif ti >= onset_of_stimuli && ti < offset_of_stimuli
        if (mod(ti*dt, .005) == 0)
            Vnoise = randn(sizeVinput)*sgmInput;
        end
        V = Vinput + Vnoise;
    else
        V = zeros(sizeVinput);
    end
    Vcourse(ti+t_stamp,:) = V;
    
    % update R, G, I
    dR = (-R(ti+t_stamp,:)' + ((winput.*V)' + a*R(ti+t_stamp,:)' + BR)./(1+wgr.*G(ti+t_stamp,:)'))/Tau(1)*dt;
    dG = (-G(ti+t_stamp,:)' + wrg*R(ti+t_stamp,:)'  + BG - D(ti+t_stamp,:)')/Tau(2)*dt;
    dD = (-D(ti+t_stamp,:)' + b*(ti >= onset_of_trigger)*R(ti+t_stamp,:)')/Tau(3)*dt;

    dwinput = eta*(R(ti+t_stamp,:).*V - winput.*(R(ti+t_stamp,:).^2))*dt;
    dwrg = eta*(G(ti+t_stamp,:)'*R(ti+t_stamp,:) - diag(G(ti+t_stamp,:).^2)*wrg)*dt;
    dwgr = eta*(G(ti+t_stamp,:).*(R(ti+t_stamp,:) - targetrate))*dt;

    R(ti+t_stamp+1,:) = R(ti+t_stamp,:) + dR' + InoiseR;
    G(ti+t_stamp+1,:) = G(ti+t_stamp,:) + dG' + InoiseG;
    D(ti+t_stamp+1,:) = D(ti+t_stamp,:) + dD';

    winput = winput + dwinput;
    wrg = wrg + dwrg;
    wgr = wgr + dwgr;
    
    % update noise
    InoiseG = InoiseG + (-InoiseG + randn(sizeVinput).*sqrt(dt).*sgmG)/tauN*dt;
    InoiseR = InoiseR + (-InoiseR + randn(sizeVinput).*sqrt(dt).*sgmR)/tauN*dt;
    % setting lower boundary, forcing neural firing rates to be non-negative
    G(ti+t_stamp+1,G(ti+t_stamp+1,:) < 0) = 0;
    D(ti+t_stamp+1,D(ti+t_stamp+1,:) < 0) = 0;
    R(ti+t_stamp+1,R(ti+t_stamp+1,:) < 0) = 0;
    winput(winput < 0) = 0;
    wrg(wrg < 0) = 0;
    wgr(wgr < 0) = 0;
    
    if mod((ti+t_stamp)*dt, 1/smplrate) == 0
        si = round((ti+t_stamp)*dt*smplrate);
        winputp(si,:) = winput;
        wrgp(si,:,:) = wrg;
        wgrp(si,:) = wgr;
    end
    % threshold detecting
    if ti > onset_of_trigger && isnan(rt)
        if max(R(ti+t_stamp,:)) >= thresh
            rt = (ti-onset_of_trigger)*dt;
            choice = find(R(ti+t_stamp,:) == max(R(ti+t_stamp,:)));
            if stoprule == 1
                break;
            else
                a = zeros(sizeVinput(2)); % reset 'a' after decision was made
                b = zeros(sizeVinput(2)); % reset 'b' after decision was made
                Vinput = zeros(sizeVinput); % reset 'Vinput' after decision was made
            end
        end
    end
end
