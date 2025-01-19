% Excitatory neurons, tuning to the inputs, binary as an example
%% setup
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionLIFPPSync';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
addpath('./utils/');
addpath('./functions/');
addpath('../utils/bluewhitered/');
Setup;
%% Generate the network from single neurons
Networkgenerator;
%% define kernels
% Kernel - Integrate and fire model + STDP sliding
Ntwk.VL = -70; % mV, resting potential
Ntwk.Vreset = -55; % mV, threshold/reset potential (Wang 2002; Burkitt et al., 2004)
Ntwk.VE = 0; % mV, excitatory synaptic potential
Ntwk.VI = -80; % mV, inhibitory synaptic potential (Vogels et al. 2011; Burkitt et al., 2004)
Ntwk.Exct.Cm = .5; % nF, membrane capacity for pyramidal neurons (Wang 2002)
Ntwk.Inhbt.Cm = .2; % nF, membrane capacity for interneurons (Wang 2002)
Ntwk.Exct.gL = 25; % nS, membrane leaky conductance for pyramidal neurons (Wang 2002)
Ntwk.Inhbt.gL = 20; % nS, membrane leaky conductance for interneurons (Wang 2002)
Ntwk.Exct.taum = Ntwk.Exct.Cm/Ntwk.Exct.gL; % .02s, membrane time constant of pyramidal neurons (Wang, 2002; Vogels et al., 2011; Burkitt et al., 2004)
Ntwk.Inhbt.taum = Ntwk.Inhbt.Cm/Ntwk.Inhbt.gL; % .01s, membrane time constant of interneurons (Wang 2002)
Ntwk.tauAMPA = .002; % s, the decay time constant of AMPA currents. (Wang 2002)
Ntwk.tauNMDA.rise = .002; % s, the rising time constant of NMDA currents. (Wang 2002)
Ntwk.tauNMDA.decay = .1; % s, the decay time constant of NMDA currents. (Wang 2002)
Ntwk.tauGABA = .005; % s, the decay time constent of GABA currents. (Wang 2002)
Ntwk.tauREF = .005; % s, refractory period (Vogels et al., 2011)
Ntwk.Exct.tauREF = .002; % s, refractory period for pyramidal neurons (Wang 2002)
Ntwk.Inhbt.tauREF = .001; % s, refractory period for interneurons (Wang 2002)

V0 = Ntwk.VL;
V = V0;
for t = 1:timesteps
    dV = (-gL*(V(t) - VL) - Isyn(t))*dt/Cm;
    Isyn(t) = IrecNMDA(t) + IrecGABA(t);
    IrecNMDA(t) = gNMDA(V(t) - VE)*sum(w*sNMDA(:,t));
    IrecGABA(t) = gGABA*(V(t) - VI)*sum(sGABA(:,t));
    dsNMDA(t) = -sNMDA(t)/tauNMDA.decay + alpha* x;
end

