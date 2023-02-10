%% Test potentiation on inhibitory connection weights (simulating GABAergic agonist)
% This coded compared three different candidates, i.e., LDDM (or LcD in the
% code), RNM (or Wong06 in the code) and AsymW. Only the LDDM and RNM are
% reported and discussed in Shen et al., 2023
addpath('../Perturbation');
%% visualize the dynamics under control and inhibitory potentiation condition
Models_TimeCourse;

%% RT and choice under different levels of inhibitory weights
Models_GABA_Activation_ACC_RT;

%% Examine a wider parameter regime of the LDDM
% This step requires GPU computation
LcD_GABA_Widerspace_2mtrx_Dsktp;