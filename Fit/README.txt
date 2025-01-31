Version description:
--
FitBhrv, LDDM series
FitBhvr: Fitting reaction time and choice only;
FitDynamic: Fitting the temporal dynamic of neural firing rates
FitDynamicBhvr: Fitting RT, choice and the dynamic of neural firing rates at the same time
NParams: number of free parameters

I,II,III, IV version notes: I & II is similar, assume pre-stimuli period as 0s, with gap of input 90ms (controlled by parameter ndt). During the gap, neural firing rates drop from initial value, which is set as 35(Hz) in versionII and calculated based on Roitman's data at versionI. A decision is made when one side of neural firing rates reach 70(Hz) in version II and a specific number calculated based on Roitman's data in versionI.
Version III assumes the accumulation process starts after the initial dip and recovery, which is 90ms(start of dip) + ~100ms(end of recover) after the onset of stimuli according to Roitman's paper. The initial value at this start point is set as 42(Hz) based on Roitman's data. Threshold is set as 70Hz. For RT, another 30ms is added to the neural dynamic, which captures the motor output delay reported in Roitman's paper. 
Version IV assumes the accumulation process starts from the bottom of initial dip, which is 90ms later than the onset of stimuli according to RoitmanShadlen2002. The initial values were set as 32 Hz.

Number of parameters: 4-parameter versions kept free parameters as alpha, beta, noise, and input scaling. tauR, tauG, and tauI were all set as 100ms. 7-parameter versions kept free parameters as alpha, beta, noise, input scaling, tauR, tauG, and tauI.

Cut1,2, and 3:In order to clean the neural dynamic in Dynamic fitting, there are many versions of try.
Most of the file, the target function is defined to minimize the curve of histogram between empirical and predicted RT distribution. However, this method is thought having problem because it is not considering the weighting of observation number within each histogram bin. As an outcome, the tail and head of the distribution, which has less numbers of trials than the body of distribution is overweighted, causing a problem of suboptimal fitting on meanRT and choice.
Algorithm QMLE (Quantile Maximum likelihood estimation) is based on reference: Heathcote & Australia, and Mewhort, 2002.
4cpParams means 4 free parameters (alpha, beta, noise sigma, and scale) + 6 coherence levels.

FitBhvr7ParamsV: Free parameters: a, b, sgm, scale, tauR, tauG, tauI; thresh = 70, initial values = [42, 42; 2*42, 2*42; 0, 0]. Non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
FitBhvr6ParamsVI: Free parameters: a, b, sgm,       tauR, tauG, tauI; thresh = 70, initial values = [42, 42; (2w-b)*42, (2w-b)*42; b*42, b*42]. scale = (2w-b)*eqlb.^2 + (1-a).*eqlb; Non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
FitBhvr7ParamsVII: Free parameters: a, b, sgm, B0,  tauR, tauG, tauI; thresh = 70, initial values = [42, 42; (2w-b)*42, (2w-b)*42; b*42, b*42]. scale = (2w-b)*eqlb.^2 + (1-a).*eqlb; Non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
FitBhvr7ParamsVIII: Free parameters: a, b, sgm, B0,  tauR, tauG, tauI; thresh = 70, initial values = [32, 32; (2w-b)*32, (2w-b)*32; b*32, b*32]. scale = (2w-b)*eqlb.^2 + (1-a).*eqlb; Non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
FitBhvr7ParamsIX: Free parameters: a, b, sgm, Star, tauR, tauG, tauI; thresh = 70, fitted initial values, scale = max([5, (((2*mean(w,'all') - params(2)))*eqlb.^2 + (1-a(1)).*eqlb)]); Non-decision time = 90 + 30 ms; sims = 1024*5; iters = 20*8. 
FitBhvr7ParamsX: Free parameters: a, b, sgmInput, scale, tauR, tauG, tauI; thresh = 70, sgm = .01; initial values = [32, 32; 2*32, 2*32; 0, 0]; Non-decision time = 90 + 30 ms; sims = 10240; Updated BADS package to 2022 version.
FitDynmc7Params: fitting the neural dynamics only, using OLS. Free parameters: a, b, sgm, tauR, tauG, tauD, threshold. Rstar = 43.1026, initial values = Rstar*[1,1; 2-b, 2-b; b, b]. Scale = (2w-b)*Rstar.^2 + (1-a).* Rstar, non-decision time = 190 + 30 ms. Fitting the dynamics for two parts: 1. Dynamics sorting to the onset of stimuli starting 190 ms after the onset, until median RT, excluding any trace within 100ms of saccade; 2. Dynamics sorting to saccade, leaving 30ms right before saccade and trace back until median RT, exclude any trace within 200ms from stimuli onset.
FitBhvr7Params: equivalent to FitDynmc7Params, but fit behavior only.

--
AsymW series is for the model with asymmetric w weights (Li's model)
AsymWFitBhvr4pars: assumes the accumulation process start 90ms after the stimuli onset. The initial values on R were set for 32 Hz based on Roitman&Shadlen's data. Estimated parameters are alpha, w_ii, noise amplitute, and scaling of input. w_ij was set as 1. non-decision time was considered at the side of action for 30ms. With the gap in the beginning, there was in total 100ms non-decision time.
AsymWFitBhvr6pars: similar to the above. Set the time constant for R and G as free parameters.

--
WW06, Wong & Wang 2006 series of fitting
WW06_5Params: free parameters: Self excitation and lateral inhibition in the matrix of JN, baseline input I0, noise parameter sigma, and input scaling miu0. Other parameters kept the same as in the paper of wong&wang 2006 (referred as the paper hereafter). Initial values were set as H = [2, 2], the same as in the paper. Threshold was set as 15Hz, the same as in the paper. Non-decision time before initial dip and after hitting the threshold were considered as 90 and 30 ms, respectively.
WW06_7Params: Similar to the 5-parameter version mentioned above. Set the tauNMDA and tauAMPA as free parameters. Originally in the paper, tauNMDA = 100ms, and tauAMPA = 2 ms.
WW06_5ParamsII: similar to WW06_5Params but made the following changes: The initial values were set as [32, 32] Hz for the two excitatory pools. The threshold was set as 70Hz. These changes were aimed to match the empirical evidence of Roitman & Shadlen, 2002. The parameter constraints for miu0 was relaxed from 0-60 to 0-120. Number of iteration in bads fitting was reduced from 40*4 = 160 to 20*4 = 80 because of computational power draining.
WW06_5ParamsII: similar to WW06_7Params but made the following changes: The initial values were set as [32, 32] Hz for the two excitatory pools. The threshold was set as 70Hz. These changes were aimed to match the empirical evidence of Roitman & Shadlen, 2002. The parameter constraints for miu0 was relaxed from 0-60 to 0-120. Number of iteration in bads fitting was reduced from 20*8 = 160 to 20*4 = 80 because of computational power draining.
WW06_6ParamsIII: Free parameters: JNn, JNp, sgm, I0, miu0, tauNMDA; fixed parameters: tauAMPA = .002, gamma = .64. Initial value H0 = 2, S0 based on H0. Threshold = 15Hz, non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
WW06_7ParamsIV: Free parameters: JNn, JNp, sgm, I0, miu0, gamma, tauNMDA; fixed parameters: tauAMPA = .002. Initial value H0 = 42/70*15, S0 based on H0. Threshold = 15Hz, non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
WW06_7ParamsV: Free parameters: JNn, JNp, sgm, I0, miu0, gamma, tauNMDA; fixed parameters: tauAMPA = .002. Initial value H0 = 32/70*15, S0 based on H0. Threshold = 15Hz, non-decision time = 90 + 30 ms. Strings = 10240, iterations = 20*16.
