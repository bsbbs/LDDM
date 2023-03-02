% Fig. 6 - Fit LDDM, RNM, and LCA to Roitman & Shadlen, 2002
%%%%%%%%%%%%%
% In order to implement the model fitting code, you need to 
% 1. adapt the directory settings on the top of each code to your local
% environment;
% 2. download the optimization algorithm "bads" from https://github.com/acerbilab/bads
% 3. make sure the matlab you are using has equiped with gpu computation,
% to test this, type "gpuDevice" in the command window of matlab. If gpu is
% in function, it will list the information of the gpu card.

cd ../Fit;
%% Fit LDDM
main_FitBhvr7ParamsIV_QMLE_svr;
LLSpace_LDDM7IVG0;
PrmtrsRcvry_LDDM7IV;

%% Fit RNM
main_WW06FitBhvr8ParamsVI_QMLE_svr;

%% Fit LCA
main_FitLCA;