dataBhvr = load_data("behavData_dam.mat");
version = "MonkeyD_accuracy_v1";
condition = "accuracy";

alpha = 34.695025;
beta = 0.978248;
sigma = 16.577909;
S = 350.064641;
tauR = 0.008438;
tauG = 0.112152;
tauD = 0.155874;
params = [alpha, beta, sigma, S, tauR, tauG, tauD];

QP_plot(params, dataBhvr, condition, version);

