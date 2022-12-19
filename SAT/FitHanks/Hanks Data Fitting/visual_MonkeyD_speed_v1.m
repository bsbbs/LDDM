dataBhvr = load_data("behavData_dam.mat");
version = "MonkeyD_speed_v1";
condition = "speed";

alpha = 33.842539;
beta = 1.821676;
sigma = 35.723187;
S = 2368.22093;
tauR = 0.33628;
tauG = 0.001002;
tauD = 0.226118;
params = [alpha, beta, sigma, S, tauR, tauG, tauD];

QP_plot(params, dataBhvr, condition, version);