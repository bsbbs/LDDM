monkeyD_data = load_data("behavData_dam.mat");
QQ_plot(monkeyD_data, "monkeyD");

monkeyE_data = load_data("behavData_eli.mat");
QQ_plot(monkeyE_data, "monkeyE");

simulate_data = load("PlotData.mat");
empirical_data = load_data("behavData_dam.mat");