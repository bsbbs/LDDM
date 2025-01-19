outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/General';
Setup;
dt = .001;
tau_pos = .02;
tau_neg = .02;

% kernel - PP
h = figure; hold on;
filename = 'Kernel_PP';
t_pos = -3*tau_pos/dt:0;
y1 = exp(t_pos*dt / tau_pos);
t_neg = 0:3*tau_neg/dt;
y2 = exp(-t_neg*dt / tau_neg);
plot([t_pos, t_neg], [y1, y2]*100, '-', 'Color',OKeeffe(10,:), 'LineWidth',2);
plot([0,0],[-100, 100], 'k--');
plot([-60,60],[0,0], 'k--');
ylim([-100, 100]);
ylabel('Modification (%)');
xlabel('t_{pre} - t_{post} (ms)', 'FontAngle', 'italic');
text(-50,0.6*100,'Pre-Post');
text(25,0.6*100,'Post-Pre');
mysavefig(h, filename, outdir, 14, [3, 2], 2);

% kernel - PN
h = figure; hold on;
filename = 'Kernel_PN';
t_pos = -3*tau_pos/dt:0;
y1 = exp(t_pos*dt / tau_pos);
t_neg = 0:3*tau_neg/dt;
y2 = -exp(-t_neg*dt / tau_neg);
plot([t_pos, t_neg], [y1, y2]*100, '-', 'Color',OKeeffe(10,:), 'LineWidth',2);
plot([0,0],[-100, 100], 'k--');
plot([-60,60],[0,0], 'k--');
ylim([-100, 100]);
ylabel('Modification (%)');
xlabel('t_{pre} - t_{post} (ms)', 'FontAngle', 'italic');
text(-50,0.6*100,'Pre-Post');
text(25,0.6*100,'Post-Pre');
mysavefig(h, filename, outdir, 14, [3, 2], 2);


% kernel - PP
h = figure; hold on;
filename = 'Kernel_PP2';
t_pos = -3*tau_pos/dt:0;
y1 = exp(t_pos*dt / tau_pos);
t_neg = 0:3*tau_neg/dt;
y2 = exp(-t_neg*dt / tau_neg);
plot([t_pos, t_neg], [y1, y2]*100, '-', 'Color',OKeeffe(10,:), 'LineWidth',2);
plot([0,0],[-100, 100], 'k--');
plot([-60,60],[0,0], 'k--');
ylim([-100, 100]);
ylabel('Modification (%)');
xlabel('t_{post} - t_{pre} (ms)', 'FontAngle', 'italic');
text(-50,0.6*100,'Post-Pre');
text(25,0.6*100,'Pre-Post');
mysavefig(h, filename, outdir, 14, [3, 2], 2);

% kernel - PN
h = figure; hold on;
filename = 'Kernel_PN2';
t_pos = -3*tau_pos/dt:0;
y1 = -exp(t_pos*dt / tau_pos);
t_neg = 0:3*tau_neg/dt;
y2 = exp(-t_neg*dt / tau_neg);
plot([t_pos, t_neg], [y1, y2]*100, '-', 'Color',OKeeffe(10,:), 'LineWidth',2);
plot([0,0],[-100, 100], 'k--');
plot([-60,60],[0,0], 'k--');
ylim([-100, 100]);
ylabel('Modification (%)');
xlabel('t_{post} - t_{pre} (ms)', 'FontAngle', 'italic');
text(-50,0.6*100,'Post-Pre');
text(25,0.6*100,'Pre-Post');
mysavefig(h, filename, outdir, 14, [3, 2], 2);


