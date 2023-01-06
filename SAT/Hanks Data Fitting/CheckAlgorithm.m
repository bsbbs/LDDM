% Check on the optimization algorithm
dat_dir = '/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks';
addpath('~/Documents/LDDM/utils');

h = figure;
% MonkeyE - speed
% Old Algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyE_speed_Svr/CollectRslts104901676.mat');
VslzHanks(Collect, h, 1);
% New algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyE_speed_Svr/MnkyE_Speed_CollectRslts13301600.mat');
VslzHanks(Collect,  h, 2);

% MonkeyE - accuracy
% Old Algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyE_accuracy_Svr/CollectRslts85528855.mat');
VslzHanks(Collect,  h, 3);
% New algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyE_accuracy_Svr/MnkyE_Speed_CollectRslts103311583.mat');
VslzHanks(Collect,  h, 4);


% MonkeyD - speed

% Old Algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyD_speed_Svr/CollectRslts77407851.mat');
VslzHanks(Collect,  h, 5);
% New algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyD_speed_Svr/MnkyD_Speed_CollectRslts45089019.mat');
VslzHanks(Collect,  h, 6);


% MonkeyD - accuracy
% Old Algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyD_accuracy_Svr/CollectRslts98948255.mat');
VslzHanks(Collect,  h, 7);
% New algorithm
load('/Volumes/GoogleDrive/My Drive/LDDM/SAT/Hanks/monkeyD_accuracy_Svr/MnkyD_Accuracy_CollectRslts84289691.mat');
VslzHanks(Collect,  h, 8);

%% Visualization
function h = VslzHanks(Collect, h, hi)
figure(h);
subplot(4,2,hi); hold on;
Pls = nan(8,160);
for i = 1:160
    Pls(1:7,i) = Collect(i).xest;
    Pls(8,i) = Collect(i).fval;
end
% h = figure; hold on;
scatter3(Pls(1,:), Pls(2,:), Pls(8,:), 190, Pls(8,:), 'Marker','.');
loc = find(Pls(8,:) == min(Pls(8,:)));
scatter3(Pls(1,loc), Pls(2,loc), Pls(8,loc), 190, 'r', 'Marker','*');
xlabel('alpha');
ylabel('beta');
zlabel('nLL');
grid on;
view([20, 10]);
end


