% Single input situation
%% Setpath code
gnrloutdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/General';
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionSngl2';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
outdirPP = fullfile(outdir, 'PP');
if ~exist(outdirPP,'dir')
    mkdir(outdirPP);
end
outdirPN = fullfile(outdir, 'PN');
if ~exist(outdirPN,'dir')
    mkdir(outdirPN);
end
Setup;
show = 1;
%% Generating a random network
Networkgenerator;

%% Input tuning
Ninput = 1;
Inputdir = outdir;
InputTuning;

%% Input value(s) and sequence(s)
dt = .001; % time precision for simulation, in unit of second
Ntrial = 258;
rng(2023);
Seq = CreateEvents(Ntrial, dt);
Seq = Seq(:, 1:Ntwk.Input.N);
value = 20*ones(1, Ntwk.Input.N);
time = [0:(size(Seq, 1)-1)]*dt/60; % unit in Mins
% time = [0:1000*60]*dt/60;
smpl_time = 5*60+[10000:30000]*dt; % unit in Secs
inpt_time = [0:40000]*dt;
InputSeq;

%% Simulations - PP
testdir = outdirPP;
iSTDPsign = 1;
SavedResult = fullfile(testdir, 'SavedResult.mat');
Simulator;
 Ntwk.wEE_initial = wEE_initial;
    Ntwk.wEI_initial = wEI_initial;
    Ntwk.wIE_initial = wIE_initial;
Visualization;
%% Simulations - PN
testdir = outdirPN;
iSTDPsign = -1;
SavedResult = fullfile(outdirPN, 'SavedResult.mat');
Simulator;
Visualization;

