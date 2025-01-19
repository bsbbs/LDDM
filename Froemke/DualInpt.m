% Single input situation
%% Setpath code
gnrloutdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/General';
outdir = '/Users/bs3667/Library/CloudStorage/GoogleDrive-bs3667@nyu.edu/My Drive/LDDM/Froemke/DynamicVersionDual';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
outdirSync = fullfile(outdir, 'PPSync');
if ~exist(outdirSync,'dir')
    mkdir(outdirSync);
end
outdirAsync = fullfile(outdir, 'PPAsync');
if ~exist(outdirAsync,'dir')
    mkdir(outdirAsync);
end
Setup;
show = 1;
%% Generating a random network
Networkgenerator;

%% Input tuning
Ninput = 2;
Inputdir = outdir;
InputTuning;

%% commen parameters
iSTDPsign = 1;
dt = .001; % time precision for simulation, in unit of second
Ntrial = 258;
value = 20*ones(1, Ntwk.Input.N);
rng(2023);
Seqraw = CreateEvents(Ntrial, dt);
time = [0:(size(Seq, 1)-1)]*dt/60; % unit in Mins
% time = [0:1000*60]*dt/60;
smpl_time = 5*60+[10000:30000]*dt; % unit in Secs
inpt_time = [0:40000]*dt;
%% Input sequence - Sync 
testdir = outdirSync;
Seq = Seqraw;
Seq(:,2) = Seq(:,1);
Inputdir = testdir;
InputSeq;

% Simulations - Sync
SavedResult = fullfile(testdir, 'SavedResult.mat');
Simulator;
Ntwk.wEE_initial = wEE_initial;
Ntwk.wEI_initial = wEI_initial;
Ntwk.wIE_initial = wIE_initial;
Visualization;

%% Input sequence - Async 
testdir = outdirAsync;
Seq = Seqraw;
Inputdir = testdir;
InputSeq;

%% Simulations - Async
SavedResult = fullfile(testdir, 'SavedResult.mat');
Simulator;
Ntwk.wEE_initial = wEE_initial;
Ntwk.wEI_initial = wEI_initial;
Ntwk.wIE_initial = wIE_initial;
Visualization;


