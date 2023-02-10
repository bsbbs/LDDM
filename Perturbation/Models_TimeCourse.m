% Causal manipulation of GABAnergic based signal transfermation
addpath('../../../CoreFunctions');
outdir = './graphics';
if ~exist(outdir,'dir')
    mkdir(outdir);
end
lwd = 2.0;
mksz = 18;
fontsize = 14;
aspect1 = [3.9,2.2]; % 16:9, for wide temporal dynamic
aspect2 = [3 3]; % for temporal dynamic
aspect3 = [2.8 2.54];
GABA = [1, 3.8];
%%
for sgmi = 1
    %% the local disinhibition model
    name = 'LcD';
    timeCourse;
    
    %ValueRepresentation;
    %% Asymetric weigthing
    name = 'AsymW';
    timeCourse;
    
    %% XJ's two-variable model
    name = 'XJ';
    timeCourse;
    %ValueRepresentation;
end