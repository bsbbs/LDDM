%% transforming the data from the paper Rafiei and Rahnev - 2021 - Qualitative speed-accuracy tradeoff effects that...
datadir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/Speed-accuracy-tradeoff/Rafiei_Rahnev2021';

Data.Subject = [];
Data.condition = [];
Data.rt = [];
Data.correct = [];
Data.contrast = [];
for subi = 1:30
    dfolder = fullfile(datadir, sprintf('Subject %i', subi));
    Nsessions = dir(fullfile(dfolder, 'Session*'));
    for sessi = 1:numel(Nsessions)
        load(fullfile(dfolder, Nsessions(sessi).name));
        for i = 1:4
            for j = 1:5
                Data.Subject = [Data.Subject; ones(size(p.data{i,j}.rt'))*subi];
                Data.condition = [Data.condition; ones(size(p.data{i,j}.rt'))*p.condition(i,j)];
                Data.rt = [Data.rt; p.data{i,j}.rt'];
                Data.correct = [Data.correct; p.data{i,j}.correct'];
                Data.contrast = [Data.contrast; p.data{i,j}.contrast'];
            end
        end
    end
end
T = struct2table(Data);
writetable(T, fullfile(datadir, 'Transfromed.txt'));