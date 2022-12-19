function [sortNum, myCluster] = RndCtrl(numNode)
myCluster  = parcluster();
t = datenum(clock)*10^10 - floor(datenum(clock)*100)*10^8;
Myhost = [myCluster.Host '_' num2str(t)];
JobStorageLocation = sprintf('/gpfs/data/glimcherlab/BoShen/RecurrentModel/coordinate/%s/%s',date, Myhost);
mkdir(JobStorageLocation);
myCluster.JobStorageLocation = JobStorageLocation;
pause(.5);
filelist = dir(sprintf('/gpfs/data/glimcherlab/BoShen/RecurrentModel/coordinate/%s',date));
filenames = {filelist.name};
filenames = filenames(3:end);
while length(filenames) < numNode
    pause(.10);
end
sortNum = find(strcmp(filenames, Myhost));
fprintf('sortNum%i \n', sortNum);
fprintf('Host%s \n', Myhost);
