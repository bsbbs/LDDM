function [sortNum, myCluster] = RndCtrl2(numNode)
myCluster  = parcluster();
JobStorageLocation = sprintf('/gpfs/data/glimcherlab/BoShen/RecurrentModel/coordinate2/%s/%s',date,myCluster.Host);
mkdir(JobStorageLocation);
myCluster.JobStorageLocation = JobStorageLocation;
pause(.5);
filelist = dir(sprintf('/gpfs/data/glimcherlab/BoShen/RecurrentModel/coordinate2/%s',date));
filenames = {filelist.name};
filenames = filenames(3:end);
while length(filenames) < numNode
    pause(.10);
end
sortNum = find(strcmp(filenames,myCluster.Host));
fprintf('sortNum%i \n', sortNum);
fprintf('Host%s \n', myCluster.Host);
