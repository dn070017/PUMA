% 2017-04-08
% marieke kuijjer
% load lioness networks and print to txt file

outtag='output'
cd(outtag)
mats=dir;
mats=mats(3:length(mats)) % the first dir are "." and ".."
for(i=1:length(mats))
  load(['lioness.' num2str(i) '.mat'])
  dlmwrite(['../txt/lioness.', num2str(i), '.txt'], PredNet, '\t');  
end
