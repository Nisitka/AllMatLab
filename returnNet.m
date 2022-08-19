function out = returnNet

Dir = 'D:\Наука\Programs\Programs MatLab\Build NetWorks\SegmentationNet'; %путь к проекту

data = load(fullfile(Dir, 'Net.mat')); %ипорт сохраненного WorkSpace

net = data.net; 
clear data;
clear Dir;

out = net;