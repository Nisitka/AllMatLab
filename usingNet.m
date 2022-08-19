%__________________________программа для тестирования сети_________________

Dir = 'D:\Наука\Programs\Programs MatLab\Build NetWorks\SegmentationNet'; %путь к проекту

data = load(fullfile(Dir, 'Net.mat')); %ипорт сохраненного WorkSpace

net = data.net; 
clear data;

newClassNames = ["road","sky","vehicle","pedestrian","background"];

TestImagesDir = fullfile(Dir, 'TestImages');  
TestImages = imageDatastore(TestImagesDir);

number = 5;
I = imread(TestImages.Files{number});
figure
imshow(I);
I = imresize(I, [720, 960]);
figure
imshow(I);

C = semanticseg(I,net);    %сеть возвращает категориальный массив 
cmap = jet(numel(newClassNames));
B = labeloverlay(I,C,'Colormap',cmap);
figure
imshow(B)

N = numel(newClassNames);
ticks = 1/(N*2):1/N:1;
colorbar('TickLabels',cellstr(newClassNames),'Ticks',ticks,'TickLength',0,'TickLabelInterpreter','none');
colormap(cmap)
