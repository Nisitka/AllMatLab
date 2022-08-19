%Загрузка данных тренировки.
dataSetDir = fullfile('C:\Program Files\Polyspace\R2020b\toolbox\vision\visiondata\Street segmentation'); 
imageDir = fullfile(dataSetDir,'trainingImages', '701_StillsRaw_full');
labelDir = fullfile(dataSetDir,'trainingLabels');

%Создайте pixelLabelDatastore для заземляющих пиксельных меток истины.

classNames = [...
    "Animal", ...
    "Archway", ...
    "Bicyclist", ...
    "Bridge", ...
    "Building", ...
    "Car", ...
    "CartLuggagePram", ...
    "Child", ...
    "Column_Pole", ...
    "Fence", ...
    "LaneMkgsDriv", ...
    "LaneMkgsNonDriv", ...
    "Misc_Text", ...
    "MotorcycleScooter", ...
    "OtherMoving", ...
    "ParkingBlock", ...
    "Pedestrian", ...
    "Road", ...
    "RoadShoulder", ...
    "Sidewalk", ...
    "SignSymbol", ...
    "Sky", ...
    "SUVPickupTruck", ...
    "TrafficCone", ...
    "TrafficLight", ...
    "Train", ...
    "Tree", ...
    "Truck_Bus", ...
    "Tunnel", ...
    "VegetationMisc", ...
    "Wall"];
labelIDs = [...
    064 128 064; ... % "Animal"
    192 000 128; ... % "Archway"
    000 128 192; ... % "Bicyclist"
    000 128 064; ... % "Bridge"
    128 000 000; ... % "Building"
    064 000 128; ... % "Car"
    064 000 192; ... % "CartLuggagePram"
    192 128 064; ... % "Child"
    192 192 128; ... % "Column_Pole"
    064 064 128; ... % "Fence"
    128 000 192; ... % "LaneMkgsDriv"
    192 000 064; ... % "LaneMkgsNonDriv"
    128 128 064; ... % "Misc_Text"
    192 000 192; ... % "MotorcycleScooter"
    128 064 064; ... % "OtherMoving"
    064 192 128; ... % "ParkingBlock"
    064 064 000; ... % "Pedestrian"
    128 064 128; ... % "Road"
    128 128 192; ... % "RoadShoulder"
    000 000 192; ... % "Sidewalk"
    192 128 128; ... % "SignSymbol"
    128 128 128; ... % "Sky"
    064 128 192; ... % "SUVPickupTruck"
    000 000 064; ... % "TrafficCone"
    000 064 064; ... % "TrafficLight"
    192 064 128; ... % "Train"
    128 128 000; ... % "Tree"
    192 128 192; ... % "Truck_Bus"
    064 000 064; ... % "Tunnel"
    192 192 000; ... % "VegetationMisc"
    064 192 000];    % "Wall"

%Визуализация пиксельных меток для одного из изображений CamVid.
labels = imread(fullfile(labelDir,'0001TP_006690_L.png'));
figure
imshow(labels)

% Add colorbar to show class to color mapping.
N = numel(classNames);
ticks = 1/(N*2):1/N:1;
colorbar('TickLabels',cellstr(classNames),'Ticks',ticks,'TickLength',0,'TickLabelInterpreter','none');
colormap(labelIDs./255)

%Создание imageDatastore, чтобы загрузить изображения CamVid:
imds = imageDatastore(fullfile(dataSetDir,'trainingImages','701_StillsRaw_full'));
%Создание pixelLabelDatastore, чтобы загрузить пиксельные метки CamVid:
pxds = pixelLabelDatastore(labelDir,classNames,labelIDs);

I = readimage(imds,10);
C = readimage(pxds,10);

B = labeloverlay(I,C,'Colormap',labelIDs./255);
figure
imshow(B)

% Add a colorbar.
N = numel(classNames);
ticks = 1/(N*2):1/N:1;
colorbar('TickLabels',cellstr(classNames),'Ticks',ticks,'TickLength',0,'TickLabelInterpreter','none');
colormap(labelIDs./255)

%---Неопределенные или пустые метки
%Пикселю маркированные наборы данных свойственно включать "неопределенные" или "пустые" метки. 
%Они используются, чтобы определять пиксели, которые не были маркированы. Например, в CamVid, метка ID [0 0 0] используется, чтобы определять "пустой" класс. 
%Учебные алгоритмы и алгоритмы оценки, как ожидают, не будут включать эти метки ни в какие вычисления.
%"Пустой" класс нельзя явным образом назвать при использовании pixelLabelDatastore. 
%Любая метка ID, которая не сопоставлена с именем класса, автоматически маркирована "неопределенной" и исключена из вычислений. 
%Чтобы видеть неопределенные пиксели, используйте isundefined, чтобы создать маску и затем отобразить его сверху изображения.

undefinedPixels = isundefined(C);
B = labeloverlay(I,undefinedPixels);
figure
imshow(B)
title('Undefined Pixel Labels')

%---Объединение классов
newClassNames = ["road","sky","vehicle","pedestrian","background"];

groupedLabelIDs = {
    % road
    [
    128 064 128; ... % "Road"
    128 000 192; ... % "LaneMkgsDriv"
    192 000 064; ... % "LaneMkgsNonDriv"
    000 000 192; ... % "Sidewalk" 
    064 192 128; ... % "ParkingBlock"
    128 128 192; ... % "RoadShoulder"
    ]
   
    % "sky"
    [
    128 128 128; ... % "Sky"
    ]
    
    % "vehicle"
    [
    064 000 128; ... % "Car"
    064 128 192; ... % "SUVPickupTruck"
    192 128 192; ... % "Truck_Bus"
    192 064 128; ... % "Train"
    000 128 192; ... % "Bicyclist"
    192 000 192; ... % "MotorcycleScooter"
    128 064 064; ... % "OtherMoving"
    ]
     
    % "pedestrian"
    [
    064 064 000; ... % "Pedestrian"
    192 128 064; ... % "Child"
    064 000 192; ... % "CartLuggagePram"
    064 128 064; ... % "Animal"
    ]
    
    % "background"      
    [
    128 128 000; ... % "Tree"
    192 192 000; ... % "VegetationMisc"    
    192 128 128; ... % "SignSymbol"
    128 128 064; ... % "Misc_Text"
    000 064 064; ... % "TrafficLight"  
    064 064 128; ... % "Fence"
    192 192 128; ... % "Column_Pole"
    000 000 064; ... % "TrafficCone"
    000 128 064; ... % "Bridge"
    128 000 000; ... % "Building"
    064 192 000; ... % "Wall"
    064 000 064; ... % "Tunnel"
    192 000 128; ... % "Archway"
    ]
    };

%Создание pixelLabelDatastore с помощью нового класса и метки IDs:
pxds = pixelLabelDatastore(labelDir,newClassNames,groupedLabelIDs);

%Считывание 10-ти пиксельное изображение метки и отображение его сверху изображения.
C = readimage(pxds,10);
cmap = jet(numel(newClassNames));
B = labeloverlay(I,C,'Colormap',cmap);
figure
imshow(B)

% add colorbar
N = numel(newClassNames);
ticks = 1/(N*2):1/N:1;
colorbar('TickLabels',cellstr(newClassNames),'Ticks',ticks,'TickLength',0,'TickLabelInterpreter','none');
colormap(cmap)

%------------------------------Создание нейронной сети---------------------

%Инцелизация входного слоя изображения:
%%%Семантическая сеть сегментации запускается с imageInputLayer, 
%который задает самый маленький размер изображения, который может обработать сеть. 
%Большинство семантических сетей сегментации является полностью сверточным, что означает, 
%что они могут обработать изображения, которые больше, чем заданный входной размер. 
%Здесь, размер изображения [720 960 3] используется для сети к процессу 64x64 изображения RGB.

inputSize = [720 960 3];

numFilters = 35;
filterSize = 3;
numClasses = numel(newClassNames);

layers = segnetLayers(inputSize, numClasses, 3);

opts = trainingOptions('sgdm', ...
    'InitialLearnRate',1e-3, ...
    'MaxEpochs',25, ...
    'MiniBatchSize',2, ...
    'Verbose',false, ...
    'ExecutionEnvironment','cpu', ...
	'Plots','training-progress');

%Создайте пиксельный datastore метки изображений, который содержит данные тренировки.
trainingData = pixelLabelImageDatastore(imds,pxds);

net = trainNetwork(trainingData,layers,opts);


