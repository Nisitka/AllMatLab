%_________________retraining net___________________________________________

Dir = 'D:\Наука\Programs MatLab\Build NetWorks\SegmentationNet'; %путь к проекту

%Загрузка данных тренировки.
dataSetDir = fullfile('C:\Program Files\Polyspace\R2020b\toolbox\vision\visiondata\Street segmentation'); 
imageDir = fullfile(dataSetDir,'trainingImages', '701_StillsRaw_full');
labelDir = fullfile(dataSetDir,'trainingLabels');

%Создание imageDatastore, чтобы загрузить изображения CamVid:
imds = imageDatastore(fullfile(dataSetDir,'trainingImages','701_StillsRaw_full'));

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

%------------------------------Создание нейронной сети---------------------

inputSize = [720 960 3];
numClasses = numel(newClassNames);

layers = resnet50('Weights','none');  %импортируем слои без коофициентов сети resnet50

%удаление не нужных слоев
layers = removeLayers(layers, {'input_1'...
                               'conv1'...
                               'avg_pool' ...
                               'ClassificationLayer_fc1000'...
                               'fc1000_softmax'...
                               'fc1000'...
                               'max_pooling2d_1'});                                    
%добавление необходимых слоев
filterSize = [3 3];
numFilters = 64;
 
layers = addLayers(layers, [imageInputLayer(inputSize, 'Name', 'input'), ...
                            convolution2dLayer(filterSize, numFilters, 'Name', 'conv1', 'Padding','same')]);

layers = addLayers(layers, [maxUnpooling2dLayer('Name','unpool'), ...
                            convolution2dLayer(3, numClasses,'Padding','same', 'Name', 'convOutput'), ...
                            softmaxLayer('Name', 'softmax'), ...
                            pixelClassificationLayer('Name', 'pixelClass')]);
                        
layers = addLayers(layers, maxPooling2dLayer(3, 'Name', 'maxPool', 'Stride', 3, 'HasUnpoolingOutputs', true));

%соединения добавленных слоев
layers = connectLayers(layers, 'conv1', 'bn_conv1');
layers = connectLayers(layers, 'activation_1_relu', 'maxPool');
layers = connectLayers(layers, 'maxPool/out', 'res2a_branch1');
layers = connectLayers(layers, 'maxPool/out', 'res2a_branch2a');

layers = connectLayers(layers, 'activation_49_relu', 'unpool/in'); 

layers = connectLayers(layers,'maxPool/indices','unpool/indices');
layers = connectLayers(layers,'maxPool/size','unpool/size');
                        
deepNetworkDesigner(layers);

%deepNetworkDesigner(segnetLayers([720 960 3], 11, 1))

opts = trainingOptions('sgdm', ...
    'InitialLearnRate',1e-3, ...
    'MaxEpochs',25, ...
    'MiniBatchSize',2, ...
    'Verbose',false, ...
	'Plots','training-progress');

%Создайте пиксельный datastore метки изображений, который содержит данные тренировки.
trainingData = pixelLabelImageDatastore(imds,pxds);

%net = trainNetwork(trainingData, layers, opts);


