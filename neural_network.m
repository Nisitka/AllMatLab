function neural_network 
path = 'C:\Program Files\Polyspace\R2020b\toolbox\nnet\nndemos\nndatasets\DigitDataset';

imds = imageDatastore(path, 'IncludeSubfolders', true, 'LabelSource','foldernames');
imds.countEachLabel

[train, test] = imds.splitEachLabel(0.6, 'randomize');


end