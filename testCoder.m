sizeImage = [720, 960];
in = imread('peppers.png');
in = imresize(in, sizeImage); 

cfg = coder.config('DLL');
cfg.TargetLang = 'C++';
cfg.DeepLearningConfig = coder.DeepLearningConfig('cudnn'); 
codegen -args {ones(720,960,3,'uint8')} -config cfg myNET -report;



