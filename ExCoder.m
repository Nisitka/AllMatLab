in = imread('peppers.png');
in = imresize(in,[224,224]);

cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.DeepLearningConfig = coder.DeepLearningConfig('mkldnn'); 
codegen -args {ones(224,224,3,'uint8')} -config cfg myVGG16 -report;

predict_scores = myVGG16_mex(in);
[scores,indx] = sort(predict_scores, 'descend');
net = coder.loadDeepLearningNetwork('vgg16.mat');
classNames = net.Layers(end).Classes;
disp(classNames(indx(1:5))); 