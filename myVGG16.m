function out = myVGG16(in)  %объявление функции

persistent mynet;  %переменная существующая только в этой функции 
if isempty(mynet)  %если переменная не инцилизирована 
    mynet = coder.loadDeepLearningNetwork('vgg16.mat', 'myVGGnet');
end

out = predict(mynet,in);
