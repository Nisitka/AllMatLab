function out = myNET(in) %#codegen

persistent mynet;

if isempty(mynet)
    mynet = coder.loadDeepLearningNetwork('sigNet.mat', 'myNET');
end

out = predict(mynet,in);   