cfg = coder.gpuConfig('dll');
cfg.TargetLang = 'C++';
cfg.DeepLearningConfig = coder.DeepLearningConfig('cudnn');
cfg.DeepLearningConfig.AutoTuning = true;
cfg.DeepLearningConfig.DataType = 'fp32';

codegen -config cfg myNET -args {ones(720,960,3,'uint8')} -report