close all
clear
rng(5);

norm=1;
IntLength=5;
to_plot=1;
Fs=25; 
load('Holter_timings.mat');
%%
%sex_vec=logical([subjData.sex]);
%%
for i=1:size(subjData,2)
 [before{i},after{i},donation{i}]=extract_timings_needle(i,norm, IntLength,subjData);
%   [before{i},after{i},donation{i}]=extract_timings_15(i,norm, IntLength,subjData);
 %  [before{i},after{i},donation{i}]=extract_timings_15(i,norm, IntLength,subjData);

 %[NCbefore{i},NCafter{i},NCdonation{i}]=NC_analysis(i, IntLength,subjData);
end
before([2,42,63,65])=[];
after([2,42,63,65])=[];

X = [before,after];  % Each is [22500×1]


for i = 1:numel(X)
    if size(X{i}, 2) == 1
     %   plot(X{i})
        X{i} = X{i}';  % make sure it's [22500 × 1]
    end
end

xxx=randperm(size(X,2)/2,size(X,2)/2);
Y = categorical([zeros(1,size(before,2)), ones(1,size(after,2))]);  % 0 = pre, 1 = post


Y=Y';

inputSize = 1;
numHiddenUnits = 64;
numClasses = 2;

layers = [
    sequenceInputLayer(inputSize)

    convolution1dLayer(5, 32, 'Padding', 'causal', 'DilationFactor', 1)
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.2)

    convolution1dLayer(5, 64, 'Padding', 'causal', 'DilationFactor', 2)
    batchNormalizationLayer
    reluLayer
    dropoutLayer(0.2)


    globalAveragePooling1dLayer

    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs', 30, ...
    'MiniBatchSize', 8, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'Verbose', false);


    net = trainNetwork(X, Y, layers, options);

    YPred = classify(net, X);
accuracy = mean(YPred == Y);
fprintf('Accuracy: %.2f%%\n', accuracy * 100);

confusionchart(Y,YPred)