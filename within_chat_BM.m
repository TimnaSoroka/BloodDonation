clear; close all; rng(10);


Fs = 25;  % sampling rate
binMinutes = 1.5;  % 2-minute bins
%winLen = Fs * winLenSec;  % 3000 samples
% overlap = 0.5;  % 50% overlap
% hopLen = round(winLen * (1 - overlap));
%hopLen = winLen / 2;

load('Holter_timings.mat');  % assumes `subjData` is available
norm = 1; 
IntLength = 10;

%subjData([subjData.Weight]>72)=[];

totalSubjects = numel(subjData);
subjectTrue = [];
subjectPred = [];

for sbj = 1:totalSubjects
    if ismember(sbj, [ 29    38]) %    25    29    38 all;[16,19,28]>80,[16,24]>75
        continue
    end


  %[before, after, ~] = extract_timings_needle_walk_in_chair(sbj, norm, IntLength, subjData); %_walk_in_chair
  Data=subjData(sbj).Data;
  D=table2array(Data(:,[3,4]));
  sum_resp=sum(D,2);
resp=zscore(sum_resp);

if numel(resp)/(60*25)<40
    Fs=6;
else
    Fs=25;
end

if ~isempty(subjData(sbj).walk)
before=resp(1:subjData(sbj).walk);
else
    before=resp(1:IntLength*Fs*60);
end

if ~isempty(subjData(sbj).in_chair)
after=resp(subjData(sbj).in_chair:end);
else
    after=resp(subjData(sbj).out+3*60*Fs:end);
end


    % Convert to row vector
    if size(before,2) == 1, before = before'; end
    if size(after,2) == 1, after = after'; end

    % % Sliding windows
    % binsBefore = sliding_window(before, winLen, hopLen);
    % binsAfter = sliding_window(after, winLen, hopLen);
L=60*Fs*binMinutes;



binsBefore = mat2cell(before,1, [repmat(L,1,floor(numel(before)/L)) rem(numel(before),L)]);
if isempty(binsBefore{end}), binsBefore(end)=[]; end     % drop empty tail if exact multiple

binsAfter = mat2cell(after,1, [repmat(L,1,floor(numel(after)/L)) rem(numel(after),L)]);
if isempty(binsAfter{end}), binsAfter(end)=[]; end     % drop empty tail if exact multiple

minLen=60*Fs;
binsBefore = binsBefore(cellfun(@numel, binsBefore) >= minLen);
binsAfter = binsAfter(cellfun(@numel, binsAfter) >= minLen);

    % Combine bins
    X = [binsBefore, binsAfter];
        N = numel(X);

        y_score = nan(N,1);   % score/probability for 'after' class
y_true  = nan(N,1);   % true labels for the test bins (categorical)

    %X_BM = cellfun(@breathmetrics_feats, X, 'UniformOutput', false);
X_BM = cellfun(@(x) breathmetrics_feats(x, Fs), X, 'UniformOutput', false);
        S = vertcat(X_BM{:}); 
        T=struct2table(S);

        m = varfun(@(x) mean(ismissing(x)), T, 'OutputFormat','uniform');
%T(:, m > 0.3) = [];
X=table2array(T);

    Y = categorical([zeros(numel(binsBefore),1); ones(numel(binsAfter),1)]);  % 0 = before, 1 = after

    % Leave-One-Bin-Out CV
    binPred = zeros(N,1);

    for i = 1:N
        trainIdx = setdiff(1:N, i);
        testIdx = i;

        X_train = X(trainIdx,:);
        Y_train = Y(trainIdx);
        X_test  = X(testIdx,:);
        Y_test  = Y(testIdx);

        %                 mdl = fitclinear(X_train,Y_train,'Learner','logistic');
[h,p] = ttest2(X_train(Y_train=="0",:),X_train(Y_train=="1",:));  % one-way ANOVA, no figure


% if sum(h)<5
   [~, ord] = sort(p, 'ascend');
    sel = ord(1:5);
% else
%     sel=find(h);
% end

% 
%         template = templateTree(...
%     'MaxNumSplits', 8, ...
%     'NumVariablesToSample', 'all');
% mdl = fitcensemble(...
%     X_train(:,sel), ...
%     Y_train, ...
%     'Method', 'Bag', ...
%     'NumLearningCycles', 30, ...
%     'Learners', template, ...
%     'ClassNames', [0,1]);

subspaceDimension = max(1, min(2, width(X_train(:,sel)) - 1));
mdl = fitcensemble(...
    X_train(:,sel), ...
    Y_train, ...
    'Method', 'Subspace', ...
    'NumLearningCycles', 30, ...
    'Learners', 'knn', ...
    'NPredToSample', subspaceDimension, ...
    'ClassNames', [0,1]);

    % y_pred(i,sbj) = predict(mdl,X_test(:,sel));
[y_pred(i), score] = predict(mdl, X_test(:,sel));               % score is 1xK
posIdx = find(mdl.ClassNames == 1);       
y_score(i) = score(1,posIdx);                          % scalar score for this test bin
end

% Test on held-out test set
    cm=confusionmat(Y,categorical(y_pred));
    acc(sbj)=sum(diag(cm))./sum(sum(cm));
    spec(sbj)=cm(2,2)/(cm(2,1)+cm(2,2));
    [FPR, TPR, ~, AUC(sbj)] = perfcurve(Y, y_score, categorical(1));
    % Vote for subject label

    clearvars y_pred
end

% Final accuracy
finalAcc = mean(acc);
finalAUC = mean(AUC);

fprintf('✅ Within-subject LOO accuracy: %.2f%%\n', finalAcc * 100);






function [segments] = sliding_window(sig, winLen, hopLen)
    segments = {};
    N = length(sig);
    for startIdx = 1:hopLen:(N - winLen + 1)
        segments{end+1,1} = sig(startIdx:startIdx + winLen - 1);
    end
end

% function feats = breathmetrics_feats(x, Fs)
% % BREATHMETRICS_FEATS
% % Compute per-window BreathMetrics secondary features for an airflow signal x.
% % x: vector (airflow), Fs: sampling rate (Hz).
% 
%     x = x(:);  % column
% 
%     % Build BreathMetrics object (human airflow mode)
%     bmObj = breathmetrics(x, Fs, 'humanAirflow');
% 
%     % Compute all features; sliding=1, plotting=0 (matches your snippet)
%     bmObj.estimateAllFeatures(0, 'simple', 1, 0);
% 
%     % Desired secondary-feature names (same order you provided)
%     variableNames = { ...
%         'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', ...
%         'AverageInhaleDuration', 'AverageInhalePauseDuration', 'AverageInhaleVolume', ...
%         'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', 'AveragePeakInspiratoryFlow', ...
%         'AverageTidalVolume', 'BreathingRate', ...
%         'CoefficientOfVariationOfBreathVolumes', 'CoefficientOfVariationOfBreathingRate', ...
%         'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
%         'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', ...
%         'DutyCycleOfExhale', 'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', ...
%         'MinuteVentilation', 'PercentOfBrethsWithExhalePause', 'PercentOfBrethsWithInhalePause' ...
%     };
% 
%     % Pull values from BreathMetrics (secondary features)
%     vals = bmObj.secondaryFeatures.values;
% 
%     % Some BreathMetrics versions return a cell array, others a numeric row vector.
%     if iscell(vals)
%         vals = cellfun(@(v) double(v), vals);
%     else
%         vals = double(vals);
%     end
% 
%     % Assign to struct fields (missing entries become NaN if needed)
%     feats = struct();
%     for ii = 1:numel(variableNames)
%         if ii <= numel(vals)
%             feats.(variableNames{ii}) = vals(ii);
%         else
%             feats.(variableNames{ii}) = NaN;
%         end
%     end
% end

function feats = breathmetrics_feats(x, Fs)
% BREATHMETRICS_FEATS
% Compute per-window BreathMetrics secondary features for an airflow signal x.
% x: vector (airflow), Fs: sampling rate (Hz).

    x = x(:);  % column

         peaks=peaks_from_ts2(x,Fs);
feats= calculate_zz(peaks);
feats=rmfield(feats,'COV_ExhaleVolume');
    % Build BreathMetrics object (human airflow mode)

    % Desired secondary-feature names (same order you provided)
    % variableNames = { ...
    %     'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', ...
    %     'AverageInhaleDuration', 'AverageInhalePauseDuration', 'AverageInhaleVolume', ...
    %     'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', 'AveragePeakInspiratoryFlow', ...
    %     'AverageTidalVolume', 'BreathingRate', ...
    %     'CoefficientOfVariationOfBreathVolumes', 'CoefficientOfVariationOfBreathingRate', ...
    %     'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    %     'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', ...
    %     'DutyCycleOfExhale', 'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', ...
    %     'MinuteVentilation', 'PercentOfBrethsWithExhalePause', 'PercentOfBrethsWithInhalePause' ...
    % };


   % 
    % % Pull values from BreathMetrics (secondary features)
    % vals = bmObj.secondaryFeatures.values;

    % Some BreathMetrics versions return a cell array, others a numeric row vector.
    % if iscell(vals)
    %     vals = cellfun(@(v) double(v), vals);
    % else
    %     vals = double(vals);
    % end

    % Assign to struct fields (missing entries become NaN if needed)
    % feats = struct();
    % for ii = 1:numel(variableNames)
    %     if ii <= numel(vals)
    %         feats.(variableNames{ii}) = vals(ii);
    %     else
    %         feats.(variableNames{ii}) = NaN;
    %     end
    % end
end