% This code use the men for training and women for testing or fat for training and less fat for testing. FS done using
% values<0.05 in the man sample.

close all
clear
rng(50)
norm=1; %1=80%, 0=50%
IntLength=5; % 3 is the best
to_plot=1;
Fs=25; 
%
 load('Holter_timings.mat');

%load('Holter_timings_controls.mat');

disp_ttests=0;

%%
for i=1:size(subjData,2)
 [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData); %_needle_walk_in_chair
%     [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle_walk_in_chair(i,norm, IntLength,subjData); %_needle_walk_in_chair
end

 %%
keep = cellfun(@(v) size(v,1) >= IntLength*60*Fs, NCafter);
removedIdx = find(~keep);              % indices of cells you dropped
disp(removedIdx);
NCafter = NCafter(keep);                            % keep only long-enough cells
NCbefore = NCbefore(keep);  
before = before(keep);  
after = after(keep);  


  noiseThreshold=0;

 for i=1:size(NCbefore,2)
     % if ismember(i,[42,63,65]) 
     %     Fs=6;
     % else
     %     Fs=25;
     % end
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters_short(NCbefore{i},Fs,noiseThreshold,'Holter');
[Laterality_IndexA(i),AmeasureResults(i)]=NasalCycleParameters_short(NCafter{i},Fs,noiseThreshold,'Holter');
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
 end

  for i=1:size(NCbefore,2)
      if ~isempty([Laterality_IndexA(i).two])
 LI_ampB{i}=abs([Laterality_IndexB(i).two]);
 LI_B{i}=[Laterality_IndexB(i).two];
 LI_ampA{i}=abs([Laterality_IndexA(i).two]);
 LI_A{i}=[Laterality_IndexA(i).two];
      end
  end

%%

 fields=fieldnames(BmeasureResults);
for i=1:size(fields,1)
    currentfield=fields{i};
test_values = [BmeasureResults(:).(currentfield)];   
retest_values = [AmeasureResults(:).(currentfield)]; 

%[p_values(i),~,Wstat(i)] = signrank(test_values, retest_values,"method","approximate");
[~,p_values_ttest(i)] = ttest2(test_values, retest_values);

% fprintf('%s ttest p = %.2f\n',currentfield, p_values_ttest(i))
%         fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])
end

%% big differences

vals_before=calculate_before_after(before);
[vals_after,vars]=calculate_before_after(after);

fields2=fieldnames(vals_after);

fields=[fields2;fields{1};fields{3};fields{4};fields{6}];

%save('Fields.mat',"fields");

% vals_before([42,63,65])=[];
% vals_after([42,63,65])=[];

X=table2array(struct2table([vals_before,vals_after]));

% BmeasureResults([42,63,65])=[];
% AmeasureResults([42,63,65])=[];
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
x=x(:,[1,3,4,6]);

X=[X,x];

Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];


thresh=mean([subjData.Weight],'omitmissing');
sex_vec1=([subjData.Weight]<thresh);
sex_vec1=sex_vec1(keep);
% sex_vec1=logical([subjData.sex]);

% sex_vec1([42,63,65])=[];
sex_vec=[sex_vec1,sex_vec1];

X_train=X(~sex_vec,:);
Y_train=Y(~sex_vec);

X_test=X(sex_vec,:);
Y_test=Y(sex_vec);

%% Descriptives
if disp_ttests
   
for v=1:size(X,2)
test_values= X(Y==1,v);
retest_values = X(Y==2,v);
    mx = mean(test_values, 'omitnan'); sx = std(test_values, 'omitnan'); nx = numel(test_values);
    my = mean(retest_values, 'omitnan'); sy = std(retest_values, 'omitnan'); ny = numel(retest_values);

    valid = ~isnan(test_values) & ~isnan(retest_values);
        x = test_values(valid);
        y = retest_values(valid);
        [~,p,~,stats] = ttest(x, y);
        % Effect size: Cohen's dz for paired (mean diff / SD diff)
            df = stats.df; tval = stats.tstat;
varName=fields{v};
    fprintf('%s: t(%d)=%.2f, p=%.4g]\n', ...
        varName, df, tval, p);
    fprintf('   before: %0.3f \xB1 %0.3f (n=%d);  after: %0.3f \xB1 %0.3f (n=%d)\n', ...
        mx, sx, nx, my, sy, ny);

    
end  
end
%% FS
n=size(X_train,1)/2;
for i=1:size(X_train,2)
    currentfield=fields{i};
    test_values = X_train(1:n,i);  
retest_values = X_train(n+1:end,i);

[~,tp_value24(i),~,tstat24(i)] = ttest2(test_values, retest_values);
      %  fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])

end

features = fscnca(X_train, Y_train);  % fscnca selects features based on correlation with the target
[f_v,f_idx]=sort(features.FeatureWeights,'descend');

% top_features_idx = tp_value24<0.05; % example, replace with your indices
   
[a,idx]=sort(tp_value24);
  top_features_idx=idx(1:5);

X_train_sel = X_train(:, top_features_idx);
X_test_sel  = X_test(:, top_features_idx);

%% scatter3
to_plot=1;
if to_plot

cond=Y_test==1;
figure('Color',[1 1 1]);
scatter3( ...
    X_test(cond, idx(1)), ...
    X_test(cond, idx(2)), ...
    X_test(cond, idx(3)), ...
    100, 'r', 'filled'); % group 1 (condition true)
hold on;

scatter3( ...
    X_test(~cond, idx(1)), ...
    X_test(~cond, idx(2)), ...
    X_test(~cond, idx(3)), ...
    100, 'b', 'filled'); % group 2 (condition false)

xlabel(fields{idx(1)},'FontSize',15);
ylabel(fields{idx(2)},'FontSize',15);
zlabel(fields{idx(3)},'FontSize',15);
legend({'No blood loss','Blood loss'});
view(45,30); % nice viewing angle


med = nanmedian(X_test,1);
Ximp = X_test;
for j = 1:size(Ximp,2), m = isnan(Ximp(:,j)); Ximp(m,j) = med(j); end

Xz=zscore(Ximp);
[coeff,score,~,~,expl] = pca(Xz);   % score = PCs
figure; scatter3(score(Y_test==2,1),score(Y_test==2,2),score(Y_test==2,3),100,'b','filled'); hold on
scatter3(score(Y_test==1,1),score(Y_test==1,2),score(Y_test==1,3),100,'r','filled'); grid on
xlabel(sprintf('PC1 (%.1f%%)',expl(1))); ylabel(sprintf('PC2 (%.1f%%)',expl(2))); zlabel(sprintf('PC3 (%.1f%%)',expl(3)));
legend('Blood loss','No blood loss'); 

end

%% --- Prepare Data ---

% Impute missing values (median of training set)
% for f = 1:size(X_train_sel,2)
%     nan_idx_train = isnan(X_train_sel(:,f));
%     X_train_sel(nan_idx_train,f) = median(X_train_sel(:,f), 'omitnan');
% 
%     nan_idx_test = isnan(X_test_sel(:,f));
%     X_test_sel(nan_idx_test,f) = median(X_train_sel(:,f), 'omitnan'); % use training median
% end

%save('Features_selected_from_training_set.mat','top_features_idx');
%% --- Classifier Definitions ---
K = 5;
cv = cvpartition(Y_train,'KFold',K);

classifier_names = {'Logistic'};
%classifier_names = {'Logistic','kNN','DecisionTree','SVM-linear','SVM-rbf'};

cv_accuracy = zeros(length(classifier_names),1);

for c = 1:length(classifier_names)
    acc_fold = zeros(K,1);  % accuracy per fold
    
    for k = 1:K
        trIdx = training(cv,k);
        valIdx = test(cv,k);
        
        Xtr = X_train_sel(trIdx,:);
        Ytr = Y_train(trIdx);
        Xval = X_train_sel(valIdx,:);
        Yval = Y_train(valIdx);
        
        % Train classifier
        switch classifier_names{c}
            case 'SVM-linear'
                mdl = fitcsvm(Xtr,Ytr,'KernelFunction','linear');
            case 'SVM-rbf'
                mdl = fitcsvm(Xtr,Ytr,'KernelFunction','rbf');
            case 'kNN'
                mdl = fitcknn(Xtr,Ytr,'NumNeighbors',5);
            case 'DecisionTree'
                mdl = fitctree(Xtr,Ytr,    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 4, ...
    'Surrogate', 'off');
            case 'Logistic'
                mdl = fitclinear(Xtr,Ytr,'Learner','logistic');
        end
        
        % Predict on validation fold
        Ypred = predict(mdl,Xval);
        acc_fold(k) = sum(Ypred == Yval)/length(Yval);
    end
    
    cv_accuracy(c) = mean(acc_fold); % mean CV accuracy
end

% Display CV results
cv_table = table(classifier_names', cv_accuracy, 'VariableNames', {'Classifier','CV_Accuracy'});
disp(cv_table);

% Pick best classifier
[~, best_idx] = max(cv_accuracy);
best_classifier_name = classifier_names{best_idx};

% Train on full training set
switch best_classifier_name
    case 'SVM-linear'
        best_clf = fitcsvm(X_train_sel, Y_train, 'KernelFunction','linear');
    case 'SVM-gaussian'
        best_clf = fitcsvm(X_train_sel, Y_train, 'KernelFunction','gaussian');
    case 'kNN'
        best_clf = fitcknn(X_train_sel, Y_train, 'NumNeighbors',5);
    case 'DecisionTree'
        best_clf = fitctree(X_train_sel, Y_train);
    case 'Logistic'
        best_clf = fitclinear(X_train_sel, Y_train,'Learner','logistic');
end

% Test on held-out test set
[y_pred_test,y_score] = predict(best_clf, X_test_sel);
accuracY_test = sum(y_pred_test == Y_test) / length(Y_test);
cm = confusionmat(Y_test, y_pred_test);

disp(['Best classifier: ', best_classifier_name]);
disp('Confusion Matrix:');
disp(cm);
disp(['Test Accuracy: ', num2str(accuracY_test)]);


% Compute ROC
% Define the positive class value exactly as in your labels:
    posClass = 2;  % <-- change if your positive class is, e.g., true/'positive'/categorical('1')

    % Find the score column corresponding to posClass
    cls = best_clf.ClassNames;            % same order as columns in y_score
    idx = find(ismember(cls, posClass));  % works for numeric, logical, char, categorical

    if isempty(idx)
        error('Positive class not found in ClassNames. Check posClass.');
    end

    % Scores for the positive class
    posScores = y_score(:, idx);

    % ROC curve & AUC
    [fpRate, tpRate, ~, AUC] = perfcurve(Y_test, posScores, posClass);
    fprintf('AUC: %.4f\n', AUC);
    
% Plot ROC curve
figure;
plot(fpRate, tpRate, 'b-', 'LineWidth', 2); hold on;
plot([0 1], [0 1], 'k--');
xlabel('False Positive Rate (1 - Specificity)');
ylabel('True Positive Rate (Sensitivity)');
title(sprintf('LOOCV ROC - Coarse Tree (AUC = %.2f)', AUC));
grid on;
axis square;


%% correlations during donation

% for i=1:size(donation,2)
%     full_don=donation{i};
%     ints=7;
%     vals_during{i}=calculate_1min_bin(full_don,ints);  
% end
% 
% vals_during_conc = [vals_during{:}];	
% time=[repmat(1:ints,1,size(vals_during,2))];
% 
% for i=1:size(fields,1)
%     current_field=fields{i};
%         current_vals=[vals_during_conc.(current_field)];
%         [r_spearman, p_spearman] = corr(time', current_vals', 'Type', 'Spearman','rows','complete'); % Spearman correlation
%         if p_spearman<0.1
%                 figure
%             hold on
%         scatter(time, current_vals)
%         lsline
%             title([current_field ' p=' num2str(p_spearman)])
%             set(gca,'FontSize',12)
% 
%         end
%     end



% function z_breath_values=calculate_1min_bin(before1,num_bins)
% 
% if isinteger(length(before1)/num_bins)
% binned_data = reshape(before1, [], num_bins); % Create a 150x10 matrix
% else
%     N = length(before1); % Get vector length
%     remainder = mod(N, num_bins); % Find remainder when divided by 10
% 
%     % Remove the first 'remainder' elements
%     trimmed_vector = before1(remainder + 1:end); 
% 
%     % Compute new length
%     new_N = length(trimmed_vector);
% 
%     % Reshape into a 10-row matrix (each column is a bin)
%     binned_data = reshape(trimmed_vector, new_N / num_bins, num_bins);
% end
% 
% for i=1:num_bins
%  if ismember(i,[42,63,65]) 
%          Fs=6;
%      else
%          Fs=25;
%      end
% if Fs<25
%     % peaks1 = peaks_from_ts_fs(binned_data(:,i),Fs);
%     % 
%     % plot(before{i})
%     % hold on
%     % plot([peaks1.PeakLocation],[peaks1.PeakValue],'bo')
%     % 
%     % z_breath_values(i) = calculate_z_blood(peaks1);
%     continue
% else
%      bmObj=breathmetrics(binned_data(:,i),Fs,'humanAirflow');
%  bmObj.estimateAllFeatures();
% 
%  % Define the variable names as field names in the struct
% variableNames = {...
%     'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
%     'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
%     'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
%     'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
%     'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
%     'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
%     'PercentOfBreathsWithInhalePause'};
% 
% values=[bmObj.secondaryFeatures.values];
% % Assign the values to the struct fields
% for ii = 1:length(variableNames)
%     dataStruct.(variableNames{ii}) = values{ii};
% end
% 
%  z_breath_values(i)= dataStruct(:);
% end
% 
% end
% end

