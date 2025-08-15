close all
clear
rng(1000)
norm=1; %1=80%, 0=50%
IntLength=5;
to_plot=0;
Fs=25; 
load('Holter_timings.mat');
%%
subjData(91)=[]; %have short after (*technical issue)

%%
for i=1:size(subjData,2)
 [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData);
end

 %%



  noiseThreshold=0;

 for i=1:size(NCbefore,2)
     if ismember(i,[42,63,65]) 
         Fs=6;
     else
         Fs=25;
     end
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters(NCbefore{i},Fs,noiseThreshold);
[Laterality_IndexA(i),AmeasureResults(i)]=NasalCycleParameters(NCafter{i},Fs,noiseThreshold);
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
 end

  for i=1:size(NCbefore,2)
      if ~isempty([Laterality_IndexA(i).one])
 LI_ampB(i,:)=abs([Laterality_IndexB(i).one]);
 LI_B(i,:)=[Laterality_IndexB(i).one];
 LI_ampA(i,1:length(abs([Laterality_IndexA(i).one])))=abs([Laterality_IndexA(i).one]);
 LI_A(i,1:1:length(abs([Laterality_IndexA(i).one])))=[Laterality_IndexA(i).one];
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

vals_before=calculate_before_after(before,IntLength);
[vals_after,vars]=calculate_before_after(after,IntLength);

fields=fieldnames(vals_after);

vals_before([42,63,65])=[];
vals_after([42,63,65])=[];

X=table2array(struct2table([vals_before,vals_after]));

BmeasureResults([42,63,65])=[];
AmeasureResults([42,63,65])=[];
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
x=x(:,[1,4,6,9]);

X=[X,x];
Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];


sex_vec1=logical([subjData.sex]);
sex_vec1([42,63,65])=[];
sex_vec=[sex_vec1,sex_vec1];

X_train=X(~sex_vec,:);
Y_train=Y(~sex_vec);

X_test=X(sex_vec,:);
Y_test=Y(sex_vec);

%% FS
n=size(X_train,1)/2;
for i=1:size(fields,1)
    currentfield=fields{i};
    test_values = X_train(1:n,i);  
retest_values = X_train(n+1:end,i);

[~,tp_value24(i),~,tstat24(i)] = ttest2(test_values, retest_values);

                % Descriptives
    mx = mean(test_values, 'omitnan'); sx = std(test_values, 'omitnan'); nx = numel(test_values);
    my = mean(retest_values, 'omitnan'); sy = std(retest_values, 'omitnan'); ny = numel(retest_values);

    valid = ~isnan(test_values) & ~isnan(retest_values);
        x = test_values(valid);
        y = retest_values(valid);
        [~,p,~,stats] = ttest2(x, y);
        % Effect size: Cohen's dz for paired (mean diff / SD diff)
            df = stats.df; tval = stats.tstat;
varName=currentfield;
    fprintf('%s: t(%d)=%.2f, p=%.4g]\n', ...
        varName, df, tval, p);
    fprintf('   before: %0.3f \xB1 %0.3f (n=%d);  after: %0.3f \xB1 %0.3f (n=%d)\n', ...
        mx, sx, nx, my, sy, ny);
                      %  fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])

end

%% FS

features = fscnca(X_train, Y_train);  % fscnca selects features based on correlation with the target
[f_v,f_idx]=sort(features.FeatureWeights,'descend');
[a,idx]=sort(tp_value24);

top_features_idx = tp_value24<0.05; % example, replace with your indices

X_train_sel = X_train(:, top_features_idx);
X_test_sel  = X_test(:, top_features_idx);

%% --- Prepare Data ---

% Impute missing values (median of training set)
for f = 1:size(X_train_sel,2)
    nan_idx_train = isnan(X_train_sel(:,f));
    X_train_sel(nan_idx_train,f) = median(X_train_sel(:,f), 'omitnan');
    
    nan_idx_test = isnan(X_test_sel(:,f));
    X_test_sel(nan_idx_test,f) = median(X_train_sel(:,f), 'omitnan'); % use training median
end

%% --- Classifier Definitions ---
K = 5;
cv = cvpartition(Y_train,'KFold',K);

%classifier_names = {'Logistic','kNN'};
classifier_names = {'Logistic','kNN','DecisionTree','SVM-linear','SVM-rbf'};

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
                mdl = fitctree(Xtr,Ytr);
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
y_pred_test = predict(best_clf, X_test_sel);
accuracY_test = sum(y_pred_test == Y_test) / length(Y_test);
cm = confusionmat(Y_test, y_pred_test);

disp(['Best classifier: ', best_classifier_name]);
disp('Confusion Matrix:');
disp(cm);
disp(['Test Accuracy: ', num2str(accuracY_test)]);


%% Compute ROC
% [fpRate, tpRate, ~, AUC] = perfcurve(labels, y_score, 1);
% 
% % Plot ROC curve
% figure;
% plot(fpRate, tpRate, 'b-', 'LineWidth', 2); hold on;
% plot([0 1], [0 1], 'k--');
% xlabel('False Positive Rate (1 - Specificity)');
% ylabel('True Positive Rate (Sensitivity)');
% title(sprintf('LOOCV ROC - Coarse Tree (AUC = %.2f)', AUC));
% grid on;
% axis square;


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

