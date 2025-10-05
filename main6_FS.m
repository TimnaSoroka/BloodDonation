% this is the same as main5_FS but using alredy prepared
% "Donors_table"/"Control_table"

clear
rng(500)

%% donors option
%load('Donors_table.mat');
%% control option
load('Controls_table.mat');
T=T1;
%%

load('Fields.mat');

thresh=mean([T.weight],'omitmissing'); %70.9622
sex_vec=([T.weight]<thresh);

X=T.X;
%Y=T.isDonorAfter;
Y=T.session;


X_train=X(~sex_vec,:);
Y_train=Y(~sex_vec);

X_test=X(sex_vec,:);
Y_test=Y(sex_vec);

% %% Descriptives
% if disp_ttests
% 
% for v=1:size(X,2)
% test_values= X(Y==1,v);
% retest_values = X(Y==2,v);
%     mx = mean(test_values, 'omitnan'); sx = std(test_values, 'omitnan'); nx = numel(test_values);
%     my = mean(retest_values, 'omitnan'); sy = std(retest_values, 'omitnan'); ny = numel(retest_values);
% 
%     valid = ~isnan(test_values) & ~isnan(retest_values);
%         x = test_values(valid);
%         y = retest_values(valid);
%         [~,p,~,stats] = ttest(x, y);
%         % Effect size: Cohen's dz for paired (mean diff / SD diff)
%             df = stats.df; tval = stats.tstat;
% varName=fields{v};
%     fprintf('%s: t(%d)=%.2f, p=%.4g]\n', ...
%         varName, df, tval, p);
%     fprintf('   before: %0.3f \xB1 %0.3f (n=%d);  after: %0.3f \xB1 %0.3f (n=%d)\n', ...
%         mx, sx, nx, my, sy, ny);
% end  
% end
%% FS
n=size(X_train,1)/2;
for i=1:size(X_train,2)
    currentfield=fields{i};
    test_values = X_train(1:n,i);  
retest_values = X_train(n+1:end,i);

[~,tp_value24(i),~,tstat24(i)] = ttest2(test_values, retest_values);
      %  fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])

end

%[a,idx]=sort(tp_value24);

top_features_idx = tp_value24<0.05; % example, replace with your indices

X_train_sel = X_train(:, top_features_idx);
X_test_sel  = X_test(:, top_features_idx);

% % %% scatter3
% % to_plot=1;
% % if to_plot
% % 
% % cond=Y_test==1;
% % figure('Color',[1 1 1]);
% % scatter3( ...
% %     X_test(cond, idx(1)), ...
% %     X_test(cond, idx(2)), ...
% %     X_test(cond, idx(3)), ...
% %     100, 'r', 'filled'); % group 1 (condition true)
% % hold on;
% % 
% % scatter3( ...
% %     X_test(~cond, idx(1)), ...
% %     X_test(~cond, idx(2)), ...
% %     X_test(~cond, idx(3)), ...
% %     100, 'b', 'filled'); % group 2 (condition false)
% % 
% % xlabel(vars{idx(1)},'FontSize',15);
% % ylabel(vars{idx(2)},'FontSize',15);
% % zlabel(vars{idx(3)},'FontSize',15);
% % %legend({'No blood loss','Blood loss'});
% % %view(45,30); % nice viewing angle
% % 
% % end

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
                mdl = fitctree(Xtr,Ytr,    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 4, ...
    'Surrogate', 'off');
            case 'Logistic'
                mdl = fitclinear(Xtr,Ytr,'Learner','logistic');
        end
        
        % Predict on validation fold
        Ypred = predict(mdl,Xval);
        CM = confusionmat(Yval, Ypred);
        acc_fold(k) =  sum(diag(CM))/length(Yval);

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

        cm = confusionmat(Y_test, y_pred_test);
        accuracY_test =  sum(diag(cm))/length(Y_test);

disp(['Best classifier: ', best_classifier_name]);
disp('Confusion Matrix:');
disp(cm);
%confusionchart(Y_test, y_pred_test);
disp(['Test Accuracy: ', num2str(accuracY_test)]);


% Compute ROC
% Define the positive class value exactly as in your labels:
    posClass = 'after';  % <-- change if your positive class is, e.g., true/'positive'/categorical('1')

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

