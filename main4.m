close all
clear
rng(1000)
norm=1;
IntLength=5;
to_plot=0;
Fs=25; 
load('Holter_timings.mat');
%%
subjData(91)=[]; %have short after (*technical issue)

sex_vec=logical([subjData.sex]);

% subjData=subjData(~sex_vec);
%%
for i=1:size(subjData,2)
 [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData);
 %[NCbefore{i},NCafter{i},NCdonation{i}]=NC_analysis(i, IntLength,subjData);
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
delta=test_values-retest_values;
outlierMask = isoutlier(delta,"mean");      % Detect outliers 
retest_values = retest_values(~outlierMask);   % Remove outliers
test_values = test_values(~outlierMask);   % Remove outliers

%[p_values(i),~,Wstat(i)] = signrank(test_values, retest_values,"method","approximate");
[~,p_values_ttest(i)] = ttest(test_values, retest_values);
%%
% sex_vec_clean=sex_vec(~outlierMask);
%  [~,p_values_ttest_women(i)] = ttest(test_values(sex_vec_clean), retest_values(sex_vec_clean));
%   [~,p_values_ttest_men(i)] = ttest(test_values(~sex_vec_clean), retest_values(~sex_vec_clean));

%%
if to_plot
    figure;
min_=min([test_values,retest_values]);
max_=max([test_values,retest_values]);
 title(currentfield)
scatter(test_values, retest_values, 100, 'b', 'filled'); % Blue filled markers
hold on
plot([min_ max_],[min_ max_], 'k--'); % Identity line
xlim([min_ max_])
ylim([min_ max_])
set(gca,'FontSize',12)
% Perform Mann-Whitney U test (non-parametric test for paired samples)

% Format title with p-value
title(sprintf('%s (ttest p = %.2f)',currentfield, p_values_ttest(i)));

% Axis labels and limits
xlabel('Before');
ylabel('After');
        end
        fprintf('%s ttest p = %.2f\n',currentfield, p_values_ttest(i))
        fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])

end
% 


%% big differences

vals_before=calculate_before_after(before,IntLength);
[vals_after,vars]=calculate_before_after(after,IntLength);

vals_before([42,63,65])=[];
vals_after([42,63,65])=[];

X=table2array(struct2table([vals_before,vals_after]));

BmeasureResults([42,63,65])=[];
AmeasureResults([42,63,65])=[];
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
x=x(:,[1,4,6,9]);

%X=[X,x]; Adding NC parameters reduce accuracy classification
X=[X];
Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];


fields=fieldnames(vals_after);


%%
for i=1:size(fields,1)
    currentfield=fields{i};
    test_values = [vals_before(:).(currentfield)];   % Test scores for 3 subjects
retest_values = [vals_after(:).(currentfield)]; % Retest scores for the same 3 subjects

delta=test_values-retest_values;
outlierMask = isoutlier(delta,"mean");      % Detect outliers (default: median + MAD)
retest_values = retest_values(~outlierMask);   % Remove outliers
test_values = test_values(~outlierMask);   % Remove outliers

%[p_value24(i),~,Wstat24(i)] = signrank(test_values, retest_values,'method','approximate');
[~,tp_value24(i),~,tstat24(i)] = ttest2(test_values, retest_values);

 sex_vec_clean2=sex_vec(~outlierMask);
%  [~,tp_value24_women(i)] = ttest(test_values(sex_vec_clean2), retest_values(sex_vec_clean2));
%   [~,tp_value24_men(i)] = ttest(test_values(~sex_vec_clean2), retest_values(~sex_vec_clean2));

if tp_value24(i)<0.05
if to_plot
    figure;
min_=min([test_values,retest_values]);
max_=max([test_values,retest_values]);
scatter(test_values(sex_vec_clean2), retest_values(sex_vec_clean2), 100, 'r', 'filled'); % Blue filled markers
hold on
scatter(test_values(~sex_vec_clean2), retest_values(~sex_vec_clean2), 100, 'b', 'filled'); % Blue filled markers

plot([min_ max_],[min_ max_], 'k--'); % Identity line
xlim([min_ max_])
ylim([min_ max_])
% Perform Mann-Whitney U test (non-parametric test for paired samples)
set(gca,'FontSize',12)

% Format title with p-value
title(sprintf('%s (ttest p = %.4f)',vars{i},tp_value24(i)));

% Axis labels and limits
xlabel('Before');
ylabel('After');
        end
                % Descriptives
    mx = mean(test_values, 'omitnan'); sx = std(test_values, 'omitnan'); nx = numel(test_values);
    my = mean(retest_values, 'omitnan'); sy = std(retest_values, 'omitnan'); ny = numel(retest_values);

    valid = ~isnan(test_values) & ~isnan(retest_values);
        x = test_values(valid);
        y = retest_values(valid);
        [~,p,~,stats] = ttest(x, y);
        % Effect size: Cohen's dz for paired (mean diff / SD diff)
            df = stats.df; tval = stats.tstat;
varName=currentfield;
    fprintf('%s: t(%d)=%.2f, p=%.4g]\n', ...
        varName, df, tval, p);
    fprintf('   before: %0.3f \xB1 %0.3f (n=%d);  after: %0.3f \xB1 %0.3f (n=%d)\n', ...
        mx, sx, nx, my, sy, ny);
                      %  fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])
end
end

%save(['test_retest_' num2str(IntLength) '.mat'] ,'vals_before','vals_after','tp_value24');

%% classification


hiddenSize = 14;
autoenc = trainAutoencoder(X', hiddenSize, ...
    'MaxEpochs',400, 'L2WeightRegularization',0.001);
X_train_enc = encode(autoenc, X')';


%mrmr_order=fscmrmr(vars_to_class_converged,labels);
%vars_to_class=[vars_to_class(:,all)]; %,more_vars_to_class(:,[2,5,7])
% [~,acc]=fine_knn(vars_to_class,labels);
% [z]=ttest2(vars_to_class(labels==1,:),vars_to_class(labels==2,:));
% Z=[11,16,24];
%  %labels(isnan(vars_to_c     lass(1,:)))=[];
% % vars_to_class(:,isnan(vars_to_class(1,:)))=[];
% [labels_hat,acc]=L_SVM(vars_to_class(:,z==1),labels);
% [labels_hat,acc_tree,y_score] = tree(vars_to_class(:,z==1),labels);


classificationSVM = fitcsvm(...
    X_train_enc, ...
    Y, ...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'Standardize', true, ...
     'ClassNames', [1 2], ...
                   'ScoreTransform', 'none');  % important for ROC);

partitionedModel = crossval(classificationSVM, 'Leaveout','on');

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

acc_autoenc_validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% Compute ROC
% Define the positive class value exactly as in your labels:
    posClass = 2;  % <-- change if your positive class is, e.g., true/'positive'/categorical('1')

    % Find the score column corresponding to posClass
    cls = [1,2];            % same order as columns in y_score
    idx = find(ismember(cls, posClass));  % works for numeric, logical, char, categorical

    if isempty(idx)
        error('Positive class not found in ClassNames. Check posClass.');
    end

    % Scores for the positive class
    posScores = validationScores(:, idx);

    % Find rows with no NaNs in X_train_enc
validIdx = all(~isnan(X_train_enc), 2);

% Use only valid labels and scores
validLabels = Y(validIdx);

    % ROC curve & AUC
    [fpRate, tpRate, ~, AUC] = perfcurve(validLabels, posScores, posClass);
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

% [trainedClassifier, validationAccuracySVM] = gausianSVM(vars_to_class,labels);
% [trainedClassifier, validationAccuracyBT] = bagged_tree(vars_to_class,labels);

% cc=vars_to_class(:,zzz(4:13));

% figure
% confusionchart(labels,labels_hat)
% 




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



function z_breath_values=calculate_1min_bin(before1,num_bins)

if isinteger(length(before1)/num_bins)
binned_data = reshape(before1, [], num_bins); % Create a 150x10 matrix
else
    N = length(before1); % Get vector length
    remainder = mod(N, num_bins); % Find remainder when divided by 10

    % Remove the first 'remainder' elements
    trimmed_vector = before1(remainder + 1:end); 

    % Compute new length
    new_N = length(trimmed_vector);
    
    % Reshape into a 10-row matrix (each column is a bin)
    binned_data = reshape(trimmed_vector, new_N / num_bins, num_bins);
end

for i=1:num_bins
 if ismember(i,[42,63,65]) 
         Fs=6;
     else
         Fs=25;
     end
if Fs<25
    % peaks1 = peaks_from_ts_fs(binned_data(:,i),Fs);
    % 
    % plot(before{i})
    % hold on
    % plot([peaks1.PeakLocation],[peaks1.PeakValue],'bo')
    % 
    % z_breath_values(i) = calculate_z_blood(peaks1);
    continue
else
     bmObj=breathmetrics(binned_data(:,i),Fs,'humanAirflow');
 bmObj.estimateAllFeatures();

 % Define the variable names as field names in the struct
variableNames = {...
    'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
    'PercentOfBreathsWithInhalePause'};

values=[bmObj.secondaryFeatures.values];
% Assign the values to the struct fields
for ii = 1:length(variableNames)
    dataStruct.(variableNames{ii}) = values{ii};
end

 z_breath_values(i)= dataStruct(:);
end

end
end

