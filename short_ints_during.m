close all
clear
rng(1500)
Fs=25; 
load(['Holter_timings_controls.mat']);
to_plot=1;

%subjData([29,38,90])=[]; %less than 5 minutes after in chair 
%subjData([16,25,29,38,90])=[]; %less than 10 minutes after in_chair 16,25,29,38,90

BM_features_Names = {...
    'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
    'PercentOfBreathsWithInhalePause'};

N=size(subjData,2);

for i = 1:numel(subjData)
    % 5 minute interval case
    % if ismember(i,[29,38,90])
    %     continue
    % end
         % 10 minute interval case


    T = subjData(i).Data;                        % table
    resp_stereo = table2array(T(:,[3 4]));      % make sure this is defined for both branches
    Fs_raw = 25;                                % default sampling rate
    rec_minutes = size(T,1) / (60*Fs_raw);

    % Handle the 6 Hz cases by resampling to 25 Hz AND scaling in/out indices
    if rec_minutes < 20
        warning('check_sampling_rate for participant %s', subjData(i).code);
        Fs_raw = 6;
        X = sum(resp_stereo,2);
        X = zscore(resample(X, 25, 6));
        resp_stereo=resample(resp_stereo, 25, 6);
        Fs = 25;
        scale = Fs / Fs_raw;
        in_idx  = max(1, min(numel(X), round(subjData(i).in  * scale)));
        out_idx = max(1, min(numel(X), round(subjData(i).out * scale)));
                in_chair = max(1, min(numel(X), round(subjData(i).in_chair * scale)));
        walk = max(1, min(numel(X), round(subjData(i).walk * scale)));

    else
        X = zscore(sum(resp_stereo,2));
        in_idx  = subjData(i).in;
        out_idx = subjData(i).out;
       in_chair = subjData(i).in_chair;
       walk = subjData(i).walk;

    end

    During{i}=X(in_idx+30*25:out_idx-30*25);
During_in_minute(i)=numel(During{i})/(60*25);
Before{i}=X(in_idx+30*25:2.5*60*25+in_idx);
          Before_LI{i}=resp_stereo(in_idx+30*25:2.5*60*25+in_idx,:);

          After{i}=X(out_idx-2.5*60*25:out_idx-30*25);
                    After_LI{i}=resp_stereo(out_idx-2.5*60*25:out_idx-30*25,:);


Before{i}=X(in_idx+30*25:2.5*60*25+in_idx);
          Before_LI{i}=resp_stereo(in_idx+30*25:2.5*60*25+in_idx,:);

          After{i}=X(out_idx-2.5*60*25:out_idx-30*25);
                    After_LI{i}=resp_stereo(out_idx-2.5*60*25:out_idx-30*25,:);

                    After_in_minute(i)=numel(After{i})/(60*25);
Before_in_minute(i)=numel(Before{i})/(60*25);
end

histogram(During_in_minute,'NumBins',20)
After_in_minute(After_in_minute==0)=nan;
% mean(After_in_minute,'omitnan')
% sum(~isnan(After_in_minute))

After=After(~isnan(After_in_minute));

Before_in_minute(Before_in_minute==0)=nan;
% mean(Before_in_minute,'omitnan')
% sum(~isnan(Before_in_minute))
Before=Before(~isnan(Before_in_minute));


vals_before_controls=calculate_before_after(Before);
vals_after_controls=calculate_before_after(After);

for i=1:size(Before_LI,2)
     % if ismember(i,[42,63,65]) 
     %     Fs=6;
     % else
     %     Fs=25;
     % end
[Laterality_IndexB(i),BmeasureResults_con(i)]=NasalCycleParameters_short(Before_LI{i},Fs,0,'Holter');
[Laterality_IndexA(i),AmeasureResults_con(i)]=NasalCycleParameters_short(After_LI{i},Fs,0,'Holter');
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
end

load(['Holter_timings.mat']);

N=size(subjData,2);

for i = 1:numel(subjData)
    % 5 minute interval case
    % if ismember(i,[29,38,90])
    %     continue
    % end
         % 10 minute interval case


    T = subjData(i).Data;                        % table
    resp_stereo = table2array(T(:,[3 4]));      % make sure this is defined for both branches
    Fs_raw = 25;                                % default sampling rate
    rec_minutes = size(T,1) / (60*Fs_raw);

    % Handle the 6 Hz cases by resampling to 25 Hz AND scaling in/out indices
    if rec_minutes < 20
        warning('check_sampling_rate for participant %s', subjData(i).code);
        Fs_raw = 6;
        X = sum(resp_stereo,2);
        X = zscore(resample(X, 25, 6));
        resp_stereo=resample(resp_stereo, 25, 6);
        Fs = 25;
        scale = Fs / Fs_raw;
        in_idx  = max(1, min(numel(X), round(subjData(i).in  * scale)));
        out_idx = max(1, min(numel(X), round(subjData(i).out * scale)));
                in_chair = max(1, min(numel(X), round(subjData(i).in_chair * scale)));
        walk = max(1, min(numel(X), round(subjData(i).walk * scale)));

    else
        X = zscore(sum(resp_stereo,2));
        in_idx  = subjData(i).in;
        out_idx = subjData(i).out;
       in_chair = subjData(i).in_chair;
       walk = subjData(i).walk;

    end

    During{i}=X(in_idx+30*25:out_idx-30*25);
During_in_minute(i)=numel(During{i})/(60*25);
Before{i}=X(in_idx+30*25:2.5*60*25+in_idx);
          Before_LI{i}=resp_stereo(in_idx+30*25:2.5*60*25+in_idx,:);

          After{i}=X(out_idx-2.5*60*25:out_idx-30*25);
                    After_LI{i}=resp_stereo(out_idx-2.5*60*25:out_idx-30*25,:);

                    After_in_minute(i)=numel(After{i})/(60*25);
Before_in_minute(i)=numel(Before{i})/(60*25);
end

histogram(During_in_minute,'NumBins',20)
After_in_minute(After_in_minute==0)=nan;
% mean(After_in_minute,'omitnan')
% sum(~isnan(After_in_minute))

After=After(~isnan(After_in_minute));

Before_in_minute(Before_in_minute==0)=nan;
% mean(Before_in_minute,'omitnan')
% sum(~isnan(Before_in_minute))
Before=Before(~isnan(Before_in_minute));


vals_before=calculate_before_after(Before);
vals_after=calculate_before_after(After);


 for i=1:size(Before_LI,2)
     % if ismember(i,[42,63,65]) 
     %     Fs=6;
     % else
     %     Fs=25;
     % end
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters_short(Before_LI{i},Fs,0,'Holter');
[Laterality_IndexA(i),AmeasureResults(i)]=NasalCycleParameters_short(After_LI{i},Fs,0,'Holter');
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
 end

X=table2array(struct2table([vals_before,vals_after]));
X_controls=table2array(struct2table([vals_before_controls,vals_after_controls]));

%%
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
 x=x(:,[1,3,4,6]);
 x_con=table2array(struct2table([BmeasureResults_con,AmeasureResults_con]));
 x_con=x_con(:,[1,3,4,6]);
fields=fieldnames(BmeasureResults);
% fields=fields([1,4,6,9]);
 fields=fields([1,3,4,6]);
fields=[BM_features_Names,fields'];
X_train_enc=[X,x];
X_con=[X_controls,x_con];

Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];
% % 

for i=1:size(fields,2)
    currentfield=fields{i};
test_values = [X_train_enc(1:size(vals_before,2),i)];   % Test scores 
retest_values =[X_train_enc(size(vals_before,2)+1:end,i)]; % Retest scores 

delta=test_values-retest_values;
outlierMask = isoutlier(delta,"mean");      % Detect outliers (default: median + MAD)
retest_values = retest_values(~outlierMask);   % Remove outliers
test_values = test_values(~outlierMask);   % Remove outliers

%[p_value24(i),~,Wstat24(i)] = signrank(test_values, retest_values,'method','approximate');
[~,tp_value24(i),~,tstat24(i)] = ttest(test_values, retest_values);


                test_values_con = [X_con(1:size(vals_before_controls,2),i)];   % Test scores 
                retest_values_con =[X_con(size(vals_before_controls,2)+1:end,i)]; % Retest scores 
                
                delta=test_values_con-retest_values_con;
                outlierMask = isoutlier(delta,"mean");      % Detect outliers (default: median + MAD)
                retest_values_con = retest_values_con(~outlierMask);   % Remove outliers
                test_values_con = test_values_con(~outlierMask);   % Remove outliers
                
                %[p_value24(i),~,Wstat24(i)] = signrank(test_values, retest_values,'method','approximate');
                [~,tp_value24_con(i),~,tstat24_con(i)] = ttest(test_values_con, retest_values_con);


%  [~,tp_value24_women(i)] = ttest(test_values(sex_vec_clean2), retest_values(sex_vec_clean2));
%   [~,tp_value24_men(i)] = ttest(test_values(~sex_vec_clean2), retest_values(~sex_vec_clean2));

if tp_value24(i)<1
if to_plot
    figure;
min_=min([test_values;retest_values;test_values_con;retest_values_con]);
max_=max([test_values;retest_values;test_values_con;retest_values_con]);
scatter(test_values, retest_values, 100, 'filled','MarkerFaceColor',[0 0.447 0.741]); % Blue filled markers
 hold on
scatter(test_values_con, retest_values_con, 100, 'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]); % control

plot([min_ max_],[min_ max_], 'k--'); % Identity line
xlim([min_ max_])
ylim([min_ max_])
% Perform Mann-Whitney U test (non-parametric test for paired samples)
set(gca,'FontSize',12)

% Format title with p-value
title(sprintf('%s (ttest p = %.4f)',fields{i},tp_value24(i)));

% Axis labels and limits
xlabel('Before');
ylabel('After');
        end
%                 % Descriptives
%     mx = mean(test_values, 'omitnan'); sx = std(test_values, 'omitnan'); nx = numel(test_values);
%     my = mean(retest_values, 'omitnan'); sy = std(retest_values, 'omitnan'); ny = numel(retest_values);
% 
%     valid = ~isnan(test_values) & ~isnan(retest_values);
%         x = test_values(valid);
%         y = retest_values(valid);
%         [~,p,~,stats] = ttest(x, y);
%         % Effect size: Cohen's dz for paired (mean diff / SD diff)
%             df = stats.df; tval = stats.tstat;
% varName=currentfield;
%     fprintf('%s: t(%d)=%.2f, p=%.4g]\n', ...
%         varName, df, tval, p);
%     fprintf('   before: %0.3f \xB1 %0.3f (n=%d);  after: %0.3f \xB1 %0.3f (n=%d)\n', ...
%         mx, sx, nx, my, sy, ny);
%  
%  fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])
disp_descriptives(test_values,retest_values,currentfield)
disp_descriptives(test_values_con,retest_values_con,currentfield)
end
end
%%
[~,pidx]=sort(tp_value24);
%topF=pidx(1:round(sqrt(N)));
topF=pidx(1:5);
X_train_sel=X_train_enc(:,topF);

% % X_train_sel : n×p features
% % Y            : n×1 labels (e.g., [1 2] or [0 1])
% 
% cvp    = cvpartition(Y,'Leaveout');
% scores = nan(size(Y));   % prob. for positive class
% yhat   = nan(size(Y));   % predicted label
% 
% for k = 1:cvp.NumTestSets
%     tr = training(cvp,k); 
%     te = test(cvp,k);
% 
%     % ----- (optional but recommended) standardize using TRAIN stats -----
%     mu = mean(X_train_sel(tr,:),1); 
%     sd = std(X_train_sel(tr,:),0,1); sd(sd==0) = 1;
%     Xtr = (X_train_sel(tr,:) );%- mu) ./ sd;
%     Xte = (X_train_sel(te,:)) ;%- mu) ./ sd;
% 
%     ytr = Y(tr);
% 
%     % Guard: need both classes present in the training fold
%     if numel(unique(ytr)) < 2
%         % fallback: predict the majority class; prob = 0 or 1 accordingly
%         maj = mode(ytr);
%         yhat(te)   = maj;
%         scores(te) = double(maj == max(Y)); % crude prob
%         continue
%     end
% 
%     % ----- logistic regression (linear, with a bit of ridge for stability) -----
%     M = fitclinear(Xtr, ytr, ...
%         'Learner','logistic', ...
%         'Regularization','ridge', ...   % helps with collinearity
%         'Lambda',1e-3, ...              % you can also try 'Lambda','auto'
%         'ClassNames',unique(Y));        % keeps class order consistent
% 
%     [lab,sc] = predict(M, Xte);         % sc(:,2) = P(positive class)
%     yhat(te)   = lab;
%     % Identify which column is the positive class you care about:
%     posClass   = max(M.ClassNames);     % e.g., 2 if classes are [1 2]
%     posCol     = find(M.ClassNames==posClass,1,'first');
%     scores(te) = sc(:,posCol);
% end
% 
% % --- Metrics ---
% acc = mean(yhat == Y);
% [posX,posY,~,AUC] = perfcurve(Y, scores, max(Y));  % positive = max label (e.g., 2)
% 
% fprintf('LOOCV logistic: Accuracy = %.3f, AUC = %.3f\n', acc, AUC);

classificationSVM = fitcsvm(...
    X_train_sel, ...
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

%Compute ROC
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
validIdx = all(~isnan(X_train_sel), 2);

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


function disp_descriptives(test_values,retest_values,currentfield)
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
end