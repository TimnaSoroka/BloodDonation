
close all
clear
rng(50)
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


thresh=mean([subjData.Weight],'omitmissing');
sex_vec1=([subjData.Weight]<thresh);

% sex_vec1=logical([subjData.sex]);

sex_vec1([42,63,65])=[];
sex_vec=[sex_vec1,sex_vec1];

X_train=X(~sex_vec,:);
Y_train=Y(~sex_vec);

X_test=X(sex_vec,:);
Y_test=Y(sex_vec);

%% ----- Iterative downsampling to N=20 on the TRAIN set, AUC distribution -----
% Assumes you already created:
%   X_train, Y_train, X_test, Y_test, vars   (from your code above)
% Positive class in your labels is 2 (after)
posClass = 2;

% Settings
K_iter          = 10000;     % number of iterations
N_train_total   = 30;       % total N in training per iteration
min_per_class   = floor(N_train_total/2);   % target balance: ~10 vs 10
alpha_FS        = 0.05;     % univariate filter p-threshold
fallback_k      = 3;        % if no feature passes alpha, keep top-k by p
rng(123);                   % reproducibility for the iteration draws

% Classifier set (same as your list)
classifier_names = {'Logistic','kNN','DecisionTree','SVM-linear','SVM-rbf'}; %

% Storage
AUC_iter      = nan(K_iter,1);
best_name_it  = strings(K_iter,1);
nFeat_it      = nan(K_iter,1);

% Pre-compute class indices in the original TRAIN set
idx_c1_all = find(Y_train==1);
idx_c2_all = find(Y_train==2);

if numel(idx_c1_all) < min_per_class || numel(idx_c2_all) < min_per_class
    error('Not enough samples per class in X_train for N=12 balanced subsamples.');
end

for it = 1:K_iter
    % --------- (1) Subsample balanced N=20 from TRAIN ---------
    idx_c1 = randsample(idx_c1_all, min_per_class, false);
    idx_c2 = randsample(idx_c2_all, min_per_class, false);
    tr_idx = sort([idx_c1; idx_c2]);   % indices into X_train/Y_train

    Xtr_full = X_train(tr_idx,:);
    Ytr      = Y_train(tr_idx);

    % --------- (2) Univariate FS within this TRAIN subset ---------
    % Welch t-test feature-by-feature using labels (NO ordering assumptions)
    pvals = nan(1,size(Xtr_full,2));
    for f = 1:size(Xtr_full,2)
        x1 = Xtr_full(Ytr==1,f);
        x2 = Xtr_full(Ytr==2,f);
        % Handle all-NaN or constant edge cases
        if all(isnan(x1)) || all(isnan(x2)) || (nanstd(x1)==0 && nanstd(x2)==0)
            pvals(f) = 1;
        else
            [~,pvals(f)] = ttest2(x1, x2, 'Vartype','unequal'); % robust to unequal variances
        end
    end
    feat_mask = pvals < alpha_FS;

    % Fallback: if nothing passes, take top-k smallest p-values
    if ~any(feat_mask)
        [~,ord] = sort(pvals,'ascend','MissingPlacement','last');
        keep = ord(1:min(fallback_k, size(Xtr_full,2)));
        feat_mask = false(1,size(Xtr_full,2));
        feat_mask(keep) = true;
    end
    nFeat_it(it) = sum(feat_mask);

    Xtr = Xtr_full(:, feat_mask);
    Xte = X_test(:, feat_mask);

    % --------- (3) Median impute using TRAIN medians only ---------
    for f = 1:size(Xtr,2)
        medf = median(Xtr(:,f), 'omitnan');
        if isnan(medf), medf = 0; end
        % Train
        nan_tr = isnan(Xtr(:,f));
        if any(nan_tr), Xtr(nan_tr,f) = medf; end
        % Test
        nan_te = isnan(Xte(:,f));
        if any(nan_te), Xte(nan_te,f) = medf; end
    end

    % --------- (4) Pick best classifier by CV on the (small) TRAIN subset ---------
    Kcv = min(5, min(histcounts(Ytr, [0.5 1.5 2.5]))); % keep folds ≤ min class count
    if Kcv < 2, Kcv = 2; end
    cv = cvpartition(Ytr, 'KFold', Kcv);

    cv_acc = zeros(numel(classifier_names),1);
    for c = 1:numel(classifier_names)
        acc_fold = zeros(Kcv,1);
        for k = 1:Kcv
            tr = training(cv,k);  te = test(cv,k);
            Xtr_k = Xtr(tr,:);    Ytr_k = Ytr(tr);
            Xva_k = Xtr(te,:);    Yva_k = Ytr(te);

            switch classifier_names{c}
                case 'SVM-linear'
                    mdl = fitcsvm(Xtr_k, Ytr_k, 'KernelFunction','linear', 'Standardize',true);
                case 'SVM-rbf'
                    mdl = fitcsvm(Xtr_k, Ytr_k, 'KernelFunction','rbf',    'Standardize',true);
                case 'kNN'
                    mdl = fitcknn(Xtr_k, Ytr_k, 'NumNeighbors',5);
                case 'DecisionTree'
                    mdl = fitctree(Xtr_k, Ytr_k, 'SplitCriterion','gdi', 'MaxNumSplits',4, 'Surrogate','off');
                case 'Logistic'
                    mdl = fitclinear(Xtr_k, Ytr_k, 'Learner','logistic', 'Regularization','ridge');
            end
            Yhat = predict(mdl, Xva_k);
            acc_fold(k) = mean(Yhat == Yva_k);
        end
        cv_acc(c) = mean(acc_fold);
    end
    [~,best_idx]   = max(cv_acc);
    best_name      = classifier_names{best_idx};
    best_name_it(it) = best_name;

    % --------- (5) Fit best on full TRAIN subset and score TEST ---------
    switch best_name
        case 'SVM-linear'
            best_mdl = fitcsvm(Xtr, Ytr, 'KernelFunction','linear', 'Standardize',true);
        case 'SVM-rbf'
            best_mdl = fitcsvm(Xtr, Ytr, 'KernelFunction','rbf',    'Standardize',true);
        case 'kNN'
            best_mdl = fitcknn(Xtr, Ytr, 'NumNeighbors',5);
        case 'DecisionTree'
            best_mdl = fitctree(Xtr, Ytr, 'SplitCriterion','gdi', 'MaxNumSplits',4, 'Surrogate','off');
        case 'Logistic'
            best_mdl = fitclinear(Xtr, Ytr, 'Learner','logistic', 'Regularization','ridge');
    end

    [~, y_score] = predict(best_mdl, Xte);

    % Find the score column corresponding to posClass
    cls = best_mdl.ClassNames;
    idx_pos = find(ismember(cls, posClass));
    if isempty(idx_pos)
        error('Positive class not found in ClassNames. Check posClass.');
    end
    posScores = y_score(:, idx_pos);

    % --------- (6) AUC on fixed TEST set ---------
    [~,~,~,AUC] = perfcurve(Y_test, posScores, posClass);
    AUC_iter(it) = AUC;

% ---------- Significance of AUC on the fixed TEST set ----------
posClass = 2;  % your positive label

% Counts in TEST
n_pos = sum(Y_test==posClass);
n_neg = sum(Y_test~=posClass);

% (a) Hanley & McNeil normal-approx p-value for AUC > 0.5
Q1 = AUC / (2 - AUC);
Q2 = 2*AUC^2 / (1 + AUC);
varA = (AUC*(1-AUC) + (n_pos-1)*(Q1 - AUC^2) + (n_neg-1)*(Q2 - AUC^2)) / (n_pos*n_neg);
seA = sqrt(max(varA, eps));
Z  = (AUC - 0.5)/seA;
p_hanley = 1 - normcdf(Z);  % one-sided test H1: AUC > 0.5

% (b) Permutation p-value (robust to distributional assumptions)
B = 2000; 
auc_null = zeros(B,1);
labels = Y_test; scores = posScores;
for b = 1:B
    perm = labels(randperm(numel(labels)));
    [~,~,~,auc_null(b)] = perfcurve(perm, scores, posClass);
end
p_perm = mean(auc_null >= AUC);  % one-sided for AUC > 0.5

% store
AUC_iter(it)    = AUC;
p_hn_iter(it)   = p_hanley;
p_perm_iter(it) = p_perm;

end

% --------- Summary outputs ---------
fprintf('Iterations: %d | Train N per iter: %d (≈%d per class)\n', K_iter, N_train_total, min_per_class);
fprintf('AUC: mean=%.3f, SD=%.3f, median=%.3f, 2.5%%=%.3f, 97.5%%=%.3f\n', ...
    mean(AUC_iter,'omitnan'), std(AUC_iter,'omitnan'), median(AUC_iter,'omitnan'), ...
    prctile(AUC_iter,2.5), prctile(AUC_iter,97.5));

% Which classifier wins most often?
[uniq_names,~,ic] = unique(best_name_it);
freq = accumarray(ic,1);
[~,ord] = sort(freq,'descend');
disp(table(uniq_names(ord), freq(ord), 'VariableNames', {'Classifier','Wins'}));

% Optional: visualize distribution
figure('Color','w'); histogram(AUC_iter, 'BinWidth',0.02);
xlabel('AUC (Test set)'); ylabel('Count'); title('AUC across N=20 training iterations'); box off;

fprintf('AUC (N=20 training): mean=%.3f, sd=%.3f, median=%.3f, 2.5%%=%.3f, 97.5%%=%.3f\n', ...
    mean(AUC_iter,'omitnan'), std(AUC_iter,'omitnan'), median(AUC_iter,'omitnan'), ...
    prctile(AUC_iter,2.5), prctile(AUC_iter,97.5));

fprintf('Significant vs 0.5 (Hanley): %.1f%% of runs (alpha=0.05)\n', 100*mean(p_hn_iter<0.05));
fprintf('Significant vs 0.5 (Permutation): %.1f%% of runs (alpha=0.05)\n', 100*mean(p_perm_iter<0.05));
fprintf('Median p (perm) = %.4f\n', median(p_perm_iter,'omitnan'));

