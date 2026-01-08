close all
clear
rng(1500)
Fs=25; 
load('Holter_timings.mat');
corrThr = 0.8;     %
method  = 'spearman'; % or 'spearman' if you prefer rank-based robustness


 exclude = {'032','042'}; %'025','105','019','028','032','042','095','025','105'

% rmCodes  = {'001','002','003','013','014','015','016','017','018','019','020','021','022','023','024','025',...
%     '026','027','028','029','030','031','032','042','096','097','098','099','100','101','102','103','104','105','106','107'};
% allCodes = {subjData.code};
% 
% %001-003	בא"ח חיפה	02/02/2025	ניידת
% %013-031	בא"ח חיפה	11-13/02/2025	ניידת
% %096-102	קרייה	16/07/2025	ניידת
% %103-107	תלה"ש	14/08/2025	ניידת
 subjData(ismember({subjData.code}, rmCodes)) = [];

allCodes = {subjData.code};            % 1xN cell array of char/string
mask = ~ismember(allCodes, exclude);
subjData = subjData(mask);

%subjData([22,98])=[]; %less th
% an 10 minutes after in_chair 16,25,29,38,90
%subjData([16,22,27,36,88])=[];
% subjData([29,38,90])=[]; %less than 5 minutes after in chair 
%subjData([16,25,29,38,90])=[]; %less than 10 minutes after in_chair 16,25,29,38,90
%subjData(17)=[];
%subjData([subjData.Weight]>90)=[];

BM_features_Names = {...
    'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
    'PercentOfBreathsWithInhalePause'};

N=size(subjData,2);
IntLength=5;
norm=1;

for i = 1:numel(subjData)
    [Before{i}, After{i},donation{i},Before_LI{i}, After_LI{i},NCdonation{i}]=extract_timings_needle_walk_in_chair2(i,norm, IntLength,subjData);

end

vals_before=calculate_before_after(Before);
vals_after=calculate_before_after(After);

 for i=1:size(Before_LI,2)
     % if ismember(i,[42,63,65]) 
     %     Fs=6;
     % else
     %     Fs=25;
     % end
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters(Before_LI{i},Fs,0,'Holter');
[Laterality_IndexA(i),AmeasureResults(i)]=NasalCycleParameters(After_LI{i},Fs,0,'Holter');
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
 end

X=table2array(struct2table([vals_before,vals_after]));

%%
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
 x=x(:,[1,9,4,6]);
 %x=x(:,[2,5,7,10]);

fields=fieldnames(BmeasureResults);
 fields=fields([1,4,6,9]);
% fields=fields([1,3,4,6]);
fields=[BM_features_Names,fields'];
vals_before_c=calculate_before_after_catch(Before);
[vals_after_c,vars]=calculate_before_after_catch(After);
 X22=table2array(struct2table([vals_before_c,vals_after_c]));

demo=[subjData.Weight]';
demo=[demo;demo];

demo2=[subjData.sex]';
demo2=[demo2;demo2];

X=[X,X22,x,demo2,demo]; %X22,X=[X22,X,x];,,demo2

Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];

%% === LOPO with L1-logistic feature selection (train-only) ===
% Assumes X (rows = [before; after] per participant) and Y (1=Before, 2=After) already exist.
% Uses lassoglm (binomial) with inner CV to pick lambda; nonzero weights = selected features.

posClass = 2;   % AFTER
assert(sum(Y==1)==sum(Y==2), 'Expect one BEFORE and one AFTER per participant.');
nSubj = sum(Y==1);
pid   = [(1:nSubj)'; (1:nSubj)'];   % participant id per row
nFeat = size(X,2);

y_true_all = []; 
y_pred_all = []; 
score_all  = [];          % P(class==2) for ROC
fold_acc   = nan(nSubj,1);

for s = 1:nSubj
    teMask = (pid == s);
    trMask = ~teMask;

    Xtr = X(trMask,:);  Ytr = Y(trMask);
    Xte = X(teMask,:);  Yte = Y(teMask);
 C = corr(Xtr, 'Type', method, 'Rows', 'pairwise');  % 52x52

    % Greedy removal: for each highly correlated pair, drop one.
    % Heuristic: drop the one with higher mean absolute correlation to others.
    meanAbsCorr = mean(abs(C - diag(diag(C))), 2); % exclude diagonal

    toDrop = false(size(meanAbsCorr));
    [rIdx, cIdx] = find(triu(abs(C),1) >= corrThr);

    % process pairs in descending correlation
    [~, ord] = sort(abs(C(sub2ind(size(C), rIdx, cIdx))), 'descend');
    rIdx = rIdx(ord); cIdx = cIdx(ord);

    for k = 1:numel(rIdx)
        a = rIdx(k); b = cIdx(k);
        if toDrop(a) || toDrop(b), continue; end

        % drop the "more redundant" one (higher avg abs corr)
        if meanAbsCorr(a) >= meanAbsCorr(b)
            toDrop(a) = true;
        else
            toDrop(b) = true;
        end
    end

    keepIdx = ~toDrop;

    Xtr2 = Xtr(:, keepIdx);
    Xte2 = Xte(:, keepIdx);


med_train = nanmedian(Xtr,1);
XtrImp = Xtr2;  for j=1:size(XtrImp,2), m=isnan(XtrImp(:,j)); XtrImp(m,j)=med_train(j); end
XteImp = Xte2;   for j=1:size(XteImp,2), m=isnan(XteImp(:,j)); XteImp(m,j)=med_train(j); end

% 2) Standardize using TRAIN stats
mu_tr = mean(XtrImp,1,'omitnan');
sd_tr = std(XtrImp,0,1,'omitnan'); sd_tr(sd_tr==0)=1;
XtrainZ = (XtrImp - mu_tr) ./ sd_tr;   % <-- this is your *only* z-scored train
XtestZ  = (XteImp - mu_tr) ./ sd_tr;   % same scaling applied to test

% 3) PCA fit on TRAIN only; transform both with same PCA
[coeff, ~, ~, ~, expl, muPCA] = pca(XtrainZ);
k = find(cumsum(expl)>=80,1,'first');
 % k=round(sqrt(size(XtrainZ,1)));
Xtr_sel = (XtrainZ - muPCA) * coeff(:,1:k);
Xte_sel  = (XtestZ  - muPCA) * coeff(:,1:k);
% 
% Xtr_sel=XtrainZ;
% Xte_sel=XtestZ;
% mrmrL=fscmrmr(XtrainZ,Ytr);
% Xtr_sel=XtrainZ(:,mrmrL(1:37));
% Xte_sel=XtestZ(:,mrmrL(1:37));

   % Xtr_sel=Xtr(:,pidx(1:3));
   % Xte_sel=Xte(:,pidx(1:3));
% Xtr_sel(:,[17,21])=[];
% Xte_sel(:,[17,21])=[];

              mdl = fitclinear(Xtr_sel,Ytr,'Learner','logistic');
                  
       %              mdl = fitcsvm(Xtr_sel,Ytr,'KernelFunction','linear');
    % 
    %               mdl = fitctree(...
    % Xtr_sel, ...
    % Ytr, ...
    % 'SplitCriterion', 'gdi', ...
    % 'MaxNumSplits', 100, ...
    % 'Surrogate', 'off');

%     subspaceDimension = max(1, min(6, width(Xtr_sel) - 1));
% mdl = fitcensemble(...
%     Xtr_sel, ...
%     Ytr, ...
%     'Method', 'Subspace', ...
%     'NumLearningCycles', 30, ...
%     'Learners', 'discriminant', ...
%     'NPredToSample', subspaceDimension);


[y_pred,y_score] = predict(mdl, Xte_sel);

        fold_acc(s) = sum(y_pred == Yte)/length(Yte);

    y_true_all = [y_true_all; Yte];
    y_pred_all = [y_pred_all; y_pred];
    score_all  = [score_all;  y_score(:,2)];

end

% ---- Aggregate metrics ----
overall_acc = mean(y_true_all == y_pred_all);
cm = confusionmat(y_true_all, y_pred_all, 'Order', [1 2]);
[fpRate, tpRate, ~, AUC] = perfcurve(y_true_all, score_all, posClass);

fprintf('LOPO (L1-logistic): Accuracy = %.4f, AUC = %.4f\n', overall_acc, AUC);
disp('Confusion matrix (rows=true [1,2], cols=pred [1,2]):');
disp(cm);

figure; 
plot(fpRate, tpRate, 'LineWidth', 2); hold on; plot([0 1],[0 1],'k--');
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title(sprintf('LOPO ROC — AUC = %.2f', AUC)); grid on; axis square;


% v={subjData.P_Donation_Amount};
%  v{60}='367';
%  v{13}='558';
%  v{31}='498';
%  v{38}='500';
% v_num = str2double(string(v));
% sum(~isnan(v_num));
% % Keep only valid numeric entries
% 
% w=[subjData.Weight];
% 
% bloodVolume  = 70 .* w;
% fractionLost = v_num ./ w;
% 
% 
% % Correlation between weight and model score
% for iii=1:size(X,2)
% [r1, p1] = corr(fractionLost', X(1:N,iii), 'type', 'Spearman','rows','complete');
% [r2, p2] = corr(fractionLost', X(N+1:end,iii), 'type', 'Spearman','rows','complete');
% if p1<0.05
% fprintf(['r1: ' num2str(r1) ' field:' fields{iii} ' \n'])
% elseif p2<0.05
% fprintf(['r2: ' num2str(r2) ' field:' fields{iii} ' \n'])
% end
% end
% 
% [r1, p1] = corr(fractionLost', score_all(1:N), 'type', 'Spearman','rows','complete');
% [r2, p2] = corr(fractionLost', score_all(N+1:end), 'type', 'Spearman','rows','complete');
% fprintf(['r1: ' num2str(r1) ' \n'])
% fprintf(['r2: ' num2str(r2) ' \n'])
% 
% plot(w', score_to_comp,'ko')
% lsline()