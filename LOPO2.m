close all
clear
rng(1500)
Fs=25; 
load('Holter_timings.mat');

 exclude = {'019','028','032','042','095','105','027'}; %','025','105''027','078','069','105',
allCodes = {subjData.code};            % 1xN cell array of char/string
mask = ~ismember(allCodes, exclude);
subjData = subjData(mask);

%subjData([22,98])=[]; %less th
% an 10 minutes after in_chair 16,25,29,38,90
%subjData([16,22,27,36,88])=[];
% subjData([29,38,90])=[]; %less than 5 minutes after in chair 
%subjData([16,25,29,38,90])=[]; %less than 10 minutes after in_chair 16,25,29,38,90
%subjData(17)=[];
subjData([subjData.Weight]>90)=[];
%subjData([subjData.sex]==1)=[];

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


for i=1:size(subjData,2)
%[before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData); %_needle_walk_in_chair
%    [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle_lay(i,norm, IntLength,subjData); %_needle_walk_in_chair
[Before{i},After{i},donation{i},Before_LI{i},After_LI{i},NCdonation{i}]=extract_timings_needle_walk_in_chair2(i,norm, IntLength,subjData); %_needle_walk_in_chair

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
% x=x(:,[1,4,6,9]);
%  x=x(:,[1,3,4,6]);

fields=fieldnames(BmeasureResults);
% fields=fields([1,4,6,9]);
% fields=fields([1,3,4,6]);
fields=[BM_features_Names,fields'];
X=[X,x];
%fields=BM_features_Names;
% % v={subjData.P_Donation_Amount};
% %  v{60}='367';
% %  v{13}='558';
% %  v{31}='498';
% %  v{38}='500';
% % v_num = str2double(string(v));
% % 
% % w_all=[v_num,v_num]';
% % 
% % valid = ~isnan(w_all);
Y=[ones(size(vals_before,2),1);2*ones(size(vals_after,2),1)];
% % 
% % X = X(valid, :);     % keep only rows with valid weight
% % Y = Y(valid, :);
% % w = w_all(valid);    % cleaned weight
% % 
% % Z = [ones(size(w)) w];   % N × 2
% % 
% % % Fit linear model X ≈ Z * B  (all features at once)
% % B     = Z \ X;           % 2 × P   regression coefficients
% % X_hat = Z * B;           % N × P   part explained by weight
% % X_res = X - X_hat; 
% % X=X_res;
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


med_train = nanmedian(Xtr,1);
XtrImp = Xtr;  for j=1:size(XtrImp,2), m=isnan(XtrImp(:,j)); XtrImp(m,j)=med_train(j); end
XteImp = Xte;   for j=1:size(XteImp,2), m=isnan(XteImp(:,j)); XteImp(m,j)=med_train(j); end

% 2) Standardize using TRAIN stats
mu_tr = mean(XtrImp,1,'omitnan');
sd_tr = std(XtrImp,0,1,'omitnan'); sd_tr(sd_tr==0)=1;
XtrainZ = (XtrImp - mu_tr) ./ sd_tr;   % <-- this is your *only* z-scored train
XtestZ  = (XteImp - mu_tr) ./ sd_tr;   % same scaling applied to test

% 3) PCA fit on TRAIN only; transform both with same PCA
% [coeff, ~, ~, ~, expl, muPCA] = pca(XtrainZ);
% %k = find(cumsum(expl)>=90,1,'first');
% k=28;
% Xtr_sel = (XtrainZ - muPCA) * coeff(:,1:k);
% Xte_sel  = (XtestZ  - muPCA) * coeff(:,1:k);

Xtr_sel=XtrainZ;
Xte_sel=XtestZ;

    %FS
   %  [h,pt]=ttest2(Xtr(1:size(Xtr,1)/2,:),Xtr(size(Xtr,1)/2+1:end,:));
   % 
   %  [~,pidx]=sort(pt);
   % 
   % Xtr_sel=Xtr(:,pidx(1:3));
   % Xte_sel=Xte(:,pidx(1:3));
% Xtr_sel(:,[17,21])=[];
% Xte_sel(:,[17,21])=[];

               mdl = fitclinear(Xtr_sel,Ytr,'Learner','logistic');
             %                   mdl = fitcsvm(Xtr_sel,Ytr,'KernelFunction','linear');

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

v={subjData.P_Donation_Amount};
 v{60}='367';
 v{13}='558';
 v{31}='498';
 v{38}='500';
v_num = str2double(string(v));
sum(~isnan(v_num));
% Keep only valid numeric entries

w=[subjData.Weight];

bloodVolume  = 70 .* w;
fractionLost = v_num ./ w;


% Correlation between weight and model score
for iii=1:size(X,2)
[r1, p1] = corr(fractionLost', X(1:N,iii), 'type', 'Spearman','rows','complete');
[r2, p2] = corr(fractionLost', X(N+1:end,iii), 'type', 'Spearman','rows','complete');
if p1<0.05
fprintf(['r1: ' num2str(r1) ' field:' fields{iii} ' \n'])
elseif p2<0.05
fprintf(['r2: ' num2str(r2) ' field:' fields{iii} ' \n'])
end
end

[r1, p1] = corr(fractionLost', score_all(1:N), 'type', 'Spearman','rows','complete');
[r2, p2] = corr(fractionLost', score_all(N+1:end), 'type', 'Spearman','rows','complete');
fprintf(['r1: ' num2str(r1) ' \n'])
fprintf(['r2: ' num2str(r2) ' \n'])

plot(w', score_to_comp,'ko')
lsline()