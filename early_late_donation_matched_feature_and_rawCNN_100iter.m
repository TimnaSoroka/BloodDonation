%% early_late_donation_matched_donor_control_100iter.m
%
% Adapted from early_late_donation.m.
%
% Goal:
%   Build early-vs-late donation windows from raw subjData for DONORS and CONTROLS,
%   extract the same feature set, then run a matched validation:
%
%   Default donor-trained test:
%       for each repeat:
%           hold out N_controls donors
%           train on remaining donors
%           test the SAME model on:
%               1) held-out donors
%               2) all controls
%
%   This directly tests whether a model trained on donor early/late donation
%   dynamics generalizes to controls, using a matched donor test-set size.
%
% Main outputs:
%   RawWindowInfo.csv
%   FeatureTable_fromRawWindows.csv
%   MatchedFeature_Results_100.csv
%   MatchedFeature_Predictions_100.csv
%   MatchedFeature_ParticipantScores_100.csv
%   MatchedFeature_Summary_100.csv

clear; clc;
rng(42);

%% ======================== USER SETTINGS ========================

CFG = struct();

% -------- Input files --------
CFG.donorSubjDataFile   = '/Users/timnas/Documents/BloodDonation/Holter_timings.mat';
CFG.controlSubjDataFile = '/Users/timnas/Documents/BloodDonation/Holter_timings_controls.mat';

% -------- Output folder --------
CFG.baseOutDir = '/Users/timnas/Documents/BloodDonation/rawSubjData_earlyLateDonation_matched_donor_control';

% Raw channel variable names inside subjData(s).Data.
CFG.leftChannelCandidates  = {'L','left','Left','leftResp','respLeft','L_resp','channelL'};
CFG.rightChannelCandidates = {'R','right','Right','rightResp','respRight','R_resp','channelR'};

%% ================= CHANNEL / NORMALIZATION / FEATURE-SELECTION OPTIONS =================

% Channel mode:
% 'stereo' = use L and R separately + nostril coordination features
% 'sum'    = use only one summed/mean channel, no nostril coordination features
CFG.channelMode = 'sum';      % 'stereo' or 'sum'

% For sum mode:
CFG.sumChannelMethod = 'sum';   % 'mean' or 'sum'

% Raw normalization stage:
% 'none'      = no raw normalization
% 'beforeCut' = normalize full raw L/R recording before cutting donation interval
% 'afterCut'  = normalize donation segment after cutting and preprocessing
% 'perWindow' = normalize each model window separately
CFG.rawNormalizeStage = 'afterCut';
CFG.rawNormalizeMethod = 'robustZ';   % 'robustZ' or 'zscore'

% Nested feature selection. If true, selection is done inside each repeat
% using the TRAINING participants only.
CFG.useNestedFeatureSelection = false;
CFG.fsK = 20;
CFG.fsCorrPrune = true;
CFG.fsCorrThreshold = 0.90;

% Sampling rate.
CFG.defaultFs = 25;
CFG.targetFs  = 25;

% Known 6 Hz participants. These are participant CODES, not serial indices.
CFG.forceFs6DonorParticipants   = ["045","067","069"];
CFG.forceFs6ControlParticipants = ["007","020"];
CFG.excludeFs6FromFeaturePrimary = true;

% Preprocessing.
CFG.doBandpass = true;
CFG.bandpassHz = [0.05 1.20];
CFG.doDetrend = true;

% Donation edge exclusion.
CFG.excludeFirstSec = 5;
CFG.excludeLastSec  = 5;

% Windowing.
CFG.windowSec = 60;
CFG.overlap = 0.50;

% Early/late definitions using normalized donation time.
CFG.earlyRange = [0.05 0.30];
CFG.lateRange  = [0.70 0.95];

% Feature model.
CFG.featureCenteringMode = 'none';
% Options:
%   'none'
%   'participantDonationMean'
%   'participantDonationZ'
%   'participantEarlyMean'

CFG.featureModelType = 'ridgeLogistic';
% Options:
%   'ridgeLogistic'
%   'linearSVM'

CFG.lambda = 1e-3;

%% ================= MATCHED EVALUATION OPTIONS =================

CFG.nMatchedRepeats = 100;
CFG.randomSeed = 42;

% Default: donor-trained model, test held-out donors vs all controls.
% This implements: train on donors minus N_controls, test on controls or matched donors.
CFG.trainGroup = "donor";      % "donor" or "control"

% For donor-trained mode:
%   test N_controls donors and all controls each repeat.
% For control-trained mode:
%   because controls are few, hold out a fraction of controls and test the
%   same number of randomly selected donors.
CFG.controlTrain_testControlFrac = 0.30;
CFG.controlTrain_minTestControls = 3;

CFG.treatParticipantAsGroupUnit = true;
CFG.makePlots = true;

% -------------------------
% Raw CNN matched evaluation
% -------------------------
% The earlier matched script evaluated only extracted-feature classifiers.
% Set this true to run the same participant-level donor/control test directly
% on raw early/late windows using the CNN branch from early_late_donation.m.
CFG.runMatchedRawCNNPipeline = true;

% Optional: keep the feature-classifier branch too, for side-by-side comparison.
CFG.runMatchedFeaturePipeline = true;

% Participant-level train/test split for the raw CNN.
% trainGroup participants are split into trainFrac / held-out test each repeat.
CFG.rawCNNTrainFrac = 0.80;

% Opposite-group test set:
%   "all"     = test all available opposite-group participants each repeat.
%   "matched" = test a random opposite-group subset matched to the held-out same-group N.
CFG.rawCNNOppositeTestMode = "all";

CFG.excludeFs6FromCNN = true;
CFG.maxEpochs = 25;
CFG.miniBatchSize = 32;
CFG.initialLearnRate = 1e-3;


CFG.analysisTag = sprintf('%s_norm-%s_FS-%d_K%d_%sTrain_%drep', ...
    CFG.channelMode, ...
    CFG.rawNormalizeStage, ...
    CFG.useNestedFeatureSelection, ...
    CFG.fsK, ...
    char(CFG.trainGroup), ...
    CFG.nMatchedRepeats);

CFG.outDir = fullfile(CFG.baseOutDir, CFG.analysisTag);
if ~exist(CFG.outDir, 'dir')
    mkdir(CFG.outDir);
end

%% ======================== LOAD subjData ========================

fprintf('\nLoading DONOR subjData from:\n%s\n', CFG.donorSubjDataFile);
subjDataDonor = loadSubjDataStruct(CFG.donorSubjDataFile);
fprintf('Loaded donor subjData with %d participants.\n', numel(subjDataDonor));

fprintf('\nLoading CONTROL subjData from:\n%s\n', CFG.controlSubjDataFile);
subjDataControl = loadSubjDataStruct(CFG.controlSubjDataFile);
fprintf('Loaded control subjData with %d participants.\n', numel(subjDataControl));

%% ======================== BUILD RAW WINDOWS ========================

fprintf('\nBuilding raw early/late windows for donors...\n');
CFGdonor = CFG;
CFGdonor.forceFs6Participants = CFG.forceFs6DonorParticipants;
[RawX_D, Y_D, WindowInfo_D] = buildRawWindowsFromSubjData(subjDataDonor, CFGdonor);
WindowInfo_D.group = repmat("donor", height(WindowInfo_D), 1);
WindowInfo_D.participantKey = WindowInfo_D.group + "_" + string(WindowInfo_D.participant);
WindowInfo_D = movevars(WindowInfo_D, {'group','participantKey'}, 'After', 'participant');

fprintf('\nBuilding raw early/late windows for controls...\n');
CFGcontrol = CFG;
CFGcontrol.forceFs6Participants = CFG.forceFs6ControlParticipants;
[RawX_C, Y_C, WindowInfo_C] = buildRawWindowsFromSubjData(subjDataControl, CFGcontrol);
WindowInfo_C.group = repmat("control", height(WindowInfo_C), 1);
WindowInfo_C.participantKey = WindowInfo_C.group + "_" + string(WindowInfo_C.participant);
WindowInfo_C = movevars(WindowInfo_C, {'group','participantKey'}, 'After', 'participant');

RawX = [RawX_D; RawX_C];
Y = [Y_D; Y_C];
WindowInfo = [WindowInfo_D; WindowInfo_C];

fprintf('\nCreated %d total windows.\n', numel(RawX));
fprintf('Early windows: %d\n', sum(Y == "early"));
fprintf('Late windows:  %d\n', sum(Y == "late"));
fprintf('Donor participants:   %d\n', numel(unique(WindowInfo.participantKey(WindowInfo.group=="donor"))));
fprintf('Control participants: %d\n', numel(unique(WindowInfo.participantKey(WindowInfo.group=="control"))));

writetable(WindowInfo, fullfile(CFG.outDir, 'RawWindowInfo.csv'));

%% ======================== FEATURE TABLE ========================

fprintf('\n================ FEATURE EXTRACTION ================\n');

if isfield(CFG, 'excludeFs6FromFeaturePrimary') && CFG.excludeFs6FromFeaturePrimary
    keepFeature = ~WindowInfo.isForced6Hz;

    RawX_Feature = RawX(keepFeature);
    Y_Feature = Y(keepFeature);
    WindowInfo_Feature = WindowInfo(keepFeature,:);

    fprintf('FEATURE: excluding forced-6Hz participants.\n');
    fprintf('FEATURE: windows kept: %d/%d\n', numel(RawX_Feature), numel(RawX));
    fprintf('FEATURE: participants kept: %d/%d\n', ...
        numel(unique(WindowInfo_Feature.participantKey)), ...
        numel(unique(WindowInfo.participantKey)));

    fprintf('FEATURE: excluded participant keys:\n');
    disp(unique(WindowInfo.participantKey(WindowInfo.isForced6Hz)));
else
    RawX_Feature = RawX;
    Y_Feature = Y;
    WindowInfo_Feature = WindowInfo;
end

FeatureTable = computeFeatureTableFromRawWindows_Matched(RawX_Feature, Y_Feature, WindowInfo_Feature, CFG);
writetable(FeatureTable, fullfile(CFG.outDir, 'FeatureTable_fromRawWindows.csv'));

fprintf('Feature table saved. Rows = %d | Columns = %d\n', ...
    height(FeatureTable), width(FeatureTable));

%% ======================== MATCHED TRAIN / TEST: FEATURES ========================

if CFG.runMatchedFeaturePipeline

    fprintf('\n================ MATCHED FEATURE EVALUATION ================\n');
    matchedRes = runMatchedFeatureGroupEvaluation(FeatureTable, CFG);

    writetable(matchedRes.results, ...
        fullfile(CFG.outDir, sprintf('MatchedFeature_Results_%d.csv', CFG.nMatchedRepeats)));

    writetable(matchedRes.predictions, ...
        fullfile(CFG.outDir, sprintf('MatchedFeature_Predictions_%d.csv', CFG.nMatchedRepeats)));

    writetable(matchedRes.participantScores, ...
        fullfile(CFG.outDir, sprintf('MatchedFeature_ParticipantScores_%d.csv', CFG.nMatchedRepeats)));

    writetable(matchedRes.summary, ...
        fullfile(CFG.outDir, sprintf('MatchedFeature_Summary_%d.csv', CFG.nMatchedRepeats)));

    disp(matchedRes.summary);

    if CFG.makePlots
        plotMatchedAucDistributions(matchedRes.results, CFG, "Feature");
    end
end

%% ======================== MATCHED TRAIN / TEST: RAW CNN ========================

if CFG.runMatchedRawCNNPipeline

    fprintf('\n================ MATCHED RAW CNN EVALUATION ================\n');

    if CFG.excludeFs6FromCNN
        keepCNN = ~WindowInfo.isForced6Hz;

        RawX_CNN = RawX(keepCNN);
        Y_CNN = Y(keepCNN);
        WindowInfo_CNN = WindowInfo(keepCNN,:);

        fprintf('RAW CNN: excluding forced-6Hz participants. Windows kept: %d/%d\n', ...
            numel(RawX_CNN), numel(RawX));
        fprintf('RAW CNN: excluded participant keys:\n');
        disp(unique(WindowInfo.participantKey(WindowInfo.isForced6Hz)));
    else
        RawX_CNN = RawX;
        Y_CNN = Y;
        WindowInfo_CNN = WindowInfo;
    end

    rawCnnRes = runMatchedRawCNNParticipantEvaluation(RawX_CNN, Y_CNN, WindowInfo_CNN, CFG);

    writetable(rawCnnRes.results, ...
        fullfile(CFG.outDir, sprintf('MatchedRawCNN_Results_%d.csv', CFG.nMatchedRepeats)));

    writetable(rawCnnRes.predictions, ...
        fullfile(CFG.outDir, sprintf('MatchedRawCNN_Predictions_%d.csv', CFG.nMatchedRepeats)));

    writetable(rawCnnRes.participantScores, ...
        fullfile(CFG.outDir, sprintf('MatchedRawCNN_ParticipantScores_%d.csv', CFG.nMatchedRepeats)));

    writetable(rawCnnRes.summary, ...
        fullfile(CFG.outDir, sprintf('MatchedRawCNN_Summary_%d.csv', CFG.nMatchedRepeats)));

    disp(rawCnnRes.summary);

    if CFG.makePlots
        plotMatchedAucDistributions(rawCnnRes.results, CFG, "RawCNN");
    end
end

fprintf('\nDone. Outputs saved to:\n%s\n', CFG.outDir);


%% =====================================================================
%% =============== MATCHED DONOR / CONTROL FEATURE FUNCTIONS ============
%% =====================================================================

function FeatureTable = computeFeatureTableFromRawWindows_Matched(RawX, Y, WindowInfo, CFG)

    rows = cell(numel(RawX), 1);

    for i = 1:numel(RawX)

        x = RawX{i};
        Fs = WindowInfo.Fs(i);

        feat = extractRespFeaturesFromStereoWindow(x, Fs);
        T = struct2table(feat);

        T.participant = WindowInfo.participant(i);
        T.participantKey = WindowInfo.participantKey(i);
        T.group = WindowInfo.group(i);
        T.subjIndex = WindowInfo.subjIndex(i);
        T.isForced6Hz = WindowInfo.isForced6Hz(i);
        T.donationFractionMid = WindowInfo.donationFractionMid(i);
        T.label = string(Y(i));
        T.labelLate = double(Y(i) == "late");

        rows{i} = T;

        if mod(i,100)==0
            fprintf('Computed features for %d/%d windows.\n', i, numel(RawX));
        end
    end

    FeatureTable = vertcat(rows{:});

    frontVars = {'participant','participantKey','group','subjIndex','isForced6Hz', ...
        'donationFractionMid','label','labelLate'};
    otherVars = setdiff(FeatureTable.Properties.VariableNames, frontVars, 'stable');
    FeatureTable = FeatureTable(:, [frontVars otherVars]);
end

function out = runMatchedFeatureGroupEvaluation(T, CFG)

    rng(CFG.randomSeed);

    T.group = lower(string(T.group));
    T.participantKey = string(T.participantKey);

    donorPids   = unique(T.participantKey(T.group == "donor"));
    controlPids = unique(T.participantKey(T.group == "control"));

    donorPids = donorPids(~ismissing(donorPids));
    controlPids = controlPids(~ismissing(controlPids));

    fprintf('Matched setup participants: donors=%d, controls=%d\n', numel(donorPids), numel(controlPids));

    if isempty(donorPids) || isempty(controlPids)
        error('Need both donor and control participants in FeatureTable.');
    end

    switch lower(string(CFG.trainGroup))
        case "donor"
            if numel(donorPids) <= numel(controlPids)
                error('Donor-trained mode needs more donors than controls. donors=%d controls=%d', ...
                    numel(donorPids), numel(controlPids));
            end
            nTestSamePerRepeat = numel(controlPids);
            nTestOppPerRepeat  = numel(controlPids);
            fprintf('Donor-trained mode: train donors=%d, test held-out donors=%d, test controls=%d each repeat.\n', ...
                numel(donorPids)-nTestSamePerRepeat, nTestSamePerRepeat, nTestOppPerRepeat);

        case "control"
            nTestSamePerRepeat = max(CFG.controlTrain_minTestControls, ...
                round(CFG.controlTrain_testControlFrac * numel(controlPids)));
            nTestSamePerRepeat = min(nTestSamePerRepeat, numel(controlPids)-2);
            nTestOppPerRepeat  = nTestSamePerRepeat;

            if nTestSamePerRepeat < 1
                error('Too few controls for control-trained matched evaluation.');
            end

            fprintf('Control-trained mode: train controls=%d, test held-out controls=%d, matched donors=%d each repeat.\n', ...
                numel(controlPids)-nTestSamePerRepeat, nTestSamePerRepeat, nTestOppPerRepeat);

        otherwise
            error('Unknown CFG.trainGroup: %s', CFG.trainGroup);
    end

    allResults = table();
    allPred = table();
    allScores = table();

    for rep = 1:CFG.nMatchedRepeats

        switch lower(string(CFG.trainGroup))
            case "donor"
                rp = randperm(numel(donorPids));
                testSamePids = donorPids(rp(1:nTestSamePerRepeat));
                trainPids    = donorPids(rp(nTestSamePerRepeat+1:end));
                testOppPids  = controlPids;

                sameLabel = "heldout_donors";
                oppLabel  = "controls";

            case "control"
                rpC = randperm(numel(controlPids));
                testSamePids = controlPids(rpC(1:nTestSamePerRepeat));
                trainPids    = controlPids(rpC(nTestSamePerRepeat+1:end));

                rpD = randperm(numel(donorPids));
                testOppPids = donorPids(rpD(1:nTestOppPerRepeat));

                sameLabel = "heldout_controls";
                oppLabel  = "matched_donors";
        end

        modelPack = fitFeatureModelOnParticipants(T, trainPids, CFG);

        [predSame, scoreSame] = applyFeatureModelToParticipants(T, testSamePids, modelPack, CFG);
        [predOpp,  scoreOpp]  = applyFeatureModelToParticipants(T, testOppPids,  modelPack, CFG);

        resSame = summarizePredictionSet(predSame, rep, string(CFG.trainGroup), sameLabel, numel(trainPids), numel(testSamePids), modelPack.nFeatures);
        resOpp  = summarizePredictionSet(predOpp,  rep, string(CFG.trainGroup), oppLabel,  numel(trainPids), numel(testOppPids),  modelPack.nFeatures);

        allResults = [allResults; resSame; resOpp]; %#ok<AGROW>

        predSame.repeat = repmat(rep, height(predSame), 1);
        predSame.trainGroup = repmat(string(CFG.trainGroup), height(predSame), 1);
        predSame.testSet = repmat(sameLabel, height(predSame), 1);

        predOpp.repeat = repmat(rep, height(predOpp), 1);
        predOpp.trainGroup = repmat(string(CFG.trainGroup), height(predOpp), 1);
        predOpp.testSet = repmat(oppLabel, height(predOpp), 1);

        allPred = [allPred; predSame; predOpp]; %#ok<AGROW>

        scoreSame.repeat = repmat(rep, height(scoreSame), 1);
        scoreSame.trainGroup = repmat(string(CFG.trainGroup), height(scoreSame), 1);
        scoreSame.testSet = repmat(sameLabel, height(scoreSame), 1);

        scoreOpp.repeat = repmat(rep, height(scoreOpp), 1);
        scoreOpp.trainGroup = repmat(string(CFG.trainGroup), height(scoreOpp), 1);
        scoreOpp.testSet = repmat(oppLabel, height(scoreOpp), 1);

        allScores = [allScores; scoreSame; scoreOpp]; %#ok<AGROW>

        if rep == 1 || mod(rep,10)==0 || rep == CFG.nMatchedRepeats
            fprintf('Repeat %3d/%3d | %s AUC=%.3f balAcc=%.3f | %s AUC=%.3f balAcc=%.3f\n', ...
                rep, CFG.nMatchedRepeats, sameLabel, resSame.AUC, resSame.BalancedAccuracy, ...
                oppLabel, resOpp.AUC, resOpp.BalancedAccuracy);
        end
    end

    summary = summarizeMatchedResults(allResults);

    out = struct();
    out.results = allResults;
    out.predictions = allPred;
    out.participantScores = allScores;
    out.summary = summary;
end

function modelPack = fitFeatureModelOnParticipants(T, trainPids, CFG)

    excludeVars = {'participant','participantKey','group','subjIndex','isForced6Hz', ...
        'donationFractionMid','label','labelLate'};
    featureVarsAll = setdiff(T.Properties.VariableNames, excludeVars, 'stable');

    trainRows = ismember(T.participantKey, trainPids);
    Ttr = T(trainRows, :);

    ytr = Ttr.labelLate;
    Ptr = Ttr.participantKey;

    if numel(unique(ytr)) < 2
        error('Training rows have fewer than 2 classes.');
    end

    Xtr = table2array(Ttr(:, featureVarsAll));

    % Training-only feature filtering.
    goodFeat = mean(isnan(Xtr),1) < 0.50;
    Xtr = Xtr(:, goodFeat);
    featureVars = featureVarsAll(goodFeat);

    % Impute train medians before constant-feature detection.
    trainMedian = median(Xtr, 1, 'omitnan');
    trainMedian(isnan(trainMedian)) = 0;
    for j = 1:size(Xtr,2)
        Xtr(isnan(Xtr(:,j)),j) = trainMedian(j);
    end

    sd0 = std(Xtr, 0, 1, 'omitnan');
    goodFeat2 = sd0 > eps;
    Xtr = Xtr(:, goodFeat2);
    featureVars = featureVars(goodFeat2);
    trainMedian = trainMedian(goodFeat2);

    % Participant-level centering using training rows only.
    Xtr = participantCenterFeatures(Xtr, ytr, Ptr, CFG.featureCenteringMode);

    % Nested feature selection using training participants only.
    if isfield(CFG, 'useNestedFeatureSelection') && CFG.useNestedFeatureSelection
        selectedIdx = selectFeaturesByParticipantDelta( ...
            Xtr, ytr, Ptr, CFG.fsK, CFG.fsCorrPrune, CFG.fsCorrThreshold);
    else
        selectedIdx = 1:size(Xtr,2);
    end

    Xtr = Xtr(:, selectedIdx);
    featureVars = featureVars(selectedIdx);
    trainMedian = trainMedian(selectedIdx);

    mu = mean(Xtr, 1, 'omitnan');
    sig = std(Xtr, 0, 1, 'omitnan');
    sig(sig < eps | isnan(sig)) = 1;

    XtrZ = (Xtr - mu) ./ sig;
    XtrZ(isnan(XtrZ)) = 0;

    switch lower(CFG.featureModelType)
        case lower('ridgeLogistic')
            Mdl = fitclinear(XtrZ, ytr, ...
                'Learner','logistic', ...
                'Regularization','ridge', ...
                'Lambda',CFG.lambda, ...
                'ClassNames',[0 1]);

        case lower('linearSVM')
            Mdl = fitclinear(XtrZ, ytr, ...
                'Learner','svm', ...
                'Regularization','ridge', ...
                'Lambda',CFG.lambda, ...
                'ClassNames',[0 1]);

        otherwise
            error('Unknown feature model type: %s', CFG.featureModelType);
    end

    modelPack = struct();
    modelPack.model = Mdl;
    modelPack.featureVars = featureVars;
    modelPack.trainMedian = trainMedian;
    modelPack.mu = mu;
    modelPack.sig = sig;
    modelPack.nFeatures = numel(featureVars);
end

function [Pred, ParticipantScores] = applyFeatureModelToParticipants(T, testPids, modelPack, CFG)

    testRows = ismember(T.participantKey, testPids);
    Tte = T(testRows, :);

    if isempty(Tte)
        Pred = table();
        ParticipantScores = table();
        return;
    end

    Xte = table2array(Tte(:, modelPack.featureVars));
    yte = Tte.labelLate;
    Pte = Tte.participantKey;

    % Impute using TRAIN medians.
    for j = 1:size(Xte,2)
        Xte(isnan(Xte(:,j)),j) = modelPack.trainMedian(j);
    end

    Xte = participantCenterFeatures(Xte, yte, Pte, CFG.featureCenteringMode);

    XteZ = (Xte - modelPack.mu) ./ modelPack.sig;
    XteZ(isnan(XteZ)) = 0;

    [yp, score] = predict(modelPack.model, XteZ);

    Pred = table();
    Pred.participant = Tte.participant;
    Pred.participantKey = Tte.participantKey;
    Pred.group = Tte.group;
    Pred.donationFractionMid = Tte.donationFractionMid;
    Pred.trueLate = yte;
    Pred.predictedLate = yp;
    Pred.lateScore = score(:,2);

    ParticipantScores = computeParticipantScores(Pred);
end

function ResultT = summarizePredictionSet(Pred, rep, trainGroup, testSet, nTrainParticipants, nTestParticipants, nFeatures)

    ResultT = table();
    ResultT.repeat = rep;
    ResultT.trainGroup = trainGroup;
    ResultT.testSet = testSet;
    ResultT.NtrainParticipants = nTrainParticipants;
    ResultT.NtestParticipants = nTestParticipants;
    ResultT.NtestRows = height(Pred);
    ResultT.Nfeatures = nFeatures;

    if isempty(Pred) || numel(unique(Pred.trueLate)) < 2 || all(isnan(Pred.lateScore))
        ResultT.Accuracy = NaN;
        ResultT.BalancedAccuracy = NaN;
        ResultT.AUC = NaN;
        ResultT.SensitivityLate = NaN;
        ResultT.SpecificityEarly = NaN;
        return;
    end

    y = Pred.trueLate;
    yhat = Pred.predictedLate;
    score = Pred.lateScore;

    ResultT.Accuracy = mean(yhat == y, 'omitnan');
    ResultT.SensitivityLate = mean(yhat(y==1) == 1, 'omitnan');
    ResultT.SpecificityEarly = mean(yhat(y==0) == 0, 'omitnan');
    ResultT.BalancedAccuracy = mean([ResultT.SensitivityLate ResultT.SpecificityEarly], 'omitnan');

    try
        [~,~,~,auc] = perfcurve(y, score, 1);
    catch
        auc = NaN;
    end
    ResultT.AUC = auc;
end

function SummaryT = summarizeMatchedResults(R)

    if isempty(R)
        SummaryT = table();
        return;
    end

    testSets = unique(string(R.testSet), 'stable');
    rows = {};

    for i = 1:numel(testSets)
        ts = testSets(i);
        idx = string(R.testSet) == ts;

        rows(end+1,:) = { ...
            ts, sum(idx), ...
            mean(R.AUC(idx), 'omitnan'), median(R.AUC(idx), 'omitnan'), ...
            prctile(R.AUC(idx), 2.5), prctile(R.AUC(idx), 97.5), ...
            mean(R.BalancedAccuracy(idx), 'omitnan'), median(R.BalancedAccuracy(idx), 'omitnan'), ...
            prctile(R.BalancedAccuracy(idx), 2.5), prctile(R.BalancedAccuracy(idx), 97.5), ...
            mean(R.Accuracy(idx), 'omitnan'), median(R.Accuracy(idx), 'omitnan')}; %#ok<AGROW>
    end

    SummaryT = cell2table(rows, 'VariableNames', { ...
        'testSet','Nrepeats', ...
        'meanAUC','medianAUC','AUC_pct025','AUC_pct975', ...
        'meanBalancedAccuracy','medianBalancedAccuracy','BalancedAccuracy_pct025','BalancedAccuracy_pct975', ...
        'meanAccuracy','medianAccuracy'});

    % Add paired difference between the two test sets when available.
    if numel(testSets) == 2
        R1 = R(string(R.testSet)==testSets(1), :);
        R2 = R(string(R.testSet)==testSets(2), :);

        commonRep = intersect(R1.repeat, R2.repeat);
        dAUC = nan(numel(commonRep),1);
        dBal = nan(numel(commonRep),1);
        dAcc = nan(numel(commonRep),1);

        for ii = 1:numel(commonRep)
            r = commonRep(ii);
            a = R1(R1.repeat == r, :);
            b = R2(R2.repeat == r, :);

            if ~isempty(a) && ~isempty(b)
                dAUC(ii) = b.AUC(1) - a.AUC(1);
                dBal(ii) = b.BalancedAccuracy(1) - a.BalancedAccuracy(1);
                dAcc(ii) = b.Accuracy(1) - a.Accuracy(1);
            end
        end

        deltaRow = table();
        deltaRow.testSet = testSets(2) + "_minus_" + testSets(1);
        deltaRow.Nrepeats = numel(commonRep);
        deltaRow.meanAUC = mean(dAUC, 'omitnan');
        deltaRow.medianAUC = median(dAUC, 'omitnan');
        deltaRow.AUC_pct025 = prctile(dAUC, 2.5);
        deltaRow.AUC_pct975 = prctile(dAUC, 97.5);
        deltaRow.meanBalancedAccuracy = mean(dBal, 'omitnan');
        deltaRow.medianBalancedAccuracy = median(dBal, 'omitnan');
        deltaRow.BalancedAccuracy_pct025 = prctile(dBal, 2.5);
        deltaRow.BalancedAccuracy_pct975 = prctile(dBal, 97.5);
        deltaRow.meanAccuracy = mean(dAcc, 'omitnan');
        deltaRow.medianAccuracy = median(dAcc, 'omitnan');

        SummaryT = [SummaryT; deltaRow];
    end
end

function plotMatchedAucDistributions(R, CFG, plotPrefix)

    if nargin < 3 || strlength(string(plotPrefix)) == 0
        plotPrefix = "Matched";
    end

    testSets = unique(string(R.testSet), 'stable');
    if isempty(testSets)
        return;
    end

    fig = figure('Color','w','Position',[200 200 760 480]);
    hold on;
    edges = 0:0.025:1;

    for i = 1:numel(testSets)
        vals = R.AUC(string(R.testSet)==testSets(i));
        vals = vals(isfinite(vals));
        if isempty(vals)
            continue;
        end
        histogram(vals, edges, 'Normalization','probability', 'FaceAlpha',0.35);
        xline(mean(vals,'omitnan'), 'LineWidth',2);
    end

    xlabel('AUC');
    ylabel('Proportion of repeats');
    title(sprintf('Matched early-vs-late model | trainGroup = %s', CFG.trainGroup));
    legend(cellstr(testSets), 'Location','best');
    xlim([0 1]);
    grid on; box off;

    saveas(fig, fullfile(CFG.outDir, sprintf('%s_AUC_hist_%d.png', char(plotPrefix), CFG.nMatchedRepeats)));
    close(fig);
end

%% =====================================================================
%% ============================ FUNCTIONS ===============================
%% =====================================================================

function subjData = loadSubjDataStruct(matFile)
    S = load(matFile);
    names = fieldnames(S);

    if isfield(S, 'subjData')
        subjData = S.subjData;
        return;
    end

    % Otherwise find first struct that looks like subjData.
    for i = 1:numel(names)
        x = S.(names{i});
        if isstruct(x) && numel(x) > 1
            if isfield(x, 'data') && isfield(x, 'in') && isfield(x, 'out')
                subjData = x;
                fprintf('Using variable "%s" as subjData.\n', names{i});
                return;
            end
        end
    end

    error('Could not find subjData-like struct in MAT file.');
end

function [RawX, Y, WindowInfo] = buildRawWindowsFromSubjData(subjData, CFG)

    RawX = {};
    Y = strings(0,1);

    rows = {};

    for s = 1:numel(subjData)

        try
            D = subjData(s).Data;
        catch
            warning('Participant %d skipped: no .data table.', s);
            continue;
        end

        if ~istable(D)
            warning('Participant %d skipped: subjData(%d).data is not a table.', s, s);
            continue;
        end

        pid = inferParticipantID(subjData, s);

        isForced6Hz = ismember(pid, CFG.forceFs6Participants);

        try
            L = getChannelFromTable(D, CFG.leftChannelCandidates);
            R = getChannelFromTable(D, CFG.rightChannelCandidates);
        catch ME
            warning('Participant %s skipped: could not get L/R channels: %s', pid, ME.message);
            continue;
        end

        L = double(L(:));
        R = double(R(:));

        n = min(numel(L), numel(R));
        L = L(1:n);
        R = R(1:n);

        % Fill missing values before any normalization or cutting.
L = fillmissing(L, 'linear', 'EndValues','nearest');
R = fillmissing(R, 'linear', 'EndValues','nearest');

% Optional normalization before cutting donation interval.
if strcmpi(CFG.rawNormalizeStage, 'beforeCut')
    L = normalizeRawVector(L, CFG.rawNormalizeMethod);
    R = normalizeRawVector(R, CFG.rawNormalizeMethod);
end

Fs = getParticipantFs(subjData, s, CFG);

if isForced6Hz
    fprintf('Participant %s: forcing Fs = 6 Hz based on known acquisition rate.\n', pid);
    Fs = 6;
end

        inIdx  = round(double(subjData(s).in));
        outIdx = round(double(subjData(s).out));

        if isempty(inIdx) || isempty(outIdx) || isnan(inIdx) || isnan(outIdx)
            warning('Participant %s skipped: missing in/out.', pid);
            continue;
        end

        if outIdx <= inIdx || inIdx < 1 || outIdx > n
            warning('Participant %s skipped: invalid in/out. in=%g out=%g n=%g', ...
                pid, inIdx, outIdx, n);
            continue;
        end

        % Cut during donation.
        a0 = inIdx + round(CFG.excludeFirstSec * Fs);
        b0 = outIdx - round(CFG.excludeLastSec * Fs);

        if b0 <= a0
            warning('Participant %s skipped: interval too short after edge exclusion.', pid);
            continue;
        end

        Ldur = L(a0:b0);
        Rdur = R(a0:b0);

       
        % Fill missing.
        Ldur = fillmissing(Ldur, 'linear', 'EndValues','nearest');
        Rdur = fillmissing(Rdur, 'linear', 'EndValues','nearest');

        % Detrend.
        if CFG.doDetrend
            Ldur = detrend(Ldur);
            Rdur = detrend(Rdur);
        end

        % Bandpass.
        if CFG.doBandpass
            Ldur = safeBandpass(Ldur, Fs, CFG.bandpassHz);
            Rdur = safeBandpass(Rdur, Fs, CFG.bandpassHz);
        end

        % Resample all participants to target Fs.
        if abs(Fs - CFG.targetFs) > 1e-6
            fprintf('Participant %s: resampling from %.3g Hz to %.3g Hz.\n', ...
                pid, Fs, CFG.targetFs);

            Ldur = resample(Ldur(:), CFG.targetFs, Fs);
            Rdur = resample(Rdur(:), CFG.targetFs, Fs);
            Fs = CFG.targetFs;
        end

        Ldur = Ldur(:)';
        Rdur = Rdur(:)';

        % Participant-level robust normalization across donation.
        % Optional normalization after cutting donation interval.
if strcmpi(CFG.rawNormalizeStage, 'afterCut')
    Ldur = normalizeRawVector(Ldur, CFG.rawNormalizeMethod);
    Rdur = normalizeRawVector(Rdur, CFG.rawNormalizeMethod);
end


        nDur = numel(Ldur);

        winN  = round(CFG.windowSec * Fs);
        stepN = max(1, round(winN * (1 - CFG.overlap)));

        if nDur < winN
            warning('Participant %s skipped: donation interval shorter than one window.', pid);
            continue;
        end

        starts = 1:stepN:(nDur - winN + 1);

        for k = 1:numel(starts)

            wStart = starts(k);
            wEnd   = wStart + winN - 1;

            fracMid = ((wStart + wEnd)/2 - 1) / max(1, nDur - 1);

            if fracMid >= CFG.earlyRange(1) && fracMid <= CFG.earlyRange(2)
                lab = "early";
            elseif fracMid >= CFG.lateRange(1) && fracMid <= CFG.lateRange(2)
                lab = "late";
            else
                continue;
            end

Lw = Ldur(wStart:wEnd);
Rw = Rdur(wStart:wEnd);

% Optional per-window normalization.
if strcmpi(CFG.rawNormalizeStage, 'perWindow')
    Lw = normalizeRawVector(Lw, CFG.rawNormalizeMethod);
    Rw = normalizeRawVector(Rw, CFG.rawNormalizeMethod);
end

% Channel mode.
switch lower(CFG.channelMode)
    case 'stereo'
        x = [Lw; Rw];

    case 'sum'
        switch lower(CFG.sumChannelMethod)
            case 'mean'
                xsum = mean([Lw; Rw], 1, 'omitnan');
            case 'sum'
                xsum = Lw + Rw;
            otherwise
                error('Unknown CFG.sumChannelMethod: %s', CFG.sumChannelMethod);
        end

        x = xsum;

    otherwise
        error('Unknown CFG.channelMode: %s', CFG.channelMode);
end

            RawX{end+1,1} = x; %#ok<AGROW>
            Y(end+1,1) = lab; %#ok<AGROW>

            globalStartApprox = a0 + round((wStart - 1) * getParticipantFs(subjData, s, CFG) / Fs);
            globalEndApprox   = a0 + round((wEnd   - 1) * getParticipantFs(subjData, s, CFG) / Fs);

            rows(end+1,:) = { ...
    pid, s, Fs, isForced6Hz, ...
    inIdx, outIdx, ...
    globalStartApprox, globalEndApprox, ...
    wStart, wEnd, ...
    fracMid, lab}; %#ok<AGROW>

        end
    end

WindowInfo = cell2table(rows, ...
    'VariableNames', { ...
    'participant','subjIndex','Fs','isForced6Hz', ...
    'donationInSample','donationOutSample', ...
    'globalWindowStartApprox','globalWindowEndApprox', ...
    'windowStartInCut','windowEndInCut', ...
    'donationFractionMid','label'});

    WindowInfo.participant = normalizeParticipantID(WindowInfo.participant);

    Y = categorical(Y, ["early","late"]);
end

function pid = inferParticipantID(subjData, s)
% Infer participant ID using the true participant code.
% IMPORTANT: This uses subjData(s).code first, not the serial index s.

    % First priority: subjData(s).code
    if isfield(subjData, 'code')
        val = subjData(s).code;

        if ~isempty(val)
            pid = normalizeParticipantID(string(val));
            return;
        end
    end

    % Fallback options if code does not exist.
    candidates = {'participant','participantID','subject_id','subjectID','subjID','id'};

    for i = 1:numel(candidates)
        if isfield(subjData, candidates{i})
            val = subjData(s).(candidates{i});

            if ~isempty(val)
                pid = normalizeParticipantID(string(val));
                return;
            end
        end
    end

    % Last fallback only.
    warning('Using serial index %d as participant ID because no subjData.code was found.', s);
    pid = sprintf('%03d', s);
end


function Fs = getParticipantFs(subjData, s, CFG)

    candidates = {'Fs','fs','samplingRate','sampleRate'};

    Fs = CFG.defaultFs;

    for i = 1:numel(candidates)
        if isfield(subjData, candidates{i})
            val = subjData(s).(candidates{i});
            if ~isempty(val) && isnumeric(val) && isfinite(val)
                Fs = double(val);
                return;
            end
        end
    end
end

function x = getChannelFromTable(T, candidates)

    vars = T.Properties.VariableNames;

    for i = 1:numel(candidates)
        idx = find(strcmpi(vars, candidates{i}), 1);
        if ~isempty(idx)
            x = T.(vars{idx});
            return;
        end
    end

    % Fallback: if no channel names match, take numeric columns.
    numericMask = varfun(@isnumeric, T, 'OutputFormat','uniform');
    numericVars = vars(numericMask);

    if numel(numericVars) >= 2
        warning('Channel name not found. Using numeric column fallback.');
        if contains(lower(candidates{1}), 'l')
            x = T.(numericVars{1});
        else
            x = T.(numericVars{2});
        end
        return;
    end

    error('No matching channel found.');
end

function y = safeBandpass(x, Fs, bandHz)

    x = double(x(:));

    if numel(x) < round(Fs * 10)
        y = x;
        return;
    end

    nyq = Fs / 2;

    if bandHz(2) >= nyq
        bandHz(2) = nyq * 0.95;
    end

    if bandHz(1) <= 0
        bandHz(1) = 0.01;
    end

    try
        y = bandpass(x, bandHz, Fs);
    catch
        [b,a] = butter(2, bandHz ./ nyq, 'bandpass');
        y = filtfilt(b,a,x);
    end

    y = y(:);
end

function z = robustZ(x)

    x = double(x);
    med = median(x, 'omitnan');
    madVal = mad(x, 1);

    if isnan(madVal) || madVal < eps
        madVal = std(x, 0, 'omitnan');
    end

    if isnan(madVal) || madVal < eps
        madVal = 1;
    end

    z = (x - med) ./ madVal;
end

function pid = normalizeParticipantID(pid)

    pid = string(pid);
    pid = strtrim(pid);
    pid = regexprep(pid, '[^\d]', '');

    out = strings(size(pid));

    for i = 1:numel(pid)
        if strlength(pid(i)) == 0
            out(i) = missing;
        else
            out(i) = sprintf('%03d', str2double(pid(i)));
        end
    end

    pid = out;
end

%% =====================================================================
%% ======================= FEATURE PIPELINE =============================
%% =====================================================================

function FeatureTable = computeFeatureTableFromRawWindows(RawX, Y, WindowInfo, CFG)

    rows = cell(numel(RawX), 1);

    for i = 1:numel(RawX)

        x = RawX{i};
        Fs = WindowInfo.Fs(i);

feat = extractRespFeaturesFromStereoWindow(x, Fs);

        T = struct2table(feat);

        T.participant = WindowInfo.participant(i);
        T.subjIndex = WindowInfo.subjIndex(i);
        T.donationFractionMid = WindowInfo.donationFractionMid(i);
        T.label = string(Y(i));
        T.labelLate = double(Y(i) == "late");

        rows{i} = T;

        if mod(i,100)==0
            fprintf('Computed features for %d/%d windows.\n', i, numel(RawX));
        end
    end

    FeatureTable = vertcat(rows{:});

    % Move identifiers to front.
    frontVars = {'participant','subjIndex','donationFractionMid','label','labelLate'};
    otherVars = setdiff(FeatureTable.Properties.VariableNames, frontVars, 'stable');
    FeatureTable = FeatureTable(:, [frontVars otherVars]);
end

function feat = extractRespFeaturesFromStereoWindow(x, Fs)

    feat = struct();

    % x can be:
    %   2 x time = stereo L/R
    %   1 x time = summed channel

    if size(x,1) == 1

        S = double(x(1,:));

        % Sum-channel features only.
        feat = addChannelFeatures(feat, S, Fs, 'SUM');

    elseif size(x,1) == 2

        L = double(x(1,:));
        R = double(x(2,:));

        feat = addChannelFeatures(feat, L, Fs, 'L');
        feat = addChannelFeatures(feat, R, Fs, 'R');

        % Nostril coordination features.
        good = isfinite(L) & isfinite(R);

        if sum(good) > 5
            Lg = L(good);
            Rg = R(good);

            feat.NC_corrLR = corr(Lg(:), Rg(:), 'Rows','complete');

            [maxCorr, lagSec] = maxLaggedCorr(Lg, Rg, Fs, 2.0);
            feat.NC_maxXcorrLR = maxCorr;
            feat.NC_lagAtMaxCorrSec = lagSec;

            denom = abs(Lg) + abs(Rg) + eps;
            LI = (Lg - Rg) ./ denom;

            feat.NC_meanLI = mean(LI, 'omitnan');
            feat.NC_stdLI = std(LI, 0, 'omitnan');
            feat.NC_absMeanLI = mean(abs(LI), 'omitnan');

            rmsL = rms(Lg);
            rmsR = rms(Rg);
            feat.NC_rmsAmplitudeLI = (rmsL - rmsR) / (rmsL + rmsR + eps);

            ampL = prctile(Lg,95) - prctile(Lg,5);
            ampR = prctile(Rg,95) - prctile(Rg,5);
            feat.NC_robustAmplitudeLI = (ampL - ampR) / (ampL + ampR + eps);

            try
                envL = abs(hilbert(Lg));
                envR = abs(hilbert(Rg));
                feat.NC_envelopeCorrLR = corr(envL(:), envR(:), 'Rows','complete');
            catch
                feat.NC_envelopeCorrLR = nan;
            end

        else
            feat.NC_corrLR = nan;
            feat.NC_maxXcorrLR = nan;
            feat.NC_lagAtMaxCorrSec = nan;
            feat.NC_meanLI = nan;
            feat.NC_stdLI = nan;
            feat.NC_absMeanLI = nan;
            feat.NC_rmsAmplitudeLI = nan;
            feat.NC_robustAmplitudeLI = nan;
            feat.NC_envelopeCorrLR = nan;
        end

    else
        error('Expected raw window to have 1 or 2 channels. Got %d.', size(x,1));
    end
end

function feat = addChannelFeatures(feat, x, Fs, prefix)

    x = double(x(:)');
    x = fillmissing(x, 'linear', 'EndValues','nearest');

    % Basic amplitude features.
    feat.([prefix '_mean']) = mean(x, 'omitnan');
    feat.([prefix '_std']) = std(x, 0, 'omitnan');
    feat.([prefix '_rms']) = rms(x);
    feat.([prefix '_mad']) = mad(x, 1);
    feat.([prefix '_iqr']) = iqr(x);
    feat.([prefix '_range95_5']) = prctile(x,95) - prctile(x,5);
    feat.([prefix '_skewness']) = skewness(x);
    feat.([prefix '_kurtosis']) = kurtosis(x);

    % Derivative / flow-like features.
    dx = [0 diff(x)] * Fs;

    feat.([prefix '_flow_std']) = std(dx, 0, 'omitnan');
    feat.([prefix '_flow_rms']) = rms(dx);
    feat.([prefix '_peakInspiratoryProxy']) = prctile(dx, 95);
    feat.([prefix '_peakExpiratoryProxy']) = prctile(dx, 5);
    feat.([prefix '_flowRange95_5']) = prctile(dx,95) - prctile(dx,5);

    % Spectral features.
    spec = spectralFeatures(x, Fs);

    names = fieldnames(spec);
    for i = 1:numel(names)
        feat.([prefix '_' names{i}]) = spec.(names{i});
    end

    % Breath-cycle-like features.
    cyc = cycleFeatures(x, Fs);

    names = fieldnames(cyc);
    for i = 1:numel(names)
        feat.([prefix '_' names{i}]) = cyc.(names{i});
    end
end

function spec = spectralFeatures(x, Fs)

    spec = struct();

    x = x(:);
    x = x - mean(x, 'omitnan');

    if numel(x) < Fs * 5
        spec.specPeakFreq = nan;
        spec.specPeakPower = nan;
        spec.power_005_015 = nan;
        spec.power_015_040 = nan;
        spec.power_040_120 = nan;
        spec.power_005_120 = nan;
        spec.spectralEntropy = nan;
        return;
    end

    try
        [Pxx,F] = pwelch(x, [], [], [], Fs);

        bandMask = F >= 0.05 & F <= 1.20;

        if any(bandMask)
            [pk, idx] = max(Pxx(bandMask));
            Fb = F(bandMask);

            spec.specPeakFreq = Fb(idx);
            spec.specPeakPower = pk;
        else
            spec.specPeakFreq = nan;
            spec.specPeakPower = nan;
        end

        spec.power_005_015 = integrateBandPower(F, Pxx, [0.05 0.15]);
        spec.power_015_040 = integrateBandPower(F, Pxx, [0.15 0.40]);
        spec.power_040_120 = integrateBandPower(F, Pxx, [0.40 1.20]);
        spec.power_005_120 = integrateBandPower(F, Pxx, [0.05 1.20]);

        P = Pxx(bandMask);
        P = P ./ (sum(P) + eps);
        spec.spectralEntropy = -sum(P .* log(P + eps));

    catch
        spec.specPeakFreq = nan;
        spec.specPeakPower = nan;
        spec.power_005_015 = nan;
        spec.power_015_040 = nan;
        spec.power_040_120 = nan;
        spec.power_005_120 = nan;
        spec.spectralEntropy = nan;
    end
end

function bp = integrateBandPower(F, Pxx, band)

    mask = F >= band(1) & F <= band(2);

    if sum(mask) < 2
        bp = nan;
    else
        bp = trapz(F(mask), Pxx(mask));
    end
end

function cyc = cycleFeatures(x, Fs)

    cyc = struct();

    x = double(x(:)');
    x = smoothdata(x, 'movmean', max(3, round(0.25 * Fs)));
    x = robustZ(x);

    minDist = max(1, round(1.5 * Fs));

    try
        [pks, pkLocs] = findpeaks(x, ...
            'MinPeakDistance', minDist, ...
            'MinPeakProminence', 0.25);

        [trsNeg, trLocs] = findpeaks(-x, ...
            'MinPeakDistance', minDist, ...
            'MinPeakProminence', 0.25);

        trs = -trsNeg;
    catch
        pks = [];
        pkLocs = [];
        trs = [];
        trLocs = [];
    end

    nPeaks = numel(pkLocs);
    nTroughs = numel(trLocs);

    cyc.nDetectedPeaks = nPeaks;
    cyc.nDetectedTroughs = nTroughs;

    durationMin = numel(x) / Fs / 60;
    cyc.breathRatePeakBased = nPeaks / max(durationMin, eps);

    if nPeaks >= 3
        ibi = diff(pkLocs) / Fs;
        cyc.meanIBI = mean(ibi, 'omitnan');
        cyc.stdIBI = std(ibi, 0, 'omitnan');
        cyc.cvIBI = std(ibi, 0, 'omitnan') / (mean(ibi, 'omitnan') + eps);
        cyc.medianIBI = median(ibi, 'omitnan');
    else
        cyc.meanIBI = nan;
        cyc.stdIBI = nan;
        cyc.cvIBI = nan;
        cyc.medianIBI = nan;
    end

    % Rise/fall durations around peaks.
    riseDur = [];
    fallDur = [];
    cycleAmp = [];

    for i = 1:numel(pkLocs)
        prevTr = trLocs(find(trLocs < pkLocs(i), 1, 'last'));
        nextTr = trLocs(find(trLocs > pkLocs(i), 1, 'first'));

        if ~isempty(prevTr) && ~isempty(nextTr)
            riseDur(end+1) = (pkLocs(i) - prevTr) / Fs; %#ok<AGROW>
            fallDur(end+1) = (nextTr - pkLocs(i)) / Fs; %#ok<AGROW>
            cycleAmp(end+1) = pks(i) - mean([x(prevTr), x(nextTr)]); %#ok<AGROW>
        end
    end

    if ~isempty(riseDur)
        cyc.meanRiseDuration = mean(riseDur, 'omitnan');
        cyc.meanFallDuration = mean(fallDur, 'omitnan');
        cyc.riseDutyCycle = mean(riseDur ./ (riseDur + fallDur + eps), 'omitnan');
        cyc.meanCycleAmplitude = mean(cycleAmp, 'omitnan');
        cyc.cvCycleAmplitude = std(cycleAmp, 0, 'omitnan') / ...
            (mean(abs(cycleAmp), 'omitnan') + eps);
    else
        cyc.meanRiseDuration = nan;
        cyc.meanFallDuration = nan;
        cyc.riseDutyCycle = nan;
        cyc.meanCycleAmplitude = nan;
        cyc.cvCycleAmplitude = nan;
    end
end

function [maxCorr, lagSec] = maxLaggedCorr(x, y, Fs, maxLagSec)

    x = x(:) - mean(x, 'omitnan');
    y = y(:) - mean(y, 'omitnan');

    maxLag = round(maxLagSec * Fs);

    try
        [c,lags] = xcorr(x, y, maxLag, 'coeff');
        [maxCorr, idx] = max(c);
        lagSec = lags(idx) / Fs;
    catch
        maxCorr = nan;
        lagSec = nan;
    end
end

function res = runFeatureLOSO(T, CFG)

    y = T.labelLate;
    participant = T.participant;

    excludeVars = {'participant','subjIndex','donationFractionMid','label','labelLate'};
    vars = T.Properties.VariableNames;
    featureVars = setdiff(vars, excludeVars, 'stable');

    X = table2array(T(:, featureVars));

    % Remove bad rows.
    goodRow = mean(isnan(X),2) < 0.50;
    X = X(goodRow,:);
    y = y(goodRow);
    participant = participant(goodRow);
    donationFractionMid = T.donationFractionMid(goodRow);

    % Remove bad features.
    goodFeat = mean(isnan(X),1) < 0.50;
    X = X(:, goodFeat);
    featureVars = featureVars(goodFeat);

    % Remove constant features after rough imputation.
    for j = 1:size(X,2)
        med = median(X(:,j), 'omitnan');
        if isnan(med)
            med = 0;
        end
        X(isnan(X(:,j)),j) = med;
    end

    sd = std(X, 0, 1, 'omitnan');
    goodFeat = sd > eps;
    X = X(:, goodFeat);
    featureVars = featureVars(goodFeat);

    fprintf('Feature model uses %d windows and %d features.\n', size(X,1), size(X,2));

    pids = unique(participant);

    lateScore = nan(size(y));
    yhat = nan(size(y));

    for p = 1:numel(pids)

        testIdx = participant == pids(p);
        trainIdx = ~testIdx;

        Xtr = X(trainIdx,:);
        Xte = X(testIdx,:);
        ytr = y(trainIdx);

        Ptr = participant(trainIdx);
        Pte = participant(testIdx);
        YteKnown = y(testIdx);

        % Participant centering.
        Xtr = participantCenterFeatures(Xtr, ytr, Ptr, CFG.featureCenteringMode);
        Xte = participantCenterFeatures(Xte, YteKnown, Pte, CFG.featureCenteringMode);

        % Optional nested feature selection.
% Important: selection is done using training participants only.

       % Nested feature selection using training participants only.
if isfield(CFG, 'useNestedFeatureSelection') && CFG.useNestedFeatureSelection

    selectedIdx = selectFeaturesByParticipantDelta( ...
        Xtr, ytr, Ptr, CFG.fsK, CFG.fsCorrPrune, CFG.fsCorrThreshold);

    Xtr = Xtr(:, selectedIdx);
    Xte = Xte(:, selectedIdx);

    selectedFeatureNames = featureVars(selectedIdx);

    fprintf('Selected %d features for participant %s:\n', ...
        numel(selectedIdx), pids(p));
    disp(string(selectedFeatureNames(:))');

else
    selectedIdx = 1:size(Xtr,2);
end

% Train-only standardization after feature selection.
mu = mean(Xtr, 1, 'omitnan');
sig = std(Xtr, 0, 1, 'omitnan');
sig(sig < eps | isnan(sig)) = 1;

Xtr = (Xtr - mu) ./ sig;
Xte = (Xte - mu) ./ sig;

        Xtr(isnan(Xtr)) = 0;
        Xte(isnan(Xte)) = 0;

        switch lower(CFG.featureModelType)

            case lower('ridgeLogistic')
                Mdl = fitclinear(Xtr, ytr, ...
                    'Learner','logistic', ...
                    'Regularization','ridge', ...
                    'Lambda',CFG.lambda, ...
                    'ClassNames',[0 1]);

            case lower('linearSVM')
                Mdl = fitclinear(Xtr, ytr, ...
                    'Learner','svm', ...
                    'Regularization','ridge', ...
                    'Lambda',CFG.lambda, ...
                    'ClassNames',[0 1]);

            otherwise
                error('Unknown feature model type.');
        end

        [yp, score] = predict(Mdl, Xte);

        yhat(testIdx) = yp;

        % For logistic this is usually posterior-like; for SVM it is score-like.
        lateScore(testIdx) = score(:,2);

        fprintf('LOSO %d/%d | participant %s | test windows %d\n', ...
            p, numel(pids), pids(p), sum(testIdx));
    end

    [~,~,~,auc] = perfcurve(y, lateScore, 1);

    predLate = lateScore >= median(lateScore, 'omitnan');
    sens = mean(predLate(y==1) == 1, 'omitnan');
    spec = mean(predLate(y==0) == 0, 'omitnan');
    balAcc = mean([sens spec], 'omitnan');

    predictions = table(participant, donationFractionMid, y, yhat, lateScore, ...
        'VariableNames', {'participant','donationFractionMid','trueLate','predictedLate','lateScore'});

    participantScores = computeParticipantScores(predictions);

    res = struct();
    res.predictions = predictions;
    res.participantScores = participantScores;
    res.auc = auc;
    res.balancedAccuracy = balAcc;
    res.featureVars = featureVars;
end

function Xc = participantCenterFeatures(X, Y, P, mode)

    Xc = X;

    switch lower(mode)

        case 'none'
            return;

        case lower('participantDonationMean')
            pids = unique(P);
            for i = 1:numel(pids)
                idx = P == pids(i);
                mu = mean(X(idx,:), 1, 'omitnan');
                Xc(idx,:) = X(idx,:) - mu;
            end

        case lower('participantDonationZ')
            pids = unique(P);
            for i = 1:numel(pids)
                idx = P == pids(i);
                mu = mean(X(idx,:), 1, 'omitnan');
                sig = std(X(idx,:), 0, 1, 'omitnan');
                sig(sig < eps | isnan(sig)) = 1;
                Xc(idx,:) = (X(idx,:) - mu) ./ sig;
            end

        case lower('participantEarlyMean')
            pids = unique(P);
            for i = 1:numel(pids)
                idx = P == pids(i);
                earlyIdx = idx & Y == 0;

                if any(earlyIdx)
                    mu = mean(X(earlyIdx,:), 1, 'omitnan');
                else
                    mu = mean(X(idx,:), 1, 'omitnan');
                end

                Xc(idx,:) = X(idx,:) - mu;
            end

        otherwise
            error('Unknown participant centering mode.');
    end
end

%% =====================================================================
%% ======================== RAW CNN PIPELINE ============================
%% =====================================================================

function res = runRawCNNGroupedCV(RawX, Y, WindowInfo, CFG)

    participant = WindowInfo.participant;
    pids = unique(participant);

    rng(1);

    switch lower(CFG.deepCVMode)

        case lower('groupKfold')
            pidsShuf = pids(randperm(numel(pids)));
            foldByPid = zeros(numel(pidsShuf),1);

            for i = 1:numel(pidsShuf)
                foldByPid(i) = mod(i-1, CFG.K) + 1;
            end

            nFolds = CFG.K;

        otherwise
            error('Currently implemented deep CV mode: groupKfold.');
    end

    lateScore = nan(numel(Y),1);
    predLabel = strings(numel(Y),1);

    for fold = 1:nFolds

        testPids = pidsShuf(foldByPid == fold);
        testIdx = ismember(participant, testPids);
        trainIdx = ~testIdx;

        XTrain = RawX(trainIdx);
        YTrain = Y(trainIdx);

        XTest = RawX(testIdx);
        YTest = Y(testIdx); %#ok<NASGU>

        [XTrain, YTrain] = balanceClasses(XTrain, YTrain);

        fprintf('\nCNN fold %d/%d | train windows %d | test windows %d | test participants %d\n', ...
            fold, nFolds, numel(XTrain), sum(testIdx), numel(testPids));

minSeqLen = min(cellfun(@(x) size(x,2), XTrain));

numInputChannels = size(XTrain{1}, 1);
minSeqLen = min(cellfun(@(x) size(x,2), XTrain));

fprintf('CNN input: channels=%d | minSeqLen=%d\n', numInputChannels, minSeqLen);

layers = makeRawCNNLayers(numInputChannels, minSeqLen);


        options = trainingOptions('adam', ...
            'MaxEpochs', CFG.maxEpochs, ...
            'MiniBatchSize', CFG.miniBatchSize, ...
            'InitialLearnRate', CFG.initialLearnRate, ...
            'Shuffle','every-epoch', ...
            'Verbose',false, ...
            'Plots','none');

        seqLens = cellfun(@(x) size(x,2), XTrain);
nChannels = cellfun(@(x) size(x,1), XTrain);

fprintf('CNN input check: minLen=%d | maxLen=%d | unique channels=%s\n', ...
    min(seqLens), max(seqLens), mat2str(unique(nChannels)));

        net = trainNetwork(XTrain, YTrain, layers, options);

        scores = predict(net, XTest, 'MiniBatchSize', CFG.miniBatchSize);
        yPred = classify(net, XTest, 'MiniBatchSize', CFG.miniBatchSize);

        cats = categories(YTrain);
        lateCol = find(strcmp(cats, 'late'));

        if isempty(lateCol)
            error('Could not find late category.');
        end

        foldRows = find(testIdx);
        lateScore(foldRows) = scores(:, lateCol);
        predLabel(foldRows) = string(yPred);
    end

    trueLate = double(Y == "late");

    [~,~,~,auc] = perfcurve(trueLate, lateScore, 1);

    predLate = lateScore >= 0.5;
    sens = mean(predLate(trueLate==1) == 1, 'omitnan');
    spec = mean(predLate(trueLate==0) == 0, 'omitnan');
    balAcc = mean([sens spec], 'omitnan');

    predictions = table( ...
        participant, ...
        WindowInfo.donationFractionMid, ...
        trueLate, ...
        predLabel, ...
        lateScore, ...
        'VariableNames', {'participant','donationFractionMid','trueLate','predictedLabel','lateScore'});

    participantScores = computeParticipantScores(predictions);

    res = struct();
    res.predictions = predictions;
    res.participantScores = participantScores;
    res.auc = auc;
    res.balancedAccuracy = balAcc;
end

function layers = makeRawCNNLayers(numInputChannels, minSeqLen)

    layers = [
        sequenceInputLayer(numInputChannels, ...
            'Name','input', ...
            'MinLength', minSeqLen)

        convolution1dLayer(15, 16, ...
            'Padding','same', ...
            'Name','conv1')
        batchNormalizationLayer('Name','bn1')
        reluLayer('Name','relu1')
        maxPooling1dLayer(2, ...
            'Stride',2, ...
            'Name','pool1')

        convolution1dLayer(15, 32, ...
            'Padding','same', ...
            'Name','conv2')
        batchNormalizationLayer('Name','bn2')
        reluLayer('Name','relu2')
        maxPooling1dLayer(2, ...
            'Stride',2, ...
            'Name','pool2')

        convolution1dLayer(15, 64, ...
            'Padding','same', ...
            'Name','conv3')
        batchNormalizationLayer('Name','bn3')
        reluLayer('Name','relu3')

        globalAveragePooling1dLayer('Name','globalAverage')

        fullyConnectedLayer(2, 'Name','fc')
        softmaxLayer('Name','softmax')
        classificationLayer('Name','classification')
    ];
end

function [Xb, Yb] = balanceClasses(X, Y)

    idxEarly = find(Y == "early");
    idxLate  = find(Y == "late");

    n = min(numel(idxEarly), numel(idxLate));

    rng(1);

    idx = [ ...
        idxEarly(randperm(numel(idxEarly), n)); ...
        idxLate(randperm(numel(idxLate), n))];

    idx = idx(randperm(numel(idx)));

    Xb = X(idx);
    Yb = Y(idx);
end

%% =====================================================================
%% ===================== PARTICIPANT SCORES / PLOTS =====================
%% =====================================================================

function participantScores = computeParticipantScores(Pred)

    if isempty(Pred)
        participantScores = table();
        return;
    end

    if ismember('participantKey', Pred.Properties.VariableNames)
        idVar = 'participantKey';
    else
        idVar = 'participant';
    end

    pids = unique(string(Pred.(idVar)));
    pids = pids(~ismissing(pids));

    rows = {};

    for i = 1:numel(pids)

        idx = string(Pred.(idVar)) == pids(i);

        earlyScore = mean(Pred.lateScore(idx & Pred.trueLate == 0), 'omitnan');
        lateScore  = mean(Pred.lateScore(idx & Pred.trueLate == 1), 'omitnan');

        respBloodLossScore = lateScore - earlyScore;

        if ismember('group', Pred.Properties.VariableNames)
            groupVal = unique(string(Pred.group(idx)));
            groupVal = groupVal(1);
        else
            groupVal = "";
        end

        if ismember('participant', Pred.Properties.VariableNames)
            participantVal = unique(string(Pred.participant(idx)));
            participantVal = participantVal(1);
        else
            participantVal = pids(i);
        end

        rows(end+1,:) = { ...
            pids(i), ...
            participantVal, ...
            groupVal, ...
            earlyScore, ...
            lateScore, ...
            respBloodLossScore, ...
            sum(idx & Pred.trueLate == 0), ...
            sum(idx & Pred.trueLate == 1)}; %#ok<AGROW>
    end

    participantScores = cell2table(rows, ...
        'VariableNames', { ...
        'participantKey', ...
        'participant', ...
        'group', ...
        'meanLateScoreEarlyWindows', ...
        'meanLateScoreLateWindows', ...
        'respBloodLossScore', ...
        'nEarlyWindows', ...
        'nLateWindows'});
end

function corrTable = correlateScoresWithMetadata(scoreTable, CFG)

    M = readtable(CFG.metaFile, 'VariableNamingRule','preserve');
    M = standardizeMetadataParticipant(M);

    J = innerjoin(scoreTable, M, 'Keys','participant');

    rows = {};

    for i = 1:numel(CFG.outcomeVars)

        v = CFG.outcomeVars{i};

        if ~ismember(v, J.Properties.VariableNames)
            continue;
        end

        y = J.(v);

        if ~isnumeric(y)
            y = str2double(string(y));
        end

        x = J.respBloodLossScore;

        good = isfinite(x) & isfinite(y);

        if sum(good) >= 8
            [rho, pSpearman] = corr(x(good), y(good), 'Type','Spearman');
            [r, pPearson] = corr(x(good), y(good), 'Type','Pearson');
        else
            rho = nan;
            pSpearman = nan;
            r = nan;
            pPearson = nan;
        end

        rows(end+1,:) = {v, sum(good), rho, pSpearman, r, pPearson}; %#ok<AGROW>
    end

    corrTable = cell2table(rows, ...
        'VariableNames', { ...
        'outcome', ...
        'N', ...
        'SpearmanRho', ...
        'SpearmanP', ...
        'PearsonR', ...
        'PearsonP'});
end

function M = standardizeMetadataParticipant(M)

    vars = M.Properties.VariableNames;

    candidates = {'participant','subject_id','subjectID','participantID','subjID'};

    found = '';

    for i = 1:numel(candidates)
        idx = find(strcmpi(vars, candidates{i}), 1);
        if ~isempty(idx)
            found = vars{idx};
            break;
        end
    end

    if isempty(found)
        error('Could not find participant ID column in metadata.');
    end

    M.participant = normalizeParticipantID(string(M.(found)));
end

function plotROC(trueLate, lateScore, savePath, plotTitle)

    [Xroc, Yroc, ~, AUC] = perfcurve(trueLate, lateScore, 1);

    fig = figure('Color','w');
    plot(Xroc, Yroc, 'LineWidth', 1.8);
    hold on;
    plot([0 1], [0 1], '--');
    xlabel('False positive rate');
    ylabel('True positive rate');
    title(sprintf('%s ROC | AUC = %.3f', plotTitle, AUC));
    axis square;
    grid on;

    saveas(fig, savePath);
    close(fig);
end

function plotParticipantScores(scoreTable, savePath, plotTitle)

    fig = figure('Color','w');
    histogram(scoreTable.respBloodLossScore);
    xlabel('Late-window score minus early-window score');
    ylabel('Number of participants');
    title(plotTitle);
    grid on;

    saveas(fig, savePath);
    close(fig);
end

function y = normalizeRawVector(x, method)

    x = double(x);

    switch lower(method)

        case lower('robustZ')
            medVal = median(x, 'omitnan');
            madVal = mad(x, 1);

            if isnan(madVal) || madVal < eps
                madVal = std(x, 0, 'omitnan');
            end

            if isnan(madVal) || madVal < eps
                madVal = 1;
            end

            y = (x - medVal) ./ madVal;

        case lower('zscore')
            mu = mean(x, 'omitnan');
            sd = std(x, 0, 'omitnan');

            if isnan(sd) || sd < eps
                sd = 1;
            end

            y = (x - mu) ./ sd;

        otherwise
            error('Unknown raw normalization method: %s', method);
    end
end


function selectedIdx = selectFeaturesByParticipantDelta(Xtr, ytr, Ptr, fsK, doCorrPrune, corrThreshold)
% Select features using only training data.
%
% For each feature:
%   participant delta = mean(late windows) - mean(early windows)
%
% Then rank features by signrank p-value across training participants.
% This avoids treating overlapping windows as independent observations.

    pids = unique(Ptr);
    nP = numel(pids);
    nF = size(Xtr,2);

    deltaMat = nan(nP, nF);

    for p = 1:nP
        idxP = Ptr == pids(p);

        idxEarly = idxP & ytr == 0;
        idxLate  = idxP & ytr == 1;

        if any(idxEarly) && any(idxLate)
            earlyMean = mean(Xtr(idxEarly,:), 1, 'omitnan');
            lateMean  = mean(Xtr(idxLate,:), 1, 'omitnan');
            deltaMat(p,:) = lateMean - earlyMean;
        end
    end

    pVals = nan(1,nF);
    effectSize = nan(1,nF);

    for f = 1:nF
        d = deltaMat(:,f);
        d = d(isfinite(d));

        if numel(d) >= 8 && std(d,0,'omitnan') > eps
            try
                pVals(f) = signrank(d, 0);
            catch
                pVals(f) = nan;
            end

            effectSize(f) = median(d, 'omitnan') / (mad(d,1) + eps);
        end
    end

    valid = isfinite(pVals);
    candidateIdx = find(valid);

    if isempty(candidateIdx)
        warning('Feature selection found no valid features. Using all features.');
        selectedIdx = 1:nF;
        return;
    end

    rankingMat = [pVals(candidateIdx(:))', -abs(effectSize(candidateIdx(:)))'];
    [~, ord] = sortrows(rankingMat);

    rankedIdx = candidateIdx(ord);

    if doCorrPrune
        selectedIdx = [];

        for ii = 1:numel(rankedIdx)

            f = rankedIdx(ii);

            if isempty(selectedIdx)
                selectedIdx = f;
            else
                keep = true;

                for jj = 1:numel(selectedIdx)
                    g = selectedIdx(jj);

                    d1 = deltaMat(:,f);
                    d2 = deltaMat(:,g);

                    good = isfinite(d1) & isfinite(d2);

                    if sum(good) >= 8
                        r = corr(d1(good), d2(good), 'Type','Spearman');

                        if abs(r) >= corrThreshold
                            keep = false;
                            break;
                        end
                    end
                end

                if keep
                    selectedIdx(end+1) = f; %#ok<AGROW>
                end
            end

            if numel(selectedIdx) >= fsK
                break;
            end
        end

    else
        selectedIdx = rankedIdx(1:min(fsK, numel(rankedIdx)));
    end

    selectedIdx = selectedIdx(:)';
end

%% =====================================================================
%% ================= MATCHED DONOR / CONTROL RAW CNN FUNCTIONS ===========
%% =====================================================================

function out = runMatchedRawCNNParticipantEvaluation(RawX, Y, WindowInfo, CFG)

    rng(CFG.randomSeed);

    WindowInfo.group = lower(string(WindowInfo.group));
    WindowInfo.participantKey = string(WindowInfo.participantKey);

    donorPids   = unique(WindowInfo.participantKey(WindowInfo.group == "donor"));
    controlPids = unique(WindowInfo.participantKey(WindowInfo.group == "control"));

    donorPids = donorPids(~ismissing(donorPids));
    controlPids = controlPids(~ismissing(controlPids));

    fprintf('Raw-CNN matched setup participants after QC: donors=%d, controls=%d\n', ...
        numel(donorPids), numel(controlPids));

    if isempty(donorPids) || isempty(controlPids)
        error('Need both donor and control participants for raw CNN matched evaluation.');
    end

    allResults = table();
    allPred = table();
    allScores = table();

    for rep = 1:CFG.nMatchedRepeats

        switch lower(string(CFG.trainGroup))
            case "donor"
                rp = randperm(numel(donorPids));
                nTrain = max(2, floor(CFG.rawCNNTrainFrac * numel(donorPids)));
                nTrain = min(nTrain, numel(donorPids)-1);

                trainPids = donorPids(rp(1:nTrain));
                testSamePids = donorPids(rp(nTrain+1:end));

                switch lower(string(CFG.rawCNNOppositeTestMode))
                    case "all"
                        testOppPids = controlPids;
                    case "matched"
                        nOpp = min(numel(controlPids), numel(testSamePids));
                        rpC = randperm(numel(controlPids));
                        testOppPids = controlPids(rpC(1:nOpp));
                    otherwise
                        error('Unknown CFG.rawCNNOppositeTestMode: %s', CFG.rawCNNOppositeTestMode);
                end

                sameLabel = "heldout_donors";
                oppLabel  = "controls";

            case "control"
                rp = randperm(numel(controlPids));
                nTrain = max(2, floor(CFG.rawCNNTrainFrac * numel(controlPids)));
                nTrain = min(nTrain, numel(controlPids)-1);

                trainPids = controlPids(rp(1:nTrain));
                testSamePids = controlPids(rp(nTrain+1:end));

                switch lower(string(CFG.rawCNNOppositeTestMode))
                    case "all"
                        testOppPids = donorPids;
                    case "matched"
                        nOpp = min(numel(donorPids), numel(testSamePids));
                        rpD = randperm(numel(donorPids));
                        testOppPids = donorPids(rpD(1:nOpp));
                    otherwise
                        error('Unknown CFG.rawCNNOppositeTestMode: %s', CFG.rawCNNOppositeTestMode);
                end

                sameLabel = "heldout_controls";
                oppLabel  = "donors";

            otherwise
                error('Unknown CFG.trainGroup: %s', CFG.trainGroup);
        end

        modelPack = fitRawCNNModelOnParticipants(RawX, Y, WindowInfo, trainPids, CFG);

        [predSame, scoreSame] = applyRawCNNModelToParticipants(RawX, Y, WindowInfo, testSamePids, modelPack, CFG);
        [predOpp,  scoreOpp]  = applyRawCNNModelToParticipants(RawX, Y, WindowInfo, testOppPids,  modelPack, CFG);

        resSame = summarizePredictionSet(predSame, rep, string(CFG.trainGroup), sameLabel, ...
            numel(trainPids), numel(testSamePids), modelPack.nInputChannels);
        resOpp  = summarizePredictionSet(predOpp, rep, string(CFG.trainGroup), oppLabel, ...
            numel(trainPids), numel(testOppPids), modelPack.nInputChannels);

        allResults = [allResults; resSame; resOpp]; %#ok<AGROW>

        predSame.repeat = repmat(rep, height(predSame), 1);
        predSame.trainGroup = repmat(string(CFG.trainGroup), height(predSame), 1);
        predSame.testSet = repmat(sameLabel, height(predSame), 1);

        predOpp.repeat = repmat(rep, height(predOpp), 1);
        predOpp.trainGroup = repmat(string(CFG.trainGroup), height(predOpp), 1);
        predOpp.testSet = repmat(oppLabel, height(predOpp), 1);

        allPred = [allPred; predSame; predOpp]; %#ok<AGROW>

        scoreSame.repeat = repmat(rep, height(scoreSame), 1);
        scoreSame.trainGroup = repmat(string(CFG.trainGroup), height(scoreSame), 1);
        scoreSame.testSet = repmat(sameLabel, height(scoreSame), 1);

        scoreOpp.repeat = repmat(rep, height(scoreOpp), 1);
        scoreOpp.trainGroup = repmat(string(CFG.trainGroup), height(scoreOpp), 1);
        scoreOpp.testSet = repmat(oppLabel, height(scoreOpp), 1);

        allScores = [allScores; scoreSame; scoreOpp]; %#ok<AGROW>

        if rep == 1 || mod(rep,10)==0 || rep == CFG.nMatchedRepeats
            fprintf('RawCNN repeat %3d/%3d | train=%d | %s AUC=%.3f balAcc=%.3f | %s AUC=%.3f balAcc=%.3f\n', ...
                rep, CFG.nMatchedRepeats, numel(trainPids), ...
                sameLabel, resSame.AUC, resSame.BalancedAccuracy, ...
                oppLabel, resOpp.AUC, resOpp.BalancedAccuracy);
        end
    end

    summary = summarizeMatchedResults(allResults);

    out = struct();
    out.results = allResults;
    out.predictions = allPred;
    out.participantScores = allScores;
    out.summary = summary;
end

function modelPack = fitRawCNNModelOnParticipants(RawX, Y, WindowInfo, trainPids, CFG)

    trainIdx = ismember(string(WindowInfo.participantKey), string(trainPids));

    XTrain = RawX(trainIdx);
    YTrain = Y(trainIdx);

    if numel(unique(YTrain)) < 2
        error('Raw CNN training set has fewer than 2 classes.');
    end

    [XTrain, YTrain] = balanceClasses(XTrain, YTrain);

    numInputChannels = size(XTrain{1}, 1);
    minSeqLen = min(cellfun(@(x) size(x,2), XTrain));

    fprintf('RawCNN train: participants=%d | windows=%d | channels=%d | minSeqLen=%d\n', ...
        numel(trainPids), numel(XTrain), numInputChannels, minSeqLen);

    layers = makeRawCNNLayers(numInputChannels, minSeqLen);

    options = trainingOptions('adam', ...
        'MaxEpochs', CFG.maxEpochs, ...
        'MiniBatchSize', CFG.miniBatchSize, ...
        'InitialLearnRate', CFG.initialLearnRate, ...
        'Shuffle','every-epoch', ...
        'Verbose',false, ...
        'Plots','none');

    net = trainNetwork(XTrain, YTrain, layers, options);

    cats = categories(YTrain);
    lateCol = find(strcmp(cats, 'late'));
    if isempty(lateCol)
        error('Could not find late category in CNN training labels.');
    end

    modelPack = struct();
    modelPack.net = net;
    modelPack.lateCol = lateCol;
    modelPack.trainCategories = cats;
    modelPack.nInputChannels = numInputChannels;
    modelPack.minSeqLen = minSeqLen;
end

function [Pred, ParticipantScores] = applyRawCNNModelToParticipants(RawX, Y, WindowInfo, testPids, modelPack, CFG)

    testIdx = ismember(string(WindowInfo.participantKey), string(testPids));

    if ~any(testIdx)
        Pred = table();
        ParticipantScores = table();
        return;
    end

    XTest = RawX(testIdx);
    YTest = Y(testIdx);
    Wte = WindowInfo(testIdx,:);

    scores = predict(modelPack.net, XTest, 'MiniBatchSize', CFG.miniBatchSize);
    yPredCat = classify(modelPack.net, XTest, 'MiniBatchSize', CFG.miniBatchSize);

    lateScore = scores(:, modelPack.lateCol);
    trueLate = double(YTest == "late");
    predictedLate = double(string(yPredCat) == "late");

    Pred = table();
    Pred.participant = Wte.participantKey;
    Pred.participantKey = Wte.participantKey;
    Pred.group = Wte.group;
    Pred.donationFractionMid = Wte.donationFractionMid;
    Pred.trueLate = trueLate;
    Pred.predictedLate = predictedLate;
    Pred.lateScore = lateScore;

    ParticipantScores = computeParticipantScores(Pred);
end
