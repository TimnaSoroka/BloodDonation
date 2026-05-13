%% ========================================================================
% STATS FOR MatchedFeature_Results_*.csv IN LONG FORMAT
%
% Your results file has one row per repeat x testSet, e.g.:
%   repeat | trainGroup | testSet         | Accuracy | BalancedAccuracy | AUC
%      1   | donor      | heldout_donors  | ...      | ...              | ...
%      1   | donor      | controls        | ...      | ...              | ...
%
% This script compares matched test sets within each repeat:
%   referenceTestSet - comparisonTestSet
%
% For donor-trained model:
%   referenceTestSet  = "heldout_donors"
%   comparisonTestSet = "controls"
%% ========================================================================

clear; clc;

%% ========================= USER SETTINGS ================================

cfg = struct();

cfg.resultsCsv = '/Users/timnas/Documents/BloodDonation/rawSubjData_earlyLateDonation_matched_donor_control/sum_norm-afterCut_FS-0_K20_donorTrain_100rep/MatchedRawCNN_Results_100.csv';

cfg.outDir     = fileparts(cfg.resultsCsv);

% Direction of the effect you want to test:
% donor-trained specificity: held-out donors should perform better than controls.
cfg.referenceTestSet  = "heldout_donors";
cfg.comparisonTestSet = "controls";

% Metrics to summarize and compare.
cfg.metrics = [
    "AUC"
    "BalancedAccuracy"
    "Accuracy"
    "SensitivityLate"
    "SpecificityEarly"
];

% NOTE: repeats reuse participants, so this is a split-stability p-value,
% not a subject-level inferential p-value.
cfg.makePlots = true;

%% ========================= LOAD =========================================

R = readtable(cfg.resultsCsv, 'TextType','string', 'VariableNamingRule','preserve');

requiredVars = ["repeat","testSet", cfg.metrics(:)'];
missingVars = setdiff(requiredVars, string(R.Properties.VariableNames));
if ~isempty(missingVars)
    error('Missing required columns: %s', strjoin(missingVars, ', '));
end

R.testSet = string(R.testSet);

fprintf('\nLoaded: %s\n', cfg.resultsCsv);
fprintf('Rows: %d\n', height(R));
fprintf('Unique repeats: %d\n', numel(unique(R.repeat))); 
fprintf('Test sets found:\n');
disp(unique(R.testSet));

%% ========================= BUILD WIDE PER-REPEAT TABLE ===================

repeats = unique(R.repeat);
repeats = repeats(:);

Wide = table();
Wide.repeat = repeats;

for m = 1:numel(cfg.metrics)
    metric = cfg.metrics(m);

    refVals = nan(numel(repeats),1);
    cmpVals = nan(numel(repeats),1);

    for i = 1:numel(repeats)
        rep = repeats(i);

        idxRef = R.repeat == rep & R.testSet == cfg.referenceTestSet;
        idxCmp = R.repeat == rep & R.testSet == cfg.comparisonTestSet;

        if sum(idxRef) == 1
            refVals(i) = R.(metric)(idxRef);
        elseif sum(idxRef) > 1
            refVals(i) = mean(R.(metric)(idxRef), 'omitnan');
        end

        if sum(idxCmp) == 1
            cmpVals(i) = R.(metric)(idxCmp);
        elseif sum(idxCmp) > 1
            cmpVals(i) = mean(R.(metric)(idxCmp), 'omitnan');
        end
    end

    refName = matlab.lang.makeValidName(cfg.referenceTestSet + "_" + metric);
    cmpName = matlab.lang.makeValidName(cfg.comparisonTestSet + "_" + metric);
    delName = matlab.lang.makeValidName("Delta_" + cfg.referenceTestSet + "_minus_" + cfg.comparisonTestSet + "_" + metric);

    Wide.(refName) = refVals;
    Wide.(cmpName) = cmpVals;
    Wide.(delName) = refVals - cmpVals;
end

wideOut = fullfile(cfg.outDir, 'MatchedFeature_Results_WIDE_byRepeat.csv');
writetable(Wide, wideOut);
fprintf('\nSaved wide per-repeat table:\n%s\n', wideOut);

%% ========================= SUMMARY STATS =================================

Summary = table();
rowCounter = 0;

for m = 1:numel(cfg.metrics)
    metric = cfg.metrics(m);

    refName = matlab.lang.makeValidName(cfg.referenceTestSet + "_" + metric);
    cmpName = matlab.lang.makeValidName(cfg.comparisonTestSet + "_" + metric);
    delName = matlab.lang.makeValidName("Delta_" + cfg.referenceTestSet + "_minus_" + cfg.comparisonTestSet + "_" + metric);

    names = [refName; cmpName; delName];
    labels = [
        cfg.referenceTestSet + "_" + metric
        cfg.comparisonTestSet + "_" + metric
        cfg.referenceTestSet + "_minus_" + cfg.comparisonTestSet + "_" + metric
    ];

    for j = 1:numel(names)
        x = Wide.(names(j));
        x = x(isfinite(x));

        rowCounter = rowCounter + 1;
        Summary.metric(rowCounter,1) = labels(j);
        Summary.N(rowCounter,1) = numel(x);
        Summary.mean(rowCounter,1) = mean(x, 'omitnan');
        Summary.sd(rowCounter,1) = std(x, 'omitnan');
        Summary.median(rowCounter,1) = median(x, 'omitnan');
        Summary.pct025(rowCounter,1) = prctile(x, 2.5);
        Summary.pct975(rowCounter,1) = prctile(x, 97.5);

        if contains(labels(j), "_minus_")
            Summary.proportion_gt_0(rowCounter,1) = mean(x > 0, 'omitnan');
            Summary.empirical_p_one_sided_gt_0(rowCounter,1) = (sum(x <= 0) + 1) / (numel(x) + 1);
            Summary.empirical_p_two_sided_sign(rowCounter,1) = 2 * min( ...
                (sum(x <= 0) + 1) / (numel(x) + 1), ...
                (sum(x >= 0) + 1) / (numel(x) + 1));
            Summary.empirical_p_two_sided_sign(rowCounter,1) = min(Summary.empirical_p_two_sided_sign(rowCounter,1), 1);
        else
            Summary.proportion_gt_0(rowCounter,1) = NaN;
            Summary.empirical_p_one_sided_gt_0(rowCounter,1) = NaN;
            Summary.empirical_p_two_sided_sign(rowCounter,1) = NaN;
        end
    end
end

summaryOut = fullfile(cfg.outDir, 'MatchedFeature_Stats_longFormat.csv');
writetable(Summary, summaryOut);

fprintf('\nSummary stats:\n');
disp(Summary);
fprintf('\nSaved summary stats:\n%s\n', summaryOut);

%% ========================= MAIN RESULT PRINT =============================

aucDeltaName = matlab.lang.makeValidName("Delta_" + cfg.referenceTestSet + "_minus_" + cfg.comparisonTestSet + "_AUC");
if ismember(aucDeltaName, string(Wide.Properties.VariableNames))
    dAUC = Wide.(aucDeltaName);
    dAUC = dAUC(isfinite(dAUC));

    fprintf('\nMAIN AUC RESULT (%s - %s):\n', cfg.referenceTestSet, cfg.comparisonTestSet);
    fprintf('N repeats = %d\n', numel(dAUC));
    fprintf('Mean delta AUC = %.4f\n', mean(dAUC, 'omitnan'));
    fprintf('Median delta AUC = %.4f\n', median(dAUC, 'omitnan'));
    fprintf('95%% split interval = [%.4f, %.4f]\n', prctile(dAUC,2.5), prctile(dAUC,97.5));
    fprintf('Proportion delta > 0 = %.3f\n', mean(dAUC > 0, 'omitnan'));
    fprintf('Empirical one-sided p(delta > 0) = %.4f\n', (sum(dAUC <= 0) + 1) / (numel(dAUC) + 1));
end

%% ========================= PLOTS ========================================

if cfg.makePlots

    plotDir = fullfile(cfg.outDir, 'stats_plots');
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end

    for m = 1:numel(cfg.metrics)
        metric = cfg.metrics(m);
        refName = matlab.lang.makeValidName(cfg.referenceTestSet + "_" + metric);
        cmpName = matlab.lang.makeValidName(cfg.comparisonTestSet + "_" + metric);
        delName = matlab.lang.makeValidName("Delta_" + cfg.referenceTestSet + "_minus_" + cfg.comparisonTestSet + "_" + metric);

        if ~all(ismember([refName cmpName delName], string(Wide.Properties.VariableNames)))
            continue;
        end

        xRef = Wide.(refName);
        xCmp = Wide.(cmpName);
        xDel = Wide.(delName);

        % Paired box/chart-like view
        fig = figure('Color','w','Name',char(metric));
        hold on;
        g = categorical([repmat(cfg.referenceTestSet, numel(xRef), 1); repmat(cfg.comparisonTestSet, numel(xCmp), 1)]);
        y = [xRef; xCmp];
        boxchart(g, y);
        ylabel(metric, 'Interpreter','none');
        title(metric + ": " + cfg.referenceTestSet + " vs " + cfg.comparisonTestSet, 'Interpreter','none');
        grid on;
        saveas(fig, fullfile(plotDir, "box_" + metric + ".png"));

        % Delta histogram
        fig = figure('Color','w','Name',char("Delta_" + metric));
        histogram(xDel, 30);
        hold on;
        xline(0, '--k', 'LineWidth', 1.5);
        xline(mean(xDel,'omitnan'), 'LineWidth', 2);
        xlabel(cfg.referenceTestSet + " - " + cfg.comparisonTestSet + " " + metric, 'Interpreter','none');
        ylabel('Number of repeats');
        title("Delta " + metric, 'Interpreter','none');
        grid on;
        saveas(fig, fullfile(plotDir, "delta_hist_" + metric + ".png"));
    end
end

fprintf('\nDone.\n');
