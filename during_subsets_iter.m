
close all
clc
clear
rng(1500)
Fs=25; 
noiseThreshold=0.05;
to_plot=0;

load('Holter_timings_controls.mat');
subjData1=subjData;
[subjData1(:).Group]=deal('control');
load('Holter_timings.mat');
[subjData(:).Group]=deal('Donors');
subjData=rmfield(subjData,{'P_Donation_Amount'});

subjDataAll=[subjData,subjData1];
subjDataAll(isnan([subjDataAll.Weight]))=[];
 subjDataAll([subjDataAll.Weight]>90)=[];

  subjDataAll(strcmpi({subjDataAll.code},'032'))=[];

BM_features_Names = {...
    'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
    'PercentOfBreathsWithInhalePause'};

noiseThreshold=0.05;
device_name='Holter';

%featureFcn = @(seg,Fs) nc_feats(seg, Fs, noiseThreshold, device_name);% 
featureFcn = @(seg,Fs) breathmetrics_feats(seg, Fs);


N=size(subjDataAll,2);
norm=1;

for i = 1:N
            [before{i}, after{i},during{i},NCbefore{i}, NCafter{i},NCdonation{i}]=extract_timings_needle_walk_in_chair2(i,norm,10,subjDataAll);
      %      [before{i}, after{i},during{i},NCbefore{i}, NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm,5,subjDataAll);

end

 nWindows  = 10;
winFrac   = 0.1;     % each window spans 20% of that participant's phase
minWinSec = 90;      % don’t go below 45 s per window (for BM stability)
noiseThreshold = 0.05;         % adjust as needed
device_name    = 'Holter';     % your NC function expects 'Holter' or 'Mustache'



[winTableN, winIdxN] = window_phases_make_bins_fixedN( ...
    before, during, after, Fs, ...
    nWindows, 'frac', winFrac, minWinSec, featureFcn);


% [winTableN, winIdxN] = window_phases_make_bins_fixedN( ...
%     NCbefore, NCdonation, NCafter, Fs, ...
%     nWindows, 'frac', winFrac, minWinSec, featureFcn);

G=sort(repmat({subjDataAll.Group},1,nWindows*3))';
winTableN.Group = categorical(G);

meta = {'participant','Group','phase','win_index','t_start_s','t_end_s','phase_frac_start','phase_frac_end'};
%featList = setdiff(string(winTableN.Properties.VariableNames), meta);
featList={'AverageInhaleVolume', 'AveragePeakExpiratoryFlow','AveragePeakInspiratoryFlow','AverageTidalVolume','MinuteVentilation', 'PercentOfBreathsWithExhalePause',...
    'DutyCycleOfExhalePause','AverageExhalePauseDuration','BreathingRate'};
 % 'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
 %    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
 %    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
 %    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
 %    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
 %    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause',
% 

nIter    = 10000;   % number of random donor subsets
subsetN  = 30;    % number of participants per subset

donorLabel   = "Donors";    % adjust to your actual labels
controlLabel = "control";  % adjust to your actual labels

for f = 1:numel(featList)
    feat = featList{f};
 %   if ~isnumeric(winTableN.(feat)), continue; end

    % ----- DURING only, split donors vs controls -----
    isDuring  = winTableN.phase=="during";
    isDonor   = winTableN.Group==donorLabel;
    isControl = winTableN.Group==controlLabel;

    Tdon  = winTableN(isDuring & isDonor, :);
    Tctrl = winTableN(isDuring & isControl, :);

    if isempty(Tdon) || isempty(Tctrl)
        fprintf('Skipping %s: missing donor/control during data.\n', char(feat));
        continue;
    end

    Tdon.participant  = categorical(Tdon.participant);
    Tdon.win_index    = double(Tdon.win_index);
    Tctrl.participant = categorical(Tctrl.participant);
    Tctrl.win_index   = double(Tctrl.win_index);

    % ---------- 1) Donor subsets: 20 participants each ----------
    subsDon  = categories(Tdon.participant);
    nSubsDon = numel(subsDon);

    if subsetN > nSubsDon
        warning('Not enough donors with DURING data for %s', char(feat));
        continue;
    end

    slope_sub = nan(nIter,1);    % donor subset slopes

    for it = 1:nIter
        % pick 20 unique donors (order within each donor’s bins is preserved)
        idxSub    = randperm(nSubsDon, subsetN);
        chosenIDs = subsDon(idxSub);

        maskDon = isDuring & isDonor & ...
                  ismember(categorical(winTableN.participant), chosenIDs);
        Tsub = winTableN(maskDon, :);

        if isempty(Tsub), continue; end

        Tsub.participant = categorical(Tsub.participant);
        Tsub.win_index   = double(Tsub.win_index);

        outB   = mixed_effects_during(Tsub, feat, 'time','linear');
        fixedB = outB.fixed;

        b1_idx = strcmp(fixedB.Name,'win_c');
        if any(b1_idx)
            slope_sub(it) = fixedB.Estimate(b1_idx);
        end
    end

    % remove any NaNs from failed fits
    slope_sub = slope_sub(isfinite(slope_sub));

    if isempty(slope_sub)
        fprintf('No valid donor subset slopes for %s\n', char(feat));
        continue;
    end

    % ---------- 2) Control slope (all controls, n=20) ----------
    outCtrl   = mixed_effects_during(Tctrl, feat, 'time','linear');
    fixedCtrl = outCtrl.fixed;
    b1_idx    = strcmp(fixedCtrl.Name,'win_c');

    if ~any(b1_idx)
        fprintf('No win_c term in control model for %s\n', char(feat));
        continue;
    end

    slope_ctrl = fixedCtrl.Estimate(b1_idx);

    % ---------- 3) Compare control slope to donor-subset distribution ----------

    % Two-sided “donor-based” p-value:
    % How extreme is the control slope relative to donor subset slopes?
    p_ctrl_vs_don(f) = mean(abs(slope_sub) >= abs(slope_ctrl));

    fprintf('\nFeature: %s\n', char(feat));
    fprintf('  Control slope: %.4g\n', slope_ctrl);
    fprintf('  Donor subset slopes: median = %.4g, IQR = [%.4g, %.4g]\n', ...
        median(slope_sub), prctile(slope_sub,25), prctile(slope_sub,75));
    fprintf('  p_ctrl_vs_don (two-sided) = %.4f\n', p_ctrl_vs_don(f));

    % ---------- 4) Plot histogram + control line ----------
    figure('Color','w','Name',['Donor subsets vs control — %s',feat]);
    histogram(slope_sub, 'NumBins',15, 'Normalization','pdf');
    hold on;
    xline(slope_ctrl, 'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);

    xlabel('Slope during donation');
    ylabel('Density');
    % title(sprintf('%s: donor 30-subset slopes vs control (p=%.3g)', ...
    %     strrep(char(feat),'_',' '), p_ctrl_vs_don(f)), 'Interpreter','none');
    legend({'Donor 30-subset slopes','Control slope'}, 'Location','best');
end



q = mafdr(p_ctrl_vs_don,'BHFDR',true);   % Benjamini–Hochberg FDR-adjusted p-values

% Or threshold at alpha = 0.05:
alpha  = 0.05;
sigFDR = q < alpha; 
find(sigFDR)

function [winTable, winIdx] = window_phases_make_bins_fixedN( ...
        beforeCell, duringCell, afterCell, Fs, nWindows, mode, winParam, minWinSec, featureFcn)

    if nargin < 9 || isempty(featureFcn)
        featureFcn = @(x,Fs) struct('mean',mean(x,'omitnan'),'std',std(x,'omitnan'));
    end
    assert(isscalar(nWindows) && nWindows>=1 && floor(nWindows)==nWindows, 'nWindows must be a positive integer.');
    assert(ismember(mode, {'frac','abs'}), 'mode must be ''frac'' or ''abs''.');

    N = numel(beforeCell);
    phases     = {'before','during','after'};
    phaseCells = {beforeCell, duringCell, afterCell};

    minWinSamp = max(1, round(minWinSec * Fs));

    allRows = {};
    winIdx.before = cell(N,1);  winIdx.during = cell(N,1);  winIdx.after = cell(N,1);

    for p = 1:numel(phases)
        phName = phases{p};
        series = phaseCells{p};

        for i = 1:N
            x = series{i};

            % ==== CHANGED: accept vector OR 2-col matrix ====
            if isempty(x) || ~(ismatrix(x) && size(x,1)>=2 && size(x,2)>=1 && size(x,2)<=2)
                winIdx.(phName){i} = zeros(0,2);  continue;
            end

            % ==== CHANGED: finite rows require all columns finite for NC ====
            finiteMask = all(isfinite(x), 2);
            if ~any(finiteMask), winIdx.(phName){i} = zeros(0,2); continue; end
            firstFinite = find(finiteMask,1,'first');
            lastFinite  = find(finiteMask,1,'last');
            x = x(firstFinite:lastFinite, :);
            L = size(x,1);
            if L < 2, winIdx.(phName){i} = zeros(0,2); continue; end

            % === Window length (samples) ===
            switch mode
                case 'frac'
                    winFrac = winParam;
                    assert(winFrac>0 && winFrac<=1, 'winFrac must be in (0,1].');
                    winSamp = max(minWinSamp, round(winFrac * L));
                case 'abs'
                    winSec  = winParam;
                    assert(winSec>0, 'winSec must be >0.');
                    winSamp = max(minWinSamp, round(winSec * Fs));
            end
            if winSamp > L
                winIdx.(phName){i} = zeros(0,2);  continue;
            end

            if nWindows == 1
                starts = round((L - winSamp)/2) + 1;
            else
                starts = round(linspace(1, L - winSamp + 1, nWindows));
            end
            ends = starts + winSamp - 1;
            starts = max(1, min(starts, L - winSamp + 1));
            ends   = min(ends, L);

            winIdx.(phName){i} = [starts(:) ends(:)];

            % Build rows
            nW = numel(starts);
            theseRows = cell(nW,1);
            for w = 1:nW
                seg = x(starts(w):ends(w), :);        % ==== CHANGED: keep as matrix
                feats = featureFcn(seg, Fs);          % wrapper returns scalars

                row.participant       = i;
                row.phase             = categorical({phName}, phases);
                row.win_index         = w;
                row.t_start_s         = (starts(w)-1)/Fs;
                row.t_end_s           = (ends(w)-1)/Fs;
                row.phase_frac_start  = (starts(w)-1) / (L-1);
                row.phase_frac_end    = (ends(w)-1)   / (L-1);

                fns = fieldnames(feats);
                for ff = 1:numel(fns)
                    row.(fns{ff}) = feats.(fns{ff});
                end
                theseRows{w} = struct2table(row);
            end

            if ~isempty(theseRows)
                allRows{end+1} = vertcat(theseRows{:}); %#ok<AGROW>
            end
        end
    end

    if isempty(allRows)
        winTable = table();
    else
        winTable = sortrows(vertcat(allRows{:}), {'participant','phase','win_index'});
    end
end

function feats = nc_feats(seg, Fs, noiseThreshold, device_name)
% seg: [samples x 2] (Left, Right) airflow
% Returns scalar NC metrics for the window.

    % safety & formatting
    if isempty(seg) || size(seg,2) ~= 2
        feats = struct('MeanLateralityIndex2',NaN,'MeanLateralityIndex3',NaN, ...
                       'MeanAmplitudeLI2',NaN,'MeanAmplitudeLI3',NaN, ...
                       'stdLI2',NaN,'stdLI3',NaN,'Nostril_Corr_R2',NaN);
        return;
    end

    % Call your NC function (it handles Holter/Mustache column order)
    [LI, meas] = NasalCycleParameters_short(seg, Fs, noiseThreshold, device_name);

    % Map to scalars (one value per window/feature)
    feats = struct();
    feats.MeanLateralityIndex2 = safeget(meas,'MeanLateralityIndex2',NaN);
    feats.MeanLateralityIndex3 = safeget(meas,'MeanLateralityIndex3',NaN);
    feats.MeanAmplitudeLI2     = safeget(meas,'MeanAmplitudeLI2',NaN);
    feats.MeanAmplitudeLI3     = safeget(meas,'MeanAmplitudeLI3',NaN);
    feats.stdLI2               = safeget(meas,'stdLI2',NaN);
    feats.stdLI3               = safeget(meas,'stdLI3',NaN);
    feats.Nostril_Corr_R2      = safeget(meas,'Nostril_Corr_R2',NaN);

    % (Optional) If you ALSO want a robust location of the LI series in-window:
    if isfield(LI,'two') && ~isempty(LI.two)
        feats.MedianLI2 = median(LI.two,'omitnan');
    else
        feats.MedianLI2 = NaN;
    end
    if isfield(LI,'three') && ~isempty(LI.three)
        feats.MedianLI3 = median(LI.three,'omitnan');
    else
        feats.MedianLI3 = NaN;
    end
end

function out = mixed_effects_during(winTable, featureName, varargin)
% Mixed-effects analysis of trend within the "during" phase.
% time: 'linear' (default), 'cat', or 'spline' (cubic polynomial)
%
% Example:
%   out = mixed_effects_during(winTableN, "BreathingRate", 'time','linear');

    % ---- args ----
    if ischar(featureName); featureName = string(featureName); end
    validateattributes(featureName, {'string','char'}, {'nonempty'});

    p = inputParser;
    addParameter(p,'time','linear', @(s) any(strcmpi(s,{'linear','cat','spline'})));
    addParameter(p,'centerTime',true, @(x)islogical(x)||ismember(x,[0 1]));
    parse(p,varargin{:});
    modeTime   = lower(p.Results.time);
    centerTime = p.Results.centerTime;

    % ---- subset DURING ----
    T = winTable(winTable.phase=="during", :);
    if isempty(T)
        error('No rows with phase=="during" in winTable.');
    end
    if ~ismember(featureName, string(T.Properties.VariableNames))
        error('Feature "%s" not found.', featureName);
    end
    if ~isnumeric(T.(featureName))
        error('Feature "%s" must be numeric.', featureName);
    end

    % keep only needed columns (use STRING array!)
    T = T(:, ["participant","win_index", featureName]);
    T.participant = categorical(T.participant);
    T.win_index  = double(T.win_index);

    % center time (helps)
    if centerTime
        muIdx = mean(T.win_index,'omitnan');
        T.win_c = T.win_index - muIdx;
    else
        T.win_c = T.win_index;
    end

    % ---- choose formula ----
    fname = char(featureName);  % for sprintf
    switch modeTime
        case 'linear'
            % random intercept + random slope by participant
            formula = sprintf('%s ~ 1 + win_c + (1 + win_c | participant)', fname);

        case 'cat'
            T.win_cat = categorical(T.win_index);
            % random intercept per participant
            formula = sprintf('%s ~ 1 + win_cat + (1 | participant)', fname);

        case 'spline'
            % cubic polynomial in centered time
            T.win_c2 = T.win_c.^2;
            T.win_c3 = T.win_c.^3;
            formula = sprintf('%s ~ 1 + win_c + win_c2 + win_c3 + (1 + win_c | participant)', fname);
    end

    % ---- fit LME ----
    lme = fitlme(T, formula); %, 'FitMethod','REML', 'DFMethod','Kenward-Roger'

    % ---- outputs ----
    out = struct();
    out.model      = lme;
    out.anova      = anova(lme);   % fixed effect tests %,'DFMethod','Kenward-Roger'
    out.fixed      = lme.Coefficients;                        % table already
    out.designInfo = struct('mode',modeTime,'centered',centerTime,'formula',formula);

    % console summary
    % fprintf('\n=== Mixed Effects (during) — %s | time=%s ===\n', fname, modeTime);
    % disp(out.anova);
    % disp(out.fixed);
end

function val = safeget(s, field, defaultVal)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

function feats = breathmetrics_feats(x, Fs)
% BREATHMETRICS_FEATS
% Compute per-window BreathMetrics secondary features for an airflow signal x.
% x: vector (airflow), Fs: sampling rate (Hz).

    x = x(:);  % column

    % Build BreathMetrics object (human airflow mode)
    bmObj = breathmetrics(x, Fs, 'humanAirflow');

    % Compute all features; sliding=1, plotting=0 (matches your snippet)
    bmObj.estimateAllFeatures(0, 'sliding', 1, 0);

    % Desired secondary-feature names (same order you provided)
    variableNames = { ...
        'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', ...
        'AverageInhaleDuration', 'AverageInhalePauseDuration', 'AverageInhaleVolume', ...
        'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', 'AveragePeakInspiratoryFlow', ...
        'AverageTidalVolume', 'BreathingRate', ...
        'CoefficientOfVariationOfBreathVolumes', 'CoefficientOfVariationOfBreathingRate', ...
        'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
        'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', ...
        'DutyCycleOfExhale', 'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', ...
        'MinuteVentilation', 'PercentOfBreathsWithExhalePause', 'PercentOfBreathsWithInhalePause' ...
    };

    % Pull values from BreathMetrics (secondary features)
    vals = bmObj.secondaryFeatures.values;

    % Some BreathMetrics versions return a cell array, others a numeric row vector.
    if iscell(vals)
        vals = cellfun(@(v) double(v), vals);
    else
        vals = double(vals);
    end

    % Assign to struct fields (missing entries become NaN if needed)
    feats = struct();
    for ii = 1:numel(variableNames)
        if ii <= numel(vals)
            feats.(variableNames{ii}) = vals(ii);
        else
            feats.(variableNames{ii}) = NaN;
        end
    end
end