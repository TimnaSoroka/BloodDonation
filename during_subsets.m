%% during subset iterations

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

% idxBelow500 = cellfun(@(x) isnumeric(x) && ~isempty(x) && all(x < 500), {subjData.Donation_Amount});
% subjData(idxBelow500)=[];
subjDataAll=[subjData,subjData1];
% subjDataAll(isnan([subjDataAll.Weight]))=[];
%  subjDataAll([subjDataAll.Weight]>90)=[];
% 
%   subjDataAll(strcmpi({subjDataAll.code},'032'))=[];

rmCodes  = {'001','002','003','013','014','015','016','017','018','019','020','021','022','023','024','025',...
    '026','027','028','029','030','031','096','097','098','099','100','101','102','103','104','105','106','107'};
allCodes = {subjDataAll.code};

%001-003	בא"ח חיפה	02/02/2025	ניידת
%013-031	בא"ח חיפה	11-13/02/2025	ניידת
%096-102	קרייה	16/07/2025	ניידת
%103-107	תלה"ש	14/08/2025	ניידת
 subjDataAll(ismember({subjDataAll.code}, rmCodes)) = [];




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

% featureFcn = @(seg,Fs) nc_feats(seg, Fs, noiseThreshold, device_name);% 
 featureFcn = @(seg,Fs) breathmetrics_feats(seg, Fs);


N=size(subjDataAll,2);
norm=1;

for i = 1:N
            [before{i}, after{i},during{i},NCbefore{i}, NCafter{i},NCdonation{i}]=extract_timings_needle_walk_in_chair2(i,norm,10,subjDataAll);
      %      [before{i}, after{i},during{i},NCbefore{i}, NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm,5,subjDataAll);

end

 nWindows  = 10;
winFrac   = 1/10;     % each window spans 20% of that participant's phase
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
featList = setdiff(string(winTableN.Properties.VariableNames), meta);
% featList={'AverageInhaleVolume', 'AveragePeakExpiratoryFlow','AveragePeakInspiratoryFlow','AverageTidalVolume','MinuteVentilation', 'PercentOfBreathsWithExhalePause',...
%     'DutyCycleOfExhalePause','AverageExhalePauseDuration','BreathingRate'};
 % 'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
 %    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
 %    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
 %    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
 %    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
% % %  %    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause',

rng(1);

phaseName   = "during";
groupDonor  = "Donors";
groupCtrl   = "control";

 subsetSize  = 30;      % donors per subset
nIter       = 10000;    % number of subsets
 alpha       = 0.05;    % gate for donor trend

for f=1:size(featList,2)
feat = featList(f);   % <-- change to any of your 24 features
fprintf(feat);
sprintf('\n');

%% --- 1) Donors: group mean per bin + slope ---
Tdon = winTableN(winTableN.Group==groupDonor & winTableN.phase==phaseName, :);
[xDon, muDon, semDon] = group_bin_stats_table(Tdon, feat);

if numel(xDon) < 2
    error('Not enough bins for donors in phase %s.', phaseName);
end

mdlDon   = fitlm(xDon, muDon);
coefDon  = mdlDon.Coefficients;
slopeDon = coefDon.Estimate(2);
pDon     = coefDon.pValue(2);

fprintf('Donors (full group): slope = %.4g, p = %.3g\n', slopeDon, pDon);
ppsDon(f)=pDon;

%% --- 2) Controls: group mean per bin + slope (same method) ---
Tctrl = winTableN(winTableN.Group==groupCtrl & winTableN.phase==phaseName, :);
[xCtrl, muCtrl, semCtrl] = group_bin_stats_table(Tctrl, feat);

if numel(xCtrl) >= 2
    mdlCtrl   = fitlm(xCtrl, muCtrl);
    slopeCtrl = mdlCtrl.Coefficients.Estimate(2);
        pCtrl     = mdlCtrl.Coefficients.pValue(2);
    fprintf('Controls (full group): slope = %.4g ,p = %.2f\n', slopeCtrl,pCtrl);
else
    slopeCtrl = NaN;
    fprintf('Controls: not enough bins for slope.\n');
    pCtrl     = NaN;
end

%% --- 2b) Plot group mean over bins (Donors vs Controls) ---
figure('Color','w', 'Name', sprintf('Group means — %s', feat));
hold on;

% Colors
colDon  = [0 0.4470 0.7410];      % blue
colCtrl = [0.4660 0.6740 0.1880]; % green

%% --- Donors: shaded SEM + line ---
hLines   = [];      % handles for legend
legNames = {};      % corresponding names

%% --- Donors: shaded SEM + line ---
if ~isempty(xDon)
    xd  = xDon(:)';
    mud = muDon(:)';
    sed = semDon(:)';

    % Shaded area
    xFill = [xd, fliplr(xd)];
    yFill = [mud + sed, fliplr(mud - sed)];
    patch(xFill, yFill, colDon, ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');

    % Mean line
    hDon = plot(xd, mud, '-o', ...
        'Color', colDon, ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', colDon);

    hLines(end+1)   = hDon;
    legNames{end+1} = 'Donors';
end

%% --- Controls: shaded SEM + line ---
if ~isempty(xCtrl)
    xc  = xCtrl(:)';
    muc = muCtrl(:)';
    sec = semCtrl(:)';

    xFill = [xc, fliplr(xc)];
    yFill = [muc + sec, fliplr(muc - sec)];
    patch(xFill, yFill, colCtrl, ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');

    hCtrl = plot(xc, muc, '-s', ...
        'Color', colCtrl, ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', colCtrl);

    hLines(end+1)   = hCtrl;
    legNames{end+1} = 'Controls';
end

xlabel('Bin index (during intervention)');
ylabel(strrep(char(feat),'_',' '));
if ~isempty(hLines)
    legend(hLines, legNames, 'Location','best');
end
title(sprintf('Group mean %s across bins (during)', feat));
hold off;


% %% === Donors: shaded SEM ===
if ~isempty(xDon)
    figure('Color','w', 'Name', sprintf('Donors — %s', feat));
    hold on;

    % Shaded area: mean ± SEM
    xFill = [xDon, fliplr(xDon)];
    yFill = [muDon + semDon, fliplr(muDon - semDon)];
    patch(xFill, yFill, colDon, ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');

    % Mean line
    plot(xDon, muDon, '-o', ...
        'Color', colDon, ...
        'LineWidth', 1.8, ...
        'MarkerFaceColor', colDon);

    xlabel('Bin index (during intervention)');
    ylabel(strrep(char(feat),'_',' '));
    title(sprintf('%s — Donors (slope = %.3g, p = %.3g)', ...
        strrep(char(feat),'_',' '), slopeDon, pDon));
    xlim([1, max(xDon)+0.5]);
    hold off;
end
% 
%% === Controls: shaded SEM ===
if ~isempty(xCtrl) && numel(xCtrl) >= 2
    figure('Color','w', 'Name', sprintf('Controls — %s', feat));
    hold on;

    xFill = [xCtrl, fliplr(xCtrl)];
    yFill = [muCtrl + semCtrl, fliplr(muCtrl - semCtrl)];
    patch(xFill, yFill, colCtrl, ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');

    plot(xCtrl, muCtrl, '-o', ...
        'Color', colCtrl, ...
        'LineWidth', 1.8, ...
        'MarkerFaceColor', colCtrl);

    xlabel('Bin index (during intervention)');
    ylabel(strrep(char(feat),'_',' '));
    title(sprintf('%s — Controls (slope = %.3g, p = %.3g)', ...
        strrep(char(feat),'_',' '), slopeCtrl, pCtrl));
    % grid on;
    xlim([1, max(xCtrl)+0.5]);
    hold off;
end


%% --- 3) Gate: only if donor trend is significant, run subsets ---

if pDon >= alpha
    fprintf('→ Donor group trend NOT significant (p = %.3g). Skipping subset analysis.\n', pDon);
    continue;
end

fprintf('→ Donor group trend IS significant (p = %.3g). Running subset analysis...\n', pDon);

%% --- 4) Subsets of donors: each time 20 donors, recompute group curve + slope ---
subsetSlopes = donor_subset_slopes_table(Tdon, feat, subsetSize, nIter);

% permutation-style p: how extreme is control slope vs donor subset slopes?
delta_ctrl       = abs(slopeCtrl - slopeDon);          % distance of control from donor mean
delta_sub        = abs(subsetSlopes - slopeDon);       % distance of each subset from donor mean
pCtrl_vs_subsets = mean(delta_sub >= delta_ctrl);      % permutation p-value


fprintf('Subsets (donors): mean slope = %.4g, SD = %.4g\n', ...
    mean(subsetSlopes), std(subsetSlopes));
fprintf('Control slope vs donor subsets: permutation p ≈ %.3g\n', pCtrl_vs_subsets);
if pCtrl_vs_subsets == 0
    fprintf(' (report as p < %.2g; no subset as extreme as control)\n', 1/(nIter+1));
end

%% --- 5) Optional visualization ---
figure('Color','w', 'Name', sprintf('Donor subsets vs control — %s', feat));
histogram(subsetSlopes, 40, 'Normalization','pdf','FaceColor',colDon);
hold on;
yl = ylim;
if ~isnan(slopeCtrl)
    xline(slopeCtrl,'LineWidth',2, 'LineStyle','--', 'Color',colCtrl, 'Label','Controls');
end
ylim(yl);
xlabel(sprintf('Slope of %s (group mean per bin)', feat));
ylabel('Density');
title(sprintf('%s — donor subset slopes (subset=%d, nIter=%d)', feat, subsetSize, nIter));
hold off;


end

q = mafdr(ppsDon,'BHFDR',true);   % Benjamini–Hochberg FDR-adjusted p-values

% Or threshold at alpha = 0.05:
alpha  = 0.05;
sigFDR = q < alpha




%% helpers
function [xBins, mu, sem] = group_bin_stats_table(T, featureName)
% T: table already filtered to one group & one phase
% featureName: string, name of numeric column
% Returns:
%   xBins – unique win_index values
%   mu    – mean(feature) per bin
%   sem   – standard error of mean per bin

    if isempty(T)
        xBins = [];
        mu    = [];
        sem   = [];
        return;
    end

    if ~ismember(featureName, string(T.Properties.VariableNames))
        error('Feature "%s" not found.', featureName);
    end

    if ~isnumeric(T.(featureName))
        error('Feature "%s" must be numeric.', featureName);
    end

    xBins = unique(T.win_index);
    xBins = xBins(:)';

    mu  = nan(size(xBins));
    sem = nan(size(xBins));

    for i = 1:numel(xBins)
        b = xBins(i);
        rows = (T.win_index == b);
        vals = T.(featureName)(rows);

        vals = vals(~isnan(vals));
        if ~isempty(vals)
            mu(i)  = mean(vals);
            sem(i) = std(vals) / sqrt(numel(vals));
        end
    end
end


function subsetSlopes = donor_subset_slopes_table(Tdon, featureName, subsetSize, nIter)

    if isempty(Tdon)
        subsetSlopes = [];
        return;
    end

    subjIDs = unique(Tdon.participant);
    nDonors = numel(subjIDs);

    if subsetSize > nDonors
        error('subsetSize (%d) > number of donors (%d).', subsetSize, nDonors);
    end

    subsetSlopes = nan(nIter,1);

    for k = 1:nIter
        thisSubj = randsample(subjIDs, subsetSize);
        Ts = Tdon(ismember(Tdon.participant, thisSubj), :);

        [xSub, muSub] = group_bin_stats_table(Ts, featureName);  %#ok<ASGLU>
        if numel(xSub) >= 2
            mdl = fitlm(xSub, muSub);
            subsetSlopes(k) = mdl.Coefficients.Estimate(2);
        else
            subsetSlopes(k) = NaN;
        end
    end

    subsetSlopes = subsetSlopes(~isnan(subsetSlopes));
end


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
        'MinuteVentilation', 'PercentOfBrethsWithExhalePause', 'PercentOfBrethsWithInhalePause' ...
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