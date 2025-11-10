close all
clear
rng(50)

%addpath '/Users/timnas/Documents/breathmetrics-master'

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

%load('/Users/timnas/Documents/BloodDonation/BloodDonation-main/Holter_timings.mat');
load('Holter_timings_controls.mat');
%subjData([16,25,29,38,90])=[]; %less than 10 minutes after in_chair 16,25,29,38,90

% ----------------------------- EXAMPLE USAGE -----------------------------
% {
% Suppose you have 100 participants. Each phase is a cell array of vectors.
% Here we synthesize data to demonstrate:

Fs = 25;                 % 25 Hz sampling
N  = size(subjData,2);                % participants
beforeCell = cell(N,1);
duringCell = cell(N,1);
afterCell  = cell(N,1);

noiseThreshold = 0.05;         % adjust as needed
device_name    = 'Holter';     % your NC function expects 'Holter' or 'Mustache'


%Make ~15 min before (900 s), ~12 min during (variable), ~12 min after
for i = 1:N

        T = subjData(i).Data;                        
    resp_stereo = table2array(T(:,[3 4]));     


    rec_minutes = size(T,1) / (60*Fs);

    % Handle the 6 Hz cases by resampling to 25 Hz AND scaling in/out indices
    if rec_minutes < 20
        warning('check_sampling_rate for participant %s', subjData(i).code);
        Fs_raw = 6;
                Fs_new = Fs;
        X = sum(resp_stereo,2);
        X = zscore(resample(X, Fs_new, Fs_raw));
                resp_stereo=resample(resp_stereo, 25, 6);
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



    beforeCell{i} =  resp_stereo(1*60*25:in_idx-1,:);
    duringCell{i} =  resp_stereo(in_idx:out_idx,:);
    afterCell{i}  =  resp_stereo(out_idx:end,:);


        %  this option gives 80%
      % if ~isempty(walk)
      %     beforeCell{i}=X(walk-10*60*25:walk);
      % else
      %     beforeCell{i}=X(1*60*25:11*60*25);
      % end
      % if  isempty(in_chair)
      %     afterCell{i}=X(end-11*60*25:end-1*60*25);
      % else
      %     afterCell{i}=X(in_chair:in_chair+10*60*25);
      % end

end
% 
nWindows  = 5;
winFrac   = 0.2;     % each window spans 20% of that participant's phase
minWinSec = 90;      % don’t go below 45 s per window (for BM stability)

featureFcn = @(seg,Fs) nc_feats(seg, Fs, noiseThreshold, device_name);

[winTableN, winIdxN] = window_phases_make_bins_fixedN( ...
    beforeCell, duringCell, afterCell, Fs, ...
    nWindows, 'frac', winFrac, minWinSec, featureFcn);

meta = {'participant','phase','win_index','t_start_s','t_end_s','phase_frac_start','phase_frac_end'};
featList = setdiff(string(winTableN.Properties.VariableNames), meta);

results = table();


%% === PLOT ONLY THE DURING SLOPE (mean±SE + fixed-effect fit) ===
meta = {'participant','phase','win_index','t_start_s','t_end_s','phase_frac_start','phase_frac_end'};
featList = setdiff(string(winTableN.Properties.VariableNames), meta);

for f = 1:numel(featList)
    feat = featList(f);
    if ~isnumeric(winTableN.(feat)), continue; end

    % Subset DURING rows and keep only needed cols
    T = winTableN(winTableN.phase=="during", ["participant","win_index", feat]);
    if isempty(T), continue; end
    T.participant = categorical(T.participant);
    T.win_index = double(T.win_index);

    % Per-window group mean/SE
    [uIdx, ~, gid] = unique(T.win_index);
    mu = accumarray(gid, T.(feat), [], @(v) mean(v,'omitnan'));
    sd = accumarray(gid, T.(feat), [], @(v) std(v,'omitnan'));
    n  = accumarray(gid, T.(feat), [], @(v) sum(isfinite(v)));
    se = sd ./ max(sqrt(n),1);

    [uIdx, ord] = sort(uIdx); mu = mu(ord); se = se(ord);

    % Mixed-effects (same as your function) to get fixed-effect slope
    out = mixed_effects_during(winTableN, feat, 'time','linear');
    fixed = out.fixed;
    b0 = fixed.Estimate(strcmp(fixed.Name,'(Intercept)'));
    b1 = fixed.Estimate(strcmp(fixed.Name,'win_c'));
    p1 = fixed.pValue(strcmp(fixed.Name,'win_c'));

    % NOTE: mixed_effects_during centers time by mean(win_index) in DURING;
    % we recompute that mean here to build the prediction on original x.
    muIdx = mean(T.win_index,'omitnan');

    % Prediction line on a dense grid over DURING window indices
    xg = linspace(min(uIdx), max(uIdx), 200);
    yhat = b0 + b1*(xg - muIdx);

    % ---- Plot (DURING only) ----
    figure('Color','w','Name',sprintf('During slope — %s',feat));
    hold on;
    % shaded SE
    if numel(uIdx) > 1
        lo = mu - se; hi = mu + se;
        fill([uIdx; flipud(uIdx)], [lo; flipud(hi)], [0.4660 0.6740 0.1880], ...
            'FaceAlpha',0.12,'EdgeColor','none');
    end
    % mean points/line
    plot(uIdx, mu, 'o-', 'LineWidth', 1.8,  'Color',[0.4660 0.6740 0.1880] );%[0 0.447 0.741]

    % % fitted trend
    % plot(xg, yhat, '-', 'LineWidth', 2.4);
    % 
    % grid on;
    xlabel('Window index');
    ylabel(strrep(char(feat),'_',' '), 'Interpreter','none');
    title(sprintf('During: slope = %.4g, p = %.3g', b1, p1), 'Interpreter','none');

    % y-lims with a little padding
    yAll = [mu(:); yhat(:)];
    yAll = yAll(isfinite(yAll));
    if ~isempty(yAll)
        pad = 0.05*range(yAll); if pad==0, pad = 0.05*max(1e-6,abs(mean(yAll))); end
        ylim([min(yAll)-pad, max(yAll)+pad]);
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

function featureSummaries = plot_breathmetrics_by_index(winTable, varargin)
% Mean ± SE per phase, x = window index, SAME y across phases.
% Y-limits come from the plotted mean±SE envelopes (not raw values).
% Returns:
%   featureSummaries.(featureName) = table(phase, win_index, mean, se)

    % ---- params ----
    p = inputParser;
    addParameter(p,'clip',[0 100],@(v)isnumeric(v)&&numel(v)==2&&v(1)>=0&&v(2)<=100);
    addParameter(p,'semult',1.0,@(v)isnumeric(v)&&isscalar(v)&&v>0); % 1.96 ≈ 95% CI
    parse(p,varargin{:});
    qclip  = p.Results.clip;
    semult = p.Results.semult;

    phases = {'before','during','after'};
    metaCols = {'participant','phase','win_index','t_start_s','t_end_s','phase_frac_start','phase_frac_end'};
    allCols  = string(winTable.Properties.VariableNames);
    featCols = setdiff(allCols, string(metaCols));

    % stable order if present
    preferredOrder = [ ...
        "AverageExhaleDuration","AverageExhalePauseDuration","AverageExhaleVolume", ...
        "AverageInhaleDuration","AverageInhalePauseDuration","AverageInhaleVolume", ...
        "AverageInterBreathInterval","AveragePeakExpiratoryFlow","AveragePeakInspiratoryFlow", ...
        "AverageTidalVolume","BreathingRate", ...
        "CoefficientOfVariationOfBreathVolumes","CoefficientOfVariationOfBreathingRate", ...
        "CoefficientOfVariationOfExhaleDutyCycle","CoefficientOfVariationOfExhalePauseDutyCycle", ...
        "CoefficientOfVariationOfInhaleDutyCycle","CoefficientOfVariationOfInhalePauseDutyCycle", ...
        "DutyCycleOfExhale","DutyCycleOfExhalePause","DutyCycleOfInhale","DutyCycleOfInhalePause", ...
        "MinuteVentilation","PercentOfBreathsWithExhalePause","PercentOfBreathsWithInhalePause" ...
    ];
    featCols = [featCols(ismember(featCols,preferredOrder)), featCols(~ismember(featCols,preferredOrder))];
    featCols = unique(featCols,'stable');

    featureSummaries = struct();

    for f = 1:numel(featCols)
        feat = featCols(f);
        if ~isnumeric(winTable.(feat)), continue; end

        % --- per-feature summary holder ---
        groupSummary = table();

        % ---- compute per-phase mean/SE by window index ----
        S = struct(); envelopes = [];
        for pidx = 1:numel(phases)
            ph = phases{pidx};
            sub = winTable(winTable.phase==ph, ["win_index", feat]);
            if isempty(sub)
                S.(ph).idx=[]; S.(ph).mu=[]; S.(ph).se=[]; 
                continue;
            end
            [uIdx, ~, gid] = unique(sub.win_index);
            mu = accumarray(gid, sub.(feat), [], @(v) mean(v,'omitnan'));
            sd = accumarray(gid, sub.(feat), [], @(v) std(v, 'omitnan'));
            n  = accumarray(gid, sub.(feat), [], @(v) sum(isfinite(v)));
            se = (sd ./ max(sqrt(n),1)) * semult;

            [uIdx, ord] = sort(uIdx); mu = mu(ord); se = se(ord);
            S.(ph).idx = uIdx; S.(ph).mu = mu; S.(ph).se = se;

            envelopes = [envelopes; mu-se; mu+se]; %#ok<AGROW>

            % generic varnames
            T = table( ...
                repmat(categorical({ph}, phases), numel(uIdx),1), ...
                uIdx, mu, se, ...
                'VariableNames', {'phase','win_index','mean','se'});
            groupSummary = [groupSummary; T]; %#ok<AGROW>
        end

        % store in struct
        safeName = matlab.lang.makeValidName(char(feat));
        featureSummaries.(safeName) = groupSummary;

        % ---- y-lims from mean±SE envelopes only ----
        env = envelopes(isfinite(envelopes));
        if isempty(env)
            yLimFeat = [0 1];
        else
            if qclip(1)==0 && qclip(2)==100
                yLow = min(env); yHigh = max(env);
            else
                yLow  = quantile(env, qclip(1)/100);
                yHigh = quantile(env, qclip(2)/100);
            end
            if ~isfinite(yLow) || ~isfinite(yHigh) || yLow==yHigh
                pad = max(1e-6, abs(yLow)*0.05);
                yLimFeat = [yLow - pad, yHigh + pad];
            else
                pad = 0.03*(yHigh - yLow);
                yLimFeat = [yLow - pad, yHigh + pad];
            end
        end

        % ---- plot: shaded SE + mean, horizontal layout ----
        figure('Color','w','Name',char(feat)); %#ok<LFIG>
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');  % horizontal layout

        for pidx = 1:numel(phases)
            ph = phases{pidx};
            nexttile; hold on;

            idx = S.(ph).idx; mu = S.(ph).mu; se = S.(ph).se;
            if ~isempty(idx)
                lo = mu - se; hi = mu + se;
                if numel(idx)==1
                    idx = [idx; idx+0.001]; mu=[mu;mu]; lo=[lo;lo]; hi=[hi;hi];
                end

                fill([idx; flipud(idx)], [lo; flipud(hi)],[0 0.447 0.741], ...
                     'FaceAlpha',0.18,'EdgeColor','none');  % shaded SE
                plot(idx, mu, '-', 'LineWidth', 2, 'Color', [0 0.447 0.741]);
                xlim([min(idx) max(idx)]);
            else
                axis off; title(sprintf('%s (no data)', ph),'Interpreter','none');
            end

            grid on; ylim(yLimFeat);
            xlabel('Window index (1..N)');
            ylabel(strrep(char(feat),'_',' '), 'Interpreter','none');
            title(ph,'Interpreter','none'); % phase only
        end

        sgtitle(strrep(char(feat),'_',' '), 'FontWeight','bold', 'Interpreter','none');
    end
end

function plot_participant_breathmetrics(winTable, participantID, featureName)
% Plot one participant, three phases side by side (same y-limits).
%
% Inputs:
%   winTable      – table from window_phases_make_bins_fixedN
%   participantID – numeric participant index
%   featureName   – string, e.g. "BreathingRate"

    phases = {'before','during','after'};
    sub = winTable(winTable.participant==participantID,:);
    if isempty(sub)
        warning('No data for participant %d', participantID);
        return;
    end

    % === Compute unified y-limits across phases ===
    yAll = sub.(featureName);
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        yLimFeat = [0 1];
    else
        pad = 0.05 * range(yAll);
        if pad == 0, pad = max(1e-6, 0.05 * abs(mean(yAll))); end
        yLimFeat = [min(yAll)-pad, max(yAll)+pad];
    end

    % === Create figure ===
    figure('Color','w','Name',sprintf('P%02d - %s',participantID,featureName));
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    for p = 1:numel(phases)
        ph = phases{p};
        nexttile; hold on;
        dat = sub(sub.phase==ph,:);
        if isempty(dat)
            axis off; title([ph ' (no data)'],'Interpreter','none');
            continue;
        end

        x = dat.win_index;
        y = dat.(featureName);
        plot(x, y, 'o-', 'LineWidth', 1.8, 'Color', [0 0.447 0.741]);
        xlim([min(x) max(x)]);
        ylim(yLimFeat);

        xlabel('Window index');
        ylabel(strrep(featureName,'_',' '), 'Interpreter','none');
        title(ph, 'Interpreter','none');
        grid on;
    end

    sgtitle(sprintf('Participant %d — %s',participantID,featureName), ...
        'FontWeight','bold','Interpreter','none');
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
    fprintf('\n=== Mixed Effects (during) — %s | time=%s ===\n', fname, modeTime);
    disp(out.anova);
    disp(out.fixed);
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

function val = safeget(s, field, defaultVal)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end



    
