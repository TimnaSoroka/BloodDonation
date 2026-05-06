function [z_breath_values,variableNames]=catch22_feats(before)

    variableNames = {...
    'mode_5', 'mode_10', 'outlier_timing_pos', 'outlier_timing_neg', ...
    'acf_timescale', 'acf_first_min', 'low_freq_power', 'centroid_freq', ...
    'forecast_error', 'whiten_timescale', 'high_fluctuation', 'stretch_high', ...
    'stretch_decreasing', 'entropy_pairs', 'ami2', ...
    'trev', 'ami_timescale', 'transition_variance', ...
    'periodicity', 'embedding_dist', 'rs_range', 'dfa'};

%    peaks1 = peaks_from_ts2(before{i});
%    z_breath_values(i) = calculate_z_1min(peaks1,IntLength);
values=catch22(before);

 % Define the variable names as field names in the struct
% Assign the values to the struct fields
for ii = 1:length(variableNames)
    dataStruct.(variableNames{ii}) = values(ii,:);
end

 z_breath_values= dataStruct(:);
 %z_breath_values(i,:)=values;


end
