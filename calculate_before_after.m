    function [z_breath_values,vars]=calculate_before_after(before,IntLength)

for i=1:size(before,2)
%    peaks1 = peaks_from_ts2(before{i});
%    z_breath_values(i) = calculate_z_1min(peaks1,IntLength);
Fs= size(before{i},1)/(IntLength*60);
if Fs<25
    % peaks1 = peaks_from_ts_fs(before{i},Fs);
    % 
    % plot(before{i})
    % hold on
    % plot([peaks1.PeakLocation],[peaks1.PeakValue],'bo')
    % 
    % z_breath_values(i) = calculate_z_blood(peaks1);
    continue
else
     bmObj=breathmetrics(before{i},Fs,'humanAirflow');
bmObj.estimateAllFeatures(0,'sliding',1, 0);

 % Define the variable names as field names in the struct
variableNames = {...
    'AverageExhaleDuration', 'AverageExhalePauseDuration', 'AverageExhaleVolume', 'AverageInhaleDuration', ...
    'AverageInhalePauseDuration', 'AverageInhaleVolume', 'AverageInterBreathInterval', 'AveragePeakExpiratoryFlow', ...
    'AveragePeakInspiratoryFlow', 'AverageTidalVolume', 'BreathingRate', 'CoefficientOfVariationOfBreathVolumes', ...
    'CoefficientOfVariationOfBreathingRate', 'CoefficientOfVariationOfExhaleDutyCycle', 'CoefficientOfVariationOfExhalePauseDutyCycle', ...
    'CoefficientOfVariationOfInhaleDutyCycle', 'CoefficientOfVariationOfInhalePauseDutyCycle', 'DutyCycleOfExhale', ...
    'DutyCycleOfExhalePause', 'DutyCycleOfInhale', 'DutyCycleOfInhalePause', 'MinuteVentilation', 'PercentOfBreathsWithExhalePause', ...
    'PercentOfBreathsWithInhalePause'};

values=[bmObj.secondaryFeatures.values];
% Assign the values to the struct fields
for ii = 1:length(variableNames)
    dataStruct.(variableNames{ii}) = values{ii};
end

 z_breath_values(i)= dataStruct(:);
 %z_breath_values(i,:)=values;
end
end

%vars=bmObj.secondaryFeatures.keys;
load('vars_names_BM.mat')
vars=variableNames;
end
