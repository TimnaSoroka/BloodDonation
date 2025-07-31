function z = calculate_z_blood(peaks)

only_inhales=peaks([peaks.PeakValue]>0);
only_exhales=peaks([peaks.PeakValue]<0);
n=fieldnames(only_inhales);
to_exclude_inhales = false(1, numel(only_inhales));
to_exclude_exhales = false(1, numel(only_exhales));

%% Exclude inhales/exhales with duration bigger than 10
for i=1:size(n,1)
    if ~strcmp(n{i}, 'Duration')
        continue;
    end
    x=[only_inhales.(n{i})];
    y=abs([only_exhales.(n{i})]);
    
    to_exclude_inhales = to_exclude_inhales | x > 10;
    to_exclude_exhales = to_exclude_exhales | y > 10;
end


for i=1:size(n,1)
    if ismember(n{i}, {'PeakLocation', 'NumberOfPeaks'})
        continue;
    end
    x=[only_inhales(~to_exclude_inhales).(n{i})];
    y=abs([only_exhales(~to_exclude_exhales).(n{i})]);
    
    m=mean(x,'omitnan');
    s=std(x,'omitnan');
    x(x>m+3*s | x<m-3*s) = nan;
    
    m=mean(y,'omitnan');
    s=std(y,'omitnan');
    y(y>m+3*s | y<m-3*s) = nan;

    clean_inhales.(n{i})=x;
    clean_exhales.(n{i})=y;
end

%% Get time between inhales
inter_inhales=diff([only_inhales.StartTime]);
inter_inhales(inter_inhales > 15) = nan;
m=mean(inter_inhales,'omitnan');
s=std(inter_inhales,'omitnan');
inter_inhales(inter_inhales>m+3*s | inter_inhales<m-3*s)=nan;

%% Get time between exhales
inter_exhales=diff([only_exhales.StartTime]);
inter_exhales(inter_exhales > 15) = nan;
m=mean(inter_exhales,'omitnan');
s=std(inter_exhales,'omitnan');
inter_exhales(inter_exhales>m+3*s | inter_exhales<m-3*s)=nan;

%% Choose if to use inhales or exhales for time between
if mean(inter_inhales, 'omitnan') < mean(inter_exhales, 'omitnan') 
    inter_breaths = inter_inhales;
else
    inter_breaths = inter_exhales;
end
breaths_rate = 1./inter_breaths;

%% Get pauses
onsets=[peaks.StartTime];
offsets=onsets+[peaks.Duration];
[~,idx_s]=sort(onsets);
in=[peaks.PeakValue]>0;
InEx_s=in(idx_s);

onsets_sorted=onsets(idx_s);
offsets_sorted=offsets(idx_s);

pauses = onsets_sorted(2:end) - offsets_sorted(1:end-1);
inhale_pauses_indices = InEx_s(1:end-1) == 1 & InEx_s(2:end) == 0 & pauses >= 0.05;
exhale_pauses_indices = InEx_s(1:end-1) == 0 & InEx_s(2:end) == 1 & pauses >= 0.05;
inhale_pause = pauses(inhale_pauses_indices);
exhale_pause = pauses(exhale_pauses_indices);

inhale_pause(inhale_pause>mean(inhale_pause)+3*std(inhale_pause)| inhale_pause<mean(inhale_pause)-3*std(inhale_pause))=nan;
exhale_pause(exhale_pause>mean(exhale_pause)+3*std(exhale_pause)| exhale_pause<mean(exhale_pause)-3*std(exhale_pause))=nan;

%% Summary everything to output variable
% z.Inhale_Peak_Count = sum([only_inhales.NumberOfPeaks]);
% z.Exhale_Peak_Count = sum([only_exhales.NumberOfPeaks]);
% z.Inhale_Count = sum(~isnan(clean_inhales.Duration));
% z.Exhale_Count = sum(~isnan(clean_exhales.Duration));
z.AverageInhaleVolume=mean([clean_inhales.Volume],'omitnan');
z.AverageExhaleVolume=mean(abs([clean_exhales.Volume]),'omitnan');
z.AverageInhaleDuration=mean([clean_inhales.Duration],'omitnan');
z.AverageExhaleDuration=mean([clean_exhales.Duration],'omitnan');
z.AveragePeakInspiratoryFlow=mean([clean_inhales.PeakValue],'omitnan');
z.AveragePeakExpiratoryFlow=mean(abs([clean_exhales.PeakValue]),'omitnan');
z.AverageInterBreathInterval=mean(inter_breaths,'omitnan');
z_rate=1./[z.AverageInterBreathInterval];
z.BreathingRate=60*z_rate;
z.AverageTidalVolume = z.AverageInhaleVolume + z.AverageExhaleVolume;
z.MinuteVentilation = z_rate .* z.AverageTidalVolume;

z.DutyCycleOfInhale=mean([clean_inhales.Duration]./[z.AverageInterBreathInterval],'omitnan');
z.DutyCycleOfExhale=mean([clean_exhales.Duration]./[z.AverageInterBreathInterval],'omitnan');

z.CoefficientOfVariationOfInhaleDutyCycle=std([clean_inhales.Duration],'omitnan')./mean([clean_inhales.Duration],'omitnan');
z.CoefficientOfVariationOfExhaleDutyCycle=std([clean_exhales.Duration],'omitnan')./mean([clean_exhales.Duration],'omitnan');

%z.COV_Inter_breath_interval=std(inter_breaths./mean(inter_breaths,'omitnan'),'omitnan');
z.CoefficientOfVariationOfBreathingRate=std(breaths_rate./mean(breaths_rate,'omitnan'),'omitnan');

z.CoefficientOfVariationOfBreathVolumes=std([clean_inhales.Volume],'omitnan')./[z.AverageInhaleVolume];
%z.COV_ExhaleVolume=std([clean_exhales.Volume],'omitnan')./[z.Exhale_Volume];

z.AverageInhalePauseDuration=mean(inhale_pause,'omitnan');
z.AverageExhalePauseDuration=mean(exhale_pause,'omitnan');

z.CoefficientOfVariationOfInhalePauseDutyCycle=std(inhale_pause,'omitnan')./mean(inhale_pause,'omitnan');
z.CoefficientOfVariationOfExhalePauseDutyCycle=std(exhale_pause,'omitnan')./mean(exhale_pause,'omitnan');
z.DutyCycleOfInhalePause=mean(inhale_pause./[z.AverageInterBreathInterval],'omitnan');
z.DutyCycleOfExhalePause=mean(exhale_pause./[z.AverageInterBreathInterval],'omitnan');

z.PercentOfBreathsWithExhalePause=length(exhale_pause)*100./(size(peaks,1)-size(only_inhales,1));
z.PercentOfBreathsWithInhalePause=length(inhale_pause)*100./size(only_inhales,1);


end