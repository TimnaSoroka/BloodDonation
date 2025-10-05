% This code use the men for training and women for testing or fat for training and less fat for testing. FS done using
% values<0.05 in the man sample.

close all
clear
rng(50)
norm=1; %1=80%, 0=50%
IntLength=5; % 5 is the best
to_plot=0;
Fs=25; 
load('Holter_timings.mat');
disp_ttests=1;

%%
subjData(91)=[]; %have short after (*technical issue)

%%
for i=1:size(subjData,2)
 [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData);
end

 %%

  noiseThreshold=0;

 for i=1:size(NCbefore,2)
     % if ismember(i,[42,63,65]) 
     %     Fs=6;
     % else
     %     Fs=25;
     % end
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters(NCbefore{i},Fs,noiseThreshold,'Holter');
[Laterality_IndexA(i),AmeasureResults(i)]=NasalCycleParameters(NCafter{i},Fs,noiseThreshold,'Holter');
%[Laterality_IndexD(i),DmeasureResults(i)]=NasalCycleParameters(NCdonation{i},Fs,noiseThreshold);
 end

  for i=1:size(NCbefore,2)
      if ~isempty([Laterality_IndexA(i).one])
 LI_ampB(i,:)=abs([Laterality_IndexB(i).one]);
 LI_B(i,:)=[Laterality_IndexB(i).one];
 LI_ampA(i,1:length(abs([Laterality_IndexA(i).one])))=abs([Laterality_IndexA(i).one]);
 LI_A(i,1:1:length(abs([Laterality_IndexA(i).one])))=[Laterality_IndexA(i).one];
      end
  end

%%

 fields=fieldnames(BmeasureResults);
for i=1:size(fields,1)
    currentfield=fields{i};
test_values = [BmeasureResults(:).(currentfield)];   
retest_values = [AmeasureResults(:).(currentfield)]; 

%[p_values(i),~,Wstat(i)] = signrank(test_values, retest_values,"method","approximate");
[~,p_values_ttest(i)] = ttest2(test_values, retest_values);

% fprintf('%s ttest p = %.2f\n',currentfield, p_values_ttest(i))
%         fprintf(['mean+std before and after' num2str(mean(test_values,'omitnan')) '±' num2str(std(test_values,'omitnan')) ',' num2str(mean(retest_values,'omitnan')) '±' num2str(std(retest_values,'omitnan')) '\n'])
end

%% big differences

vals_before=calculate_before_after(before,IntLength);
[vals_after,vars]=calculate_before_after(after,IntLength);

fields2=fieldnames(vals_after);

fields=[fields2;fields{1};fields{4};fields{6};fields{9}];


% vals_before([42,63,65])=[];
% vals_after([42,63,65])=[];

weight_vec=[subjData.Weight];
% weight_vec([42,63,65])=[];
weight_vec_con=[weight_vec,weight_vec];

code_vec={subjData.code};
%code_vec([42,63,65])=[];

X=table2array(struct2table([vals_before,vals_after]));

% BmeasureResults([42,63,65])=[];
% AmeasureResults([42,63,65])=[];
x=table2array(struct2table([BmeasureResults,AmeasureResults]));
x=x(:,[1,4,6,9]);

X=[X,x];
Y=[repmat({'before'},size(vals_before,2),1);repmat({'after'},size(vals_after,2),1)];
%%
participant=[code_vec,code_vec];
cohort=repmat({'donors'},size(Y));
%%
T= table(X,Y,participant',weight_vec_con',cohort,'VariableNames',{'X','session','participant','weight','cohort'});
%%
save('Donors_table.mat','T');