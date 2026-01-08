


close all
clear
rng(50)
norm=1; %1=80%, 0=50%
IntLength=5; % 3 is the best
Fs=25; 
%
 load('Holter_timings.mat');
load('tempos_results.mat') %5min

%load('Holter_timings_controls.mat');

%%
for i=1:size(subjData,2)
%[before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle(i,norm, IntLength,subjData); %_needle_walk_in_chair
%    [before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle_lay(i,norm, IntLength,subjData); %_needle_walk_in_chair
[before{i},after{i},donation{i},NCbefore{i},NCafter{i},NCdonation{i}]=extract_timings_needle_walk_in_chair2(i,norm,IntLength,subjData); %_needle_walk_in_chair

end

%vals_before=calculate_before_after(before);
% vals_before_c=calculate_before_after_catch(before);
% 
noiseThreshold=0.05;

for i=1:size(NCbefore,2)
    if length(after{i})<5*Fs*60
        continue
    else
[Laterality_IndexB(i),BmeasureResults(i)]=NasalCycleParameters(NCbefore{i},Fs,noiseThreshold,'Holter');
    end
end



% vitIDs = strings(size(vitalIDs));
% for i = 1:numel(vitalIDs)
%     if isempty(vitalIDs{i})
%         vitIDs(i) = "";          % keep as empty string
%     else
%         vitIDs(i) = string(vitalIDs{i});
%     end
% end
% vitIDs = strtrim(vitIDs);

 %%
 

T=readcell('/Users/timnas/Documents/projects/BloodDonation/BloodDonationTable102.xlsx');
demo=cell2table(T([3,4,9,11,30,31,32,33,34],2:end-20)','VariableNames',T([3,4,9,11,30,31,32,33,34],1));
demo(5,:)=[];

feat=struct2table(BmeasureResults);

sx = demo.sex;
    s = lower(strtrim(string(sx)));
    sex_num = nan(height(demo),1);
    sex_num(ismember(s, ["m","זכר","1"])) = 1;
    sex_num(ismember(s, ["f","נקבה","0"])) = 0;

demo.sex = sex_num;

% --- numeric feature vars only ---
for v = 1:width(feat)
    if iscell(feat.(v))
        idx = cellfun(@isempty, feat.(v));
        feat.(v)(idx) = {nan};
    end
end


featVars = feat.Properties.VariableNames;


% ---------- age/weight correlations ----------
preds=["Age","Weight",'Spo1st','Pulse1st','Height'];
for pred =1:size(preds,2) 
    x = demo.(preds(pred));
    if iscell(x)
     x(strcmpi('NAN',x))={nan};
     x =cell2mat(x);
    end
    for j = 1:numel(featVars)
        y = feat.(featVars{j});
     y =cell2mat(y);
        [r(pred,j),p(pred,j)] = corr(x, y, 'Rows','complete');
        if p(pred,j)<0.05
            fprintf('corr between %s and %s is sig. R= %s, p=%s\n',featVars{j},preds(pred),num2str(r(pred,j)),num2str(p(pred,j)))
        end
    end
end

% ---------- sex group tests ----------
g = demo.sex;
for j = 1:numel(featVars)
    y = feat.(featVars{j});
         y =cell2mat(y);
    y0 = y(g==0);
    y1 = y(g==1);
            [~,p(pred+1,j)] = ttest2(y0, y1); % Welch
            if p(pred+1,j)<0.05
                            fprintf('ttest between sex and %s is sig. p=%s\n',featVars{j},num2str(p(pred+1,j)))
            end
end    

firstD = demo.donated_for_the_first_time;
    firstD(strcmpi('NAN',firstD))={nan};
    s = lower(strtrim(string(firstD)));
    first_time = nan(height(demo),1);
    first_time(ismember(s, ["m","לא"])) = 1;
    first_time(ismember(s, ["f","כן"])) = 0;

demo.donated_for_the_first_time = first_time;

% T = cell2table(rows, 'VariableNames', ...
%     {'Predictor','Feature','Test','Statistic','P','Effect_rankBiserial'});
% 
% % sort by p-value (smallest first)
% T = sortrows(T, 'P', 'ascend');
% 
% res = struct();
% res.table = T;
% res.features = featVars;
