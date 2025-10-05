function [before, after,donation,NCbefore, NCafter,NCdonation]=extract_timings_needle(i,norm, IntLength,subjData)

%directory='hol_data';

% if size(subjData(i).Holter,1)==1
% csv=[directory '/' subjData(i).Holter];
%  T = readtable(csv);
% else
% files=subjData(i).Holter;
% T=[];
% for f=1:size(files,1)
%     csv=[directory '/' files(f,:)];
%  TT = readtable(csv);
%  T=[T;TT];
% end
% end

T=subjData(i).Data;
resp_stereo=table2array(T(:,[3 4]));

if std(resp_stereo(:,1))<25
    warning('check if R nostril recording is valid')
     figure
    plot(resp_stereo(1:5000,:))
    title([subjData(i).code ' std: ' num2str(std(resp_stereo(:,1))) ])
end

if std(resp_stereo(:,2))<25
    warning('check if L nostril recording is valid')
        figure
    plot(resp_stereo(1:5000,:))
    title([subjData(i).code ' std: ' num2str(std(resp_stereo(:,2))) ])
end

resp_stereo(isnan(resp_stereo(:,1)),:)=[];
sum_resp=sum(resp_stereo,2);

if norm
sum_resp=zscore(sum_resp);
end

sixHzParticipants = {'045', '067','069'};  
    pid = subjData(i).code;


if ismember(pid, sixHzParticipants)
        Fs=6;
else
    Fs=25;
end
    total_time_seconds=length(sum_resp)/Fs;
total_time_minutes=total_time_seconds/60;

if total_time_minutes<40
    warning('check %s sampling rate',subjData(i).code)
end
% figure;
% plot(sum_resp)
% hold on
% xline(in,'LineWidth',1.5,'LineStyle','--');
% xline(out,'LineWidth',1.5,'LineStyle','--');
% xticks(0:5*60*25:length(sum_resp));
% xticklabels(string(0:5:total_time_minutes) + " min");
% title([ 'subj' subjN])
duration=IntLength*60*Fs;
%stop1=in;
% start1=in-(60*25)-duration+1;
 start1=1*60*Fs;
 stop1=start1+duration-1;

 stop2=subjData(i).out+duration-1;
 start2=subjData(i).out;


 if stop2>length(sum_resp)
     stop2=length(sum_resp);
    fprintf('end interval are shorter than wanted\n');
 end
 
 %  stop2=length(sum_resp);
%  start2=length(sum_resp)-duration;

% (stop1-start1)/(Fs*60)
% (stop2-start2)/(Fs*60)

before=sum_resp(start1:stop1);
after=sum_resp(start2:stop2);
donation=sum_resp(subjData(i).in:subjData(i).out);

if stop1>subjData(i).in
    fprintf('interval includes donation time \n ');
end


NCbefore=resp_stereo(start1:stop1,:); 
NCafter=resp_stereo(start2:stop2,:);
NCdonation=resp_stereo(subjData(i).in:subjData(i).out,:);


if ismember(pid, sixHzParticipants)
before=resample(before, 25, 6);
after=resample(after, 25, 6);
donation=resample(donation, 25, 6);
NCbefore=resample(NCbefore, 25, 6);
NCafter=resample(NCafter, 25, 6);
NCdonation=resample(NCdonation, 25, 6);
end

end




% % % function [before, after,donation]=NC_analysis(i,IntLength,subjData)
% % % 
% % % directory='hol_data';
% % % 
% % % if size(subjData(i).Holter,1)==1
% % % csv=[directory '/' subjData(i).Holter];
% % %  T = readtable(csv);
% % % else
% % % files=subjData(i).Holter;
% % % T=[];
% % % for f=1:size(files,1)
% % %     csv=[directory '/' files(f,:)];
% % %  TT = readtable(csv);
% % %  T=[T;TT];
% % % end
% % % end
% % % 
% % % 
% % % resp_stereo=table2array(T(:,[4 3]));
% % % 
% % % 
% % % before=resp_stereo(start1:stop1,:);
% % % after=resp_stereo(start2:stop2,:);
% % % donation=resp_stereo(subjData(i).in:subjData(i).out,:);
% % % end
