function [before, after,donation,NCbefore, NCafter,NCdonation]=extract_control_during(i,norm, IntLength,subjData)

T=subjData(i).Data;
resp_stereo=table2array(T(:,[3 4]));

resp_stereo(isnan(resp_stereo(:,1)),:)=[];
sum_resp=sum(resp_stereo,2);

if norm
sum_resp=zscore(sum_resp);
end

sixHzParticipants = {'045', '067','069', 'control007','control020'};  
    pid = subjData(i).code;


    if strcmpi('095',subjData(i).code)
sum_resp1=sum_resp(1:48030);
sum_resp2=sum_resp(48031:end);
sum_resp2_resampled=resample(sum_resp2,25,6);
sum_resp=[sum_resp1;sum_resp2_resampled];

resp1=resp_stereo(1:48030,:);
resp2=resp_stereo(48031:end,:);
resp2_resampled=resample(resp2,25,6);
resp_stereo=[resp1;resp2_resampled];
    end

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
 start1=subjData(i).in;%60*Fs;
%start1=2*60*Fs;
 stop1=start1+duration-1;

% start2=length(sum_resp)-duration+1;
 start2=subjData(i).out-60*Fs;

stop2=start2+duration-1;

  % start2=subjData(i).in_chair; %gives 85%
  % if isempty(start2)
  %     start2=subjData(i).out+2*60*Fs;
  % end
  % stop2=start2+duration-1;



  if stop2>length(sum_resp)
     stop2=length(sum_resp);
    fprintf('end interval are shorter than wanted\n');
    fprintf('interval length is %.2f minutes\n', length(sum_resp(start2:stop2))/(25*60));
  end
 
 %  stop2=length(sum_resp);
%  start2=length(sum_resp)-duration;

% (stop1-start1)/(Fs*60)
% (stop2-start2)/(Fs*60)



before=sum_resp(start1:stop1);
after=sum_resp(start2:stop2);
donation=sum_resp(subjData(i).in:subjData(i).out);


if start1<subjData(i).walk
    fprintf([subjData(i).code ':interval includes walk time \n ']);
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
