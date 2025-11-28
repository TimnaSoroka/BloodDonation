function [before, after,donation,NCbefore, NCafter,NCdonation]=extract_timings_needle_walk_in_chair2(i,norm, IntLength,subjData)

T=subjData(i).Data;
resp_stereo=table2array(T(:,[3 4]));

resp_stereo(isnan(resp_stereo(:,1)),:)=[];
sum_resp=sum(resp_stereo,2);

if norm
sum_resp=zscore(sum_resp);
end

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
    
sixHzParticipants = {'045', '067','069', 'control007','control020'};  
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

walk=subjData(i).walk;
in_chair=subjData(i).in_chair;
out_idx=subjData(i).out;

X=sum_resp;
if isempty(walk)
    before=X(2*60*Fs:2+IntLength*60*Fs);
    NCbefore=resp_stereo(2*60*Fs:(2+IntLength)*60*Fs,:);
elseif walk-IntLength*60*Fs>1
    before=X(walk-IntLength*60*Fs:walk);
    NCbefore=resp_stereo(walk-IntLength*60*Fs:walk,:);
elseif walk-IntLength*60*Fs<1
    before=X(1:walk);
    NCbefore=resp_stereo(1:walk,:);
end

       if  isempty(in_chair) 
           after=X(out_idx+2*60*Fs:out_idx+(2+IntLength)*60*Fs);
           NCafter=resp_stereo(out_idx+2*60*Fs:out_idx+(2+IntLength)*60*Fs,:);
       elseif in_chair+IntLength*60*Fs<length(X)
          after=X(in_chair:in_chair+IntLength*60*Fs);
         NCafter=resp_stereo(in_chair:in_chair+IntLength*60*Fs,:);
       else
              after=X(in_chair:end);
         NCafter=resp_stereo(in_chair:end,:);
       end

% dur=IntLength*60*Fs;
% 
%    start1=1*60*Fs;
%    stop1=start1+dur-1;

%    stop1=subjData(i).walk-1*60*Fs;
% if isempty(stop1)
%     stop1=subjData(i).in-2*60*Fs;
% end
% 
%  start1=stop1-dur;

% 
% 
%  start2=subjData(i).in_chair;
%   if isempty(start2)
%       start2=subjData(i).out+2*60*Fs;
%   end
%  stop2=start2+dur-1;
% 
%   if stop2>length(sum_resp)
%      stop2=length(sum_resp);
%     fprintf('end interval are shorter than wanted\n');
%     fprintf('interval length is %.2f minutes\n', length(sum_resp(start2:stop2))/(25*60));
%   end
% 
%  %  stop2=length(sum_resp);
% %  start2=length(sum_resp)-duration;
% 
% % (stop1-start1)/(Fs*60)
% % (stop2-start2)/(Fs*60)
% 
% before=sum_resp(start1:stop1);
% after=sum_resp(start2:stop2);
 donation=sum_resp(subjData(i).in+30*Fs:subjData(i).out);
% 
% if stop1>subjData(i).walk
%     fprintf('interval includes walk time \n ');
% end
% 
% 
% NCbefore=resp_stereo(start1:stop1,:); 
% NCafter=resp_stereo(start2:stop2,:);
 NCdonation=resp_stereo(subjData(i).in+30*Fs:subjData(i).out,:);


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
