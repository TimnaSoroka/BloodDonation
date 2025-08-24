function [Resp] = hilbert24 (Data, Fs, noiseThreshold)

% this function takes the raw-data and returns it with
% one-sample-per-minute. The raw data should be sampled in 180 ms between
% samples (1000/180 Hz)

% Updated by Maya - 15/11
% Updated by Lavi & Roni - 12/8
% Updated by Roni 10/3/14
% Updated by Lior - 17/01/17


bandPassFreqs=[0.5 5];



% Initial filter
inputDataLength=size(Data,1);
myWindow=zeros(inputDataLength,1);
max_freq=Fs/2;
df=Fs/inputDataLength;
center_freq=mean(bandPassFreqs);
filter_width=diff(bandPassFreqs);
x=0:df:max_freq;
gauss=exp(-0.5*(x-center_freq).^2);
cnt_gauss = round(center_freq/df);
flat_padd = round(filter_width/df);  % flat padding at the max value of the gaussian
padd_left = floor(flat_padd/2);
padd_right = ceil(flat_padd/2);
aux = [gauss((padd_left+1):cnt_gauss) ones(1,flat_padd) gauss((cnt_gauss+1):(end-padd_right))];
myWindow(1:length(aux))=power(aux,5);
myWindowMat = repmat(myWindow,1,size(Data,2));
fftRawEeg=fft(Data);
filtSignal=ifft(fftRawEeg.*myWindowMat,'symmetric');

Data=filtSignal;

% First - remove mean for every 10 minutes and compute hilbert transform

% L = length(Data);
T_mean = floor(5*60*Fs); % check mean every 10 minutes
K = floor(inputDataLength/T_mean);

% h1 = figure;
% h2 = figure;
%
% hold all;
clear Data_wo_mean h

% Initially we subtracted mean for battery voltage change and apparently
% this is not needed so we omitted that filter.

for j = 1:2
    Data_wo_mean{j} = [];
    for i = 1:K
        temp_vec = Data((i-1)*T_mean+1:i*T_mean,j) - mean(Data((i-1)*T_mean+1:i*T_mean,j));
        Data_wo_mean{j} = [Data_wo_mean{j} temp_vec'];
    end

%     figure(h1)
%     hold all
%     plot(Data_wo_mean{j});

data_after_hilbert = hilbert(Data_wo_mean{j});
    h(:,j) = abs(data_after_hilbert);
    Ph(:,j) = angle(data_after_hilbert);

%     figure(h2)
%     a(j) = subplot(2,1,j);
%     plot(Data_wo_mean{j});
%     hold all
%     plot(h(:,j));
%     %plot(Ph(:,j));
%
%     title(sprintf('Data #%d',j));

end

% figure properties
% figure(h2)
% linkaxes([a(2) a(1)],'x');

% figure(h1)
% legend('Data1 - mean removed','Data2 - mean removed')

% Second - get the peaks of the hilbert transform

% Notice I don't allow a peak to be close more than 10 samples to the next
% one to exclude the un-wanted jumps in values. I think this value can be
% set for every "breather" based on minimal breathing tempo (here we can go
% up to 15).
for j = 1:2
    [pks{j},locs{j}] = findpeaks(h(:,j),'MINPEAKDISTANCE',10);
    %     figure(h2)
    %     subplot(2,1,j)
    %     hold all
    %     plot(locs{j},pks{j},'m*')
end

% Third - average the peaks for 1 minutes -can play with reasonable time-frame

L = length(Data);
bin=1; % bin size in min to average peaks in
frame_length = bin*60; % (sec) check mean every "bin" minutes
T_mean = floor(frame_length*Fs);
K = floor(L/T_mean);

% h3 = figure;
NewLocs=[];

for j = 1:2
    val_per_frame{j} = [];
    Resp_freq{j} = [];
    for i = 1:K
        indices_peaks = find(locs{j} > (i-1)*T_mean & locs{j} < i*T_mean);
        if isempty(indices_peaks)
            temp_val = 0;
            no_of_pks=0;
        else
            temp_val = mean(pks{j}(indices_peaks));
            if temp_val>noiseThreshold
                no_of_pks=length(indices_peaks);
            else
                if j==1
                    tmp_j=2;
                else
                    tmp_j=1;
                end
                indices_peaks = find(locs{tmp_j} > (i-1)*T_mean & locs{tmp_j} < i*T_mean);
                no_of_pks=length(indices_peaks);
            end
        end
        val_per_frame{j} = [val_per_frame{j} temp_val];
        Resp_freq{j} = [Resp_freq{j} no_of_pks];
    end
end

NewLocs=T_mean/2:T_mean:K*T_mean-T_mean/2;

% figure(3)
% plot(NewLocs,val_per_frame{1},'.:');
% hold all;
% plot(NewLocs,val_per_frame{2},'.:');
% if ~isempty(night) stem(night,1,'LineWidth',3); end
% if ~isempty(morning)stem(morning,1,'LineWidth',3);end
% title('Peaks mean value per frame')
% legend('Right','Left');

vpr1=val_per_frame{1};
vpr2=val_per_frame{2};

noise1=vpr1<noiseThreshold;
noise2=vpr2<noiseThreshold;
noise=noise1&noise2;

Resp=[vpr1; vpr2; NewLocs; Resp_freq{1}; Resp_freq{2}];

Resp(:,noise)=[];

% h4=figure;
% plot(Resp(1,:),'.:');
% hold all;
% plot(Resp(2,:),'.:');
%rawDataFolder = fileparts(RawDataFilePath);
%ToSave=fullfile(rawDataFolder, ['hilbert-' ,SubjName, '.mat']);
%save(ToSave, 'Resp'); % To save only Resp into the file
end