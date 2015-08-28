function [MVC]=preprocessing_EMG_MCC(data,sr)
% 'preprocesso'

% bandwidth filtering-rectification-low pass filtering

% Avoid offset
[nchannels, nsamples] = size(data);
N = sr/1000;

%Butterworth 7th order
%Low pass 500 Hz 
cof=500;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'low'); %Butterworth 7th order
data_hp=filtfilt(B,A,data');     

%High pass 10Hz
cof=10;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'high'); %Butterworth 7th order
data_bp=filtfilt(B,A,data_hp);   
 
% Rectification
data_rett=abs(data_bp);
       
data_1000Hz = zeros(floor(nsamples/N),nchannels);
for i=1:nchannels
    data_1000Hz(:,i)=downsample(data_rett(:,i), N);
end
% MVC
% For each muscle selct the best burst
for j=1:size(data,1)
 
figure
hold on
plot(data_1000Hz(:,j));
disp('****************************************');
disp('*                                      *');
disp('*  ZOOM AND SELECT THE RIGHT INTERVAL  *');
disp('*                                      *');
disp('****************************************');
pause

[asc,ord]=ginput(2);
inizio=round(asc(1));
fine=round(asc(2));

close all


%% Local max
x=data_1000Hz(inizio:fine,j);   
MVC(j)=max(x);
end

end