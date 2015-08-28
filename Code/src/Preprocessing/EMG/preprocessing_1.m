function [MVC]=preprocessing_MCC(data,sr)
% 'preprocesso'

% bandwidth filtering-rectification-low pass filtering

% Avoid offset
[nchannels, nsamples] = size(data);
N = sr/1000;
data_1000Hz = zeros(nchannels,floor(nsamples/N));
for i=1:nchannels
    data_1000Hz(i,:)=decimate(data(i,:), N);
end
%Butterworth 7th order
%Highpass 500 Hz 
cof=500;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'high'); %Butterworth 7th order
data_h=filtfilt(B,A,data_1000Hz');     
 
% Rectification
data_rett=abs(data_h);

%Low pass 10Hz
cof=10;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'low'); %Butterworth 7th order
data_hrl=filtfilt(B,A,data_rett);   

       
% MVC
% For each muscle selct the best burst
for j=1:size(data,1)
 
figure
hold on
plot(data_hrl(:,j));
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
x=data_hrl(inizio:fine,j);   
MVC(j)=max(x);
end

end