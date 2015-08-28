function [data_p]=preprocessing(data,sr,MVC,nchan)
% 'preprocesso'

% bandwidth filtering-rectification-low pass filtering

% Avoid offset
for j=1:nchan
    data(j,:)=data(j,:)-mean(data(j,:));
end

%Butterworth 7th order
%Highpass 40 Hz 
cof=50;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'high'); %Butterworth 7th order
data_h=filtfilt(B,A,data');     
 
% Rectification
data_rett=abs(data_h);

%Low pass 20Hz
cof=20;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'low'); %Butterworth 7th order
data_hrl=filtfilt(B,A,data_rett);   
       
% Normalizzazione
for j=1:nchan
 data_p(:,j)=data_hrl(:,j)./MVC(j);
end

end