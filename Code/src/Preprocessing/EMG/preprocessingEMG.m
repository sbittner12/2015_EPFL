function [data_p]=preprocessingEMG(data,sr,MVC,nchan)
%%PREPROCESSINGEMG - performs preprocessing on EEG data from experiments
%  preprocessingEG()
%
%  Write a full length description
%
%  PARAMETERS:
%    exptype: string - specifies the type of experiment.  Options are:
%            'Healthy' 'IAF', 'LT', 'MJ', and 'OAF'
%
%  RETURNS:
%    specify return variables if any
%
%  @ June 2015, Sean Bittner   sbittner@andrew.cmu.edu
%
%    - adapted from file preprocessing 
%      written by Aurelie Sadaka-Stephan

[nsamples, nchannels] = size(data);
N = sr/1000;
%Butterworth 7th order
%Low pass 50 Hz 
cof=50;
Wn=(cof*2)./sr;
if Wn>1.0
    Wn=0.99;
end
[B,A] = butter(7,Wn,'low'); %Butterworth 7th order
data_hp=filtfilt(B,A,data);     

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

data_p = zeros(nsamples/N, nchannels);

% Normalizzazione
for j=1:nchannels
 data_p(:,j)=data_1000Hz(:,j)./MVC(j);
end

end