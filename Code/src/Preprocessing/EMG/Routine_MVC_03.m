clear all
clc
close all
% 
% disp('*****************************************************');
% disp('Kinematics and EMG signals ARMEO Boom-synergies - MVC');
% disp('*****************************************************');
% disp(' ');

directory_iniziale=cd;
dir_sub=input('Folder path of the subject to analyze: ','s');

name_file={'MVC.mat'};

for f=1
cd(dir_sub);
load(name_file{f});
cd(directory_iniziale);

for i=1:16-1  
    Data = getfield(EMG,'data');
    emg(:,i)=Data(:,i);
end

emg = emg';

clear Activities Data length_sec noChans samplingRate
[MVCb]=preprocessing_1(emg,3000); % MVC and the other EMG are preprocessed in the same way
clear EMG
end

cd(dir_sub);   
save MVC_proc MVCb 
  