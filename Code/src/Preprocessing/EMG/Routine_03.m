clear all
clc
close all

% disp('*****************************************************');
% disp('EMG signals EEG-EMG-glove');
% disp('*****************************************************');
% disp(' ');

directory_iniziale=cd;

dir_sub=input('Folder path of the subject to analyze: ','s');
name_file={'EE_healthy_1.mat','EE_healthy_2.mat','EE_healthy_3.mat','EE_healthy_4.mat','EE_healthy_5.mat','EE_healthy_6.mat','EE_healthy_7.mat','EE_healthy_8.mat','EE_healthy_9.mat','EE_healthy_10.mat','EE_linear_1.mat','EE_linear_2.mat','EE_linear_3.mat','EE_linear_4.mat','EE_linear_5.mat','EE_linear_6.mat','EE_linear_7.mat','EE_linear_8.mat','EE_linear_9.mat','EE_linear_10.mat','EE_MJ_1.mat','EE_MJ_2.mat','EE_MJ_3.mat','EE_MJ_4.mat','EE_MJ_5.mat','EE_MJ_6.mat','EE_MJ_7.mat','EE_MJ_8.mat','EE_MJ_9.mat','EE_MJ_10.mat','IAF_1.mat','IAF_2.mat','IAF_3.mat','IAF_4.mat','IAF_5.mat','IAF_6.mat','IAF_7.mat','IAF_8.mat','IAF_9.mat','IAF_10.mat','OAF_1.mat','OAF_2.mat','OAF_3.mat','OAF_4.mat','OAF_5.mat','OAF_6.mat','OAF_7.mat','OAF_8.mat','OAF_9.mat','OAF_10.mat','OAF_11.mat','OAF_12.mat'};
n_file1=input('How many files? ');

cd(dir_sub)
load('MVC_proc');

for n=1:n_file1 
    
%     for c=1:size(name_conf{n},2)
%      MVCval(c)=MVCb(find(strcmp(name_conf{n}{c},conf_2)==1));
%     end

     
cd(dir_sub);
load(name_file{n});
cd(directory_iniziale);
noChans = 16;

for i=1:noChans-1  
    Data = getfield(EMG,'data');
    emg(:,i)=Data(:,i);
end
clear Activities Data length_sec samplingRate
emg = emg';
[EMG_proc]=preprocessing(emg,3000,MVCb,noChans-1); % MVC and the other EMG are preprocessed in the same way

name=[name_file{n}(1:length(name_file{n})-4),'_p']; 
cd(dir_sub);  
eval(['  save ' name ' EMG_proc MVCb']);
clear emg MVCval EMG_proc i 
end
  