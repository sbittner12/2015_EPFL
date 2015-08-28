function emg = extractEMGSignal(subject, exptype, musclestr, trial)
%EXTRACTEMGBLOCK - Get preprocessed EMG signal for a particular reach
%   extractEMGBLOCK(subject, exptype, musclestr, trial)
%   This function is useful for extracting preprocessed EMG data for a
%   particular electrode or group of electrodes and a particular trial
%   Ex:
%      
%   If you choose the 'uniform length' option for the extraction_type
%   parameter, you will have to enter two additional arguements.  The first
%   is the pre move onset buffer in seconds, and the second is the post
%   move onset buffer in seconds.  In the following example, a pre move
%   buffer of 0.7 seconds and a post move buffer of 2.51 seconds are
%   specified.
%   Ex:
%   extractEMGSignal(3,'OAF','DANT',1,1,'uniform length', .70, 2.51);
%
%  PARAMETERS
%    subject - integer or number string which identifies the subject
%      Ex: '3', or 3
%    exptype - string for experiment type based on exoskeleton function
%      Ex: 'Healthy', 'IAF', 'LT', 'MJ', or 'OAF'
%    musclestr - string specifying particular EMG channel
%      Ex: 'DANT', 'BICS', etc.
%    trial - integer specifying trial number
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%
srate_new = 1000;
srate = 3000;
% Setup file I/O
if isnumeric(subject)
    subject = num2str(subject);
end
dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];
subjectEMGDir = [dataDir, 'EMG', filesep, 'Subject', subject, filesep];
subjectEventDir = [dataDir, 'Event', filesep, 'Subject_', subject,filesep]; 

emgstr = getEmgstr(exptype);
eventstr = getEventstr(exptype);

EMGfname =  [subjectEMGDir, emgstr, sprintf('_%d.mat', trial)];
if (~(exist(EMGfname, 'file') == 2))
    fprintf('Error: file %s does not exist.\n', EMGfname);
    return;
else
    load(EMGfname);
end


MVCfname =  [subjectEMGDir, 'MVC_proc_MCC.mat'];
if (~(exist(MVCfname, 'file') == 2))
    fprintf('Error: file %s does not exist.\n', MVCfname);
    return;
else
    load(MVCfname);
end

% Pre-process the EMG data
muscle_ind = getEMGChannelInd(EMG, musclestr);
emg = preprocessingEMG(EMG.data(:,muscle_ind), srate, MVCb(muscle_ind), 1);

end

