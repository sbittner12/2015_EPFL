% Generate new .mat and .bdf files for each individual target reach of each
% trial resampled to the minimum epoch sample length.

subjects = [3,4];
trials = 1:12;
exptype = 'OAF';
channels = 'all';
musclestr = 'all';

nsubjects = length(subjects);
ntrials = length(trials);
ntargets = 8;
extraction_type = 'reach duration';

baseDir = ['E:', filesep, 'Sean', filesep, 'Data', filesep];
eegDir = [baseDir, 'EEG', filesep];
emgDir = [baseDir, 'EMG', filesep];

srate = 2048;
srate_emg = 1000;

min_epoch_samples = getShortestEpoch()

for i=1:nsubjects
    subject = subjects(i);
    fprintf('subject %d\n', subject);
    EEGsubject = loadEEG(subject, exptype, channels);
    
    eegSubjectEpochDir = [eegDir, sprintf('Subject%d', subject), ...
                                  filesep, 'SingleEpochMat', filesep];
    if (exist(eegSubjectEpochDir, 'dir') ~= 7)
        mkdir(eegSubjectEpochDir);
    end
    
    emgSubjectEpochDir = [emgDir, sprintf('Subject%d', subject), ...
                                  filesep, 'SingleEpochMat', filesep];
    if (exist(emgSubjectEpochDir, 'dir') ~= 7)
        mkdir(emgSubjectEpochDir);
    end                          
    
        
    for j=1:ntrials
        trial = trials(j);
        fprintf('trial %d\n', trial);
        Event = loadEvent(subject, exptype, trial);
        eegdata = EEGsubject(trial).data';
        emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
        eegtimes = EEGsubject(trial).times';
        srate_new = EEGsubject(trial).srate;
        for target=1:ntargets
            fprintf('target %d\n', target);
            eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, extraction_type);
            emg = extractEMGSignalFast(emgdata, Event, target, extraction_type);
            times = extractEEGSignalFast(eegtimes, Event, srate_new, target, extraction_type);

            % Setup EEG file
            [epoch_samples, numeeg] = size(eeg);
            eeg_resampled = zeros(numeeg, min_epoch_samples);
            times_resampled = zeros(1, min_epoch_samples);
            for i=1:numeeg
                eeg_resampled(i,:) = ricampiona(eeg(:,i), min_epoch_samples);
            end
            times_resampled(1,:) = ricampiona(times(:,1), min_epoch_samples);
            EEG = EEGsubject(trial);
            EEG.data = eeg_resampled;
            EEG.times = times_resampled;
            EEG.srate = min_epoch_samples; % Makes it simpler for storing bdf files
            EEG.pnts = min_epoch_samples;
            EEG.setname = sprintf('Sub%d_OAF%d_%d_resampled_to_min_epoch', subject, trial, target);
            EEG.event = [];
            
            % Setup EMG file
            [epoch_samples, numemg] = size(emg);
            emg_resampled = zeros(numemg, min_epoch_samples);
            for i=1:numemg
                emg_resampled(i,:) = ricampiona(emg(:,i), min_epoch_samples);
            end
            EMG.data = emg_resampled;
            EMG.framerate = srate_emg * min_epoch_samples / epoch_samples;            
            
            % Write EEG file
            eegfnameprefix = [eegSubjectEpochDir, sprintf('OAF_%d_%d', trial, target)];
            save([eegfnameprefix, '.mat'], 'EEG');
            
            % Write EMG file
            emgfnameprefix = [emgSubjectEpochDir, sprintf('OAF_%d_%d', trial, target)];
            save([emgfnameprefix, '.mat'], 'EMG');
        end
            
    end
end