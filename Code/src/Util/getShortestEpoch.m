function [min_epoch_eeg] = getShortestEpoch()
%% finds shortest epoch time from EEG
subjects = 3:10;
trials = 1:12;
exptype = 'OAF';

nsubjects = length(subjects);
ntrials = length(trials);
ntargets = 8;

eeg_srate_new = 1000;
eeg_srate = 2048;
min_epoch_matrix = zeros(nsubjects, ntrials);

for i=1:nsubjects
    subject = subjects(i);
    for j=1:ntrials
        trial = trials(j);
        Event = loadEvent(subject, exptype, trial);
        if (isempty(Event))
            continue
        end
        triggers = parseTriggers(Event.EEG.Trigger, Event.EEG.End(2*ntargets));
        epoch_lengths = round((Event.EEG.End(2:2:(2*ntargets)) - triggers) * eeg_srate_new / eeg_srate);
        min_epoch_trial = min(epoch_lengths);
        min_epoch_matrix(i,j) = min_epoch_trial;
        if ((i==1 && j ==1) || (min_epoch_trial < min_epoch_eeg))
            min_epoch_eeg = min_epoch_trial;
        end
    end
end

%% finds shortest epoch time from EMG
exptype = 'OAF';

nsubjects = length(subjects);
ntrials = length(trials);
ntargets = 8;

emg_srate_new = 1000;
emg_srate = 3000;
min_epoch_matrix = zeros(nsubjects, ntrials);

for i=1:nsubjects
    subject = subjects(i);
    for j=1:ntrials
        trial = trials(j);
        Event = loadEvent(subject, exptype, trial);
        if (isempty(Event))
            fprintf('Warning: Event file for subject %d and trial %d does not exist.\n', subject, trial);
            continue
        end
        triggers = parseTriggers(Event.EMG.Trigger, Event.EMG.End(2*ntargets));
        epoch_lengths = round((Event.EMG.End(2:2:(2*ntargets)) - triggers) * emg_srate_new / emg_srate);
        min_epoch_trial = min(epoch_lengths);
        min_epoch_matrix(i,j) = min_epoch_trial;
        if ((i==1 && j ==1) || (min_epoch_trial < min_epoch_emg))
            min_epoch_emg = min_epoch_trial;
        end
    end
end

if (min_epoch_eeg == min_epoch_emg)
    min_epoch = min_epoch_eeg;
else
    fprintf('Error: minimum epoch lengths extracted from EEG and EMG indices do not match.\n');
    min_epoch = -1;
end