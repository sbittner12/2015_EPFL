subjects = 3:11;
exptype = 'OAF';
trials = 1:12;
ntargets = 8;
nsubjects = length(subjects);
ntrials = length(trials);

for i=1:nsubjects
    subject =subjects(i);
    for j=1:ntrials
        trial = trials(j);
        Event = loadEvent(subject, exptype, trial);
        EEGtriggers = parseTriggers(Event.EEG.Trigger, Event.EEG.End(end));
        EMGtriggers = parseTriggers(Event.EMG.Trigger, Event.EMG.End(end)); 
        trialplot = 1;
        for k=1:(ntargets-1)
            if (Event.EEG.End(2*k) > EEGtriggers(k+1))
                fprintf('Overlap in EEG between targets %d and %d in subject %d trial %d\n', ...
                         k, k+1, subject, trial);
                if (trialplot)
                    plotEvent(Event.EEG);
                    trialplot = 0;
                end
            end
            if (Event.EMG.End(2*k) > EMGtriggers(k+1))
                fprintf('Overlap in EMG between targets %d and %d in subject %d trial %d\n', ...
                         k, k+1, subject, trial);
                if (trialplot)
                    plotEvent(Event.EMG);
                    trialplot = 0;
                end
            end
        end
    end
end
