function [raw, preprocessed] = preprocessView(datatype, subject, exptype, electrode, trial, targets)
%PREPROCESSVIEW Summary of this function goes here
%   Detailed explanation goes here
    SUBJECTS = [3,4,7];
    extraction_type = 'uniform length';
    % Setup file I/O
    baseDir = ['E:', filesep, 'Sean', filesep];
    dataDir = [baseDir, 'Data', filesep];
    eventDir = [dataDir, 'Event', filesep];
    eventSubjectDir = [eventDir, 'Subject_' subject, filesep];
    switch exptype
        case 'Healthy'
            emgstr = 'EE_healthy';
            eventstr = 'HM';
        case 'IAF'
            emgstr = 'IAF';
            eventstr = 'IAF';
        case 'LT'
            emgstr = 'EE_linear';
            eventstr = 'LT';
        case 'MJ'
            emgstr = 'EE_MJ';
            eventstr = 'MJ';
        case 'OAF'
            emgstr = 'OAF';
            eventstr = 'OAF';
        otherwise
            fprintf('Error: invalid exptype.\n');
    end
    eventfname = [eventSubjectDir, sprintf('%s_%d_event', eventstr, trial)];
    load(eventfname);
    [premove_buf, postmove_buf] = getMovementBuffers(SUBJECTS, eventstr);
    [target_begin_ind, target_end_ind] = parseTargets(targets, Event);
    
    % ********************* EEG ********************************
    if (strcmp(datatype, 'EEG'))
        % Load raw data
        srate = 2048;
        
        eegDir = [dataDir, 'EEG', filesep];
        eegSubjectDir = [eegDir, 'Subject', subject, filesep];
        if (strcmp(subject, '3') || strcmp(subject, '4'))
            eegfname = [eegSubjectDir, sprintf('%s%d.bdf', exptype, trial)];
        else
            eegfname = [eegSubjectDir, sprintf('%s_%d.bdf', exptype, trial)];
        end
        if (exist(eegfname, 'file') == 2)
            eeglab;
            data = pop_biosig(eegfname);
        else
            fprintf('Error: File %s not found.\n', eegfname);
        end
        nsamples = size(data.data,2);
        % Determine index of electrode
        ch_file=[eegDir, filesep, 'channels_location(64).ced']; % to change 
        chlocs=readlocs(ch_file);
        nchanlocs = length(chlocs);
        for j=1:nchanlocs
            %fprintf('j: %d, %s\n', j, chlocs(j).labels);
            if (strcmp(electrode, chlocs(j).labels))
                channel_ind = j;
                break;
            elseif (j == nchanlocs)
              fprintf('Error: channel location %s not found in data.\n', electrode);
              return;
            end
        end
        % Determine temporal index of data
        ntargets = length(targets);
        if (ntargets == 1);
            move_onset = Event.EEG.Start(2*target_begin_ind-1);
            move_onset_t = move_onset/srate;
            EEG_begin_ind = round(move_onset - (premove_buf*srate));
            EEG_end_ind = round(move_onset + (postmove_buf*srate));
            begin_t = move_onset_t - premove_buf;
            end_t = move_onset_t + postmove_buf;
        else
            EEG_begin_ind = 1;
            EEG_end_ind = nsamples;
            begin_t = EEG_begin_ind/srate;
            end_t = EEG_end_ind/srate;
        end
        
        % Extract raw and preprocessed signals
        raw = data.data(channel_ind, EEG_begin_ind:EEG_end_ind);
        preprocessed = extractEEGSignal(subject, exptype, electrode, trial, targets, ...
                       extraction_type, premove_buf, postmove_buf);
        
    % ********************* EMG ********************************   
    elseif (strcmp(datatype, 'EMG'))
        srate = 3000;
        subjectEMGDir = [dataDir, 'EMG', filesep, 'Subject', subject, filesep];
        EMGfname =  [subjectEMGDir, emgstr, sprintf('_%d.mat', trial)];
        if (~(exist(EMGfname, 'file') == 2))
            fprintf('Error: file %s does not exist.\n', EMGfname);
            return;
        else
            load(EMGfname);
        end
        
        % Determine muscle indices
        nchanlocs = length(EMG.label);
        for i=1:nchanlocs
            if (strcmp(electrode, EMG.label{i}))
                muscle_ind = i;
                break;
            elseif (i == nchanlocs)
                fprintf('Error: could not find muscle %s.\n', musclestr);
            end
        end
        % Determine temporal index of data
        move_onset = Event.EMG.Start(2*target_begin_ind-1);
        move_onset_t = move_onset/srate;
        EMG_begin_ind = round(move_onset - (premove_buf*srate));
        EMG_end_ind = round(move_onset + (postmove_buf*srate));
        begin_t = move_onset_t - premove_buf;
        end_t = move_onset_t + postmove_buf;
        
        raw = EMG.data(EMG_begin_ind:EMG_end_ind, muscle_ind);
        preprocessed = extractEMGSignal(subject, exptype, electrode, trial, targets, ...
                       extraction_type, premove_buf, postmove_buf);
    else
        fprintf('Error: datatype parameter must be EEG or EMG.\n');
    end
    
    % Plot Raw
    nsamples = length(raw);
    figure, subplot(2,1,1), plot(linspace(begin_t, end_t, nsamples), raw);
    if (length(targets) == 1)
        title(sprintf('%s Raw: subject %s trial %d target %d %s', datatype, ...
                       subject, trial, targets, electrode));
    else
        title(sprintf('%s Raw: subject %s trial %d target %d:%d %s', datatype, ...
                       subject, trial, targets(1), targets(end), electrode));
    end
    xlabel('time (s)');
    ylabel('Amplitude (V)');
    minraw = min(raw);
    maxraw = max(raw);
    spacing = abs(maxraw - minraw)/10;
    axis([begin_t, end_t, minraw-spacing, maxraw+spacing]);

    % Plot preprocessed
    nsamples = length(preprocessed);
    subplot(2,1,2), plot(linspace(begin_t, end_t, nsamples), preprocessed), title('PREPROCESSED');
    if (length(targets) == 1)
        title(sprintf('%s Preproc: subject %s trial %d target %d %s', datatype, ...
                       subject, trial, targets, electrode));
    else
        title(sprintf('%s Preproc: subject %s trial %d target %d:%d %s', datatype, ...
                       subject, trial, targets(1), targets(end), electrode));
    end
    xlabel('time (s)');
    ylabel('Amplitude (V)');minpre = min(preprocessed);
    maxpre = max(preprocessed);
    spacing = abs(maxpre - minpre)/10;
    axis([begin_t, end_t, minpre-spacing, maxpre+spacing]);

end

