function EEG = loadEEG(subject, exptype, channels, srate)
%UNTITLED Summary of this function goes here
%   loadEEG(subject, exptype, channels)
%   Detailed explanation goes here
% Setup file I/O
    if isnumeric(subject)
        subject = num2str(subject);
    end
    dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];
    subjectEEGDir = [dataDir, 'EEG', filesep, 'Subject', subject, filesep];

    if (srate == 2048)
        fname =  [subjectEEGDir, 'Preprocessed', exptype, '_2048.mat'];
    elseif (srate == 1000)
        fname =  [subjectEEGDir, 'Preprocessed', exptype, '.mat'];
    else
        fprintf('Erro: invalid srate: must be 1000 or 2048.\n');
        EEG = -1;
        return
    end
    
    if (~(exist(fname, 'file') == 2))
        fprintf('Error: file %s does not exist.\n', fname);
        return;
    else
        load(fname);
    end

    % Determine channel indices
    chanlocs = Trials(1).chanlocs;
    if (~strcmp(channels, 'all'))
        channel_inds = getEEGChannelInds(chanlocs, channels);

        % Only keep data for specified eeg electrodes
        ntrials = length(Trials);
        for i=1:ntrials
            Trials(i).data = Trials(i).data(channel_inds,:);
        end
    end
    EEG = Trials;
end

