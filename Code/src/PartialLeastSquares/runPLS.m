function [X,Y] = runPLS(subject, exptype, channels, trials, targets)
%RUNPLS Summary of this function goes here
%   Detailed explanation goes here

T = .001; % 1 kHz
Fs = 1/T;
windowsize = 300;
noverlap = 285;
musclestr = 'all';
fftbins = (windowsize/2)+1;
freq_res = Fs/windowsize;
if isnumeric(subject)
    subject = num2str(subject);
end
ntrials = length(trials);
ntargets = length(targets);
                 
eventstr = getEventstr(exptype);

EEG = loadEEG(subject, exptype, channels);

eeg_all = [];
emg_all = [];
for i=1:ntrials
    trial = trials(i);
    fprintf('Trial %d\n', trial);
    Event = loadEvent(subject, exptype, trial);
    eegdata = EEG(trial).data';
    emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
    srate_new = EEG(trial).srate;
    for j=1:ntargets
        %tic;
        target = targets(j);
        % Extract corresponding eeg and emg signals
        extraction_type = 'reach duration';
        eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, extraction_type);
        emg = extractEMGSignalFast(emgdata, Event, target, extraction_type);
        lengthdiff = size(eeg,1) - size(emg,1);
        if (lengthdiff)
            eeg = eeg(1:end-lengthdiff, :);
        end
        assert(size(eeg,1) == size(emg,1));
        eeg_all = [eeg_all; eeg];
        emg_all = [emg_all; emg];
    end
end
size(eeg_all)
size(emg_all)
[N, Ne] = size(emg_all);
Nt = ceil((N-windowsize)/(windowsize-noverlap));
Y = zeros(Nt, Ne);
startind = 1;
stopind = windowsize;
di = windowsize - noverlap;
i = 1;
while(startind < N - windowsize)
    for j=1:Ne
        Y(i,j) = rms(emg_all(startind:stopind, j));
    end
    startind = startind + di;
    stopind = stopind + di;
    i = i + 1;
end   

X = spectrogram(double(eeg_all), windowsize, noverlap)';
X = X(:,1:26); % Keep only bands up to 50 Hz

ncomps = 10;
[P,Q,T,U] = pls_fromscratch(X,Y,ncomps);


end

