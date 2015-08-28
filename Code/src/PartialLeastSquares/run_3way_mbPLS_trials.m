%run_3way_mbPLS_trials

subject = 7;
SUBJECTS = 3:9;
exptype = 'OAF';
musclestr = 'all';
channels = 'all';
trials = 1:12;
targets = 1:8;
ntrials = length(trials);
ntargets = length(targets);
if (subject == 7)
    eegsrate = 2048;
else
    eegsrate = 1000;
end
eegsrate_new = 1000;
eeg_flag = 0;
emg_flag = 0;
correct_labels = {'TRAPS', 'TRAPM', 'DANT', 'DMED', 'DPOS', 'PEC', 'INFRA',...
                  'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};
              
% choose decomposition method
eegMode = 'eegSpectrogram'; %{eegSpectrogram, eegCoherence}
emgMode = 'emgPre'; %{imfAmp, emgAmp)

baseDir = ['E:', filesep, 'Sean', filesep];
dataDir = [baseDir, 'Data', filesep];
resDir = [baseDir, 'Results', filesep, 'PartialLeastSquares', filesep];

eegDir = [dataDir, 'EEG', filesep];
eegSubjectDir = [eegDir, sprintf('Subject%d', subject), filesep];
spectfname = [eegSubjectDir, 'Spectrogram_AllTrials.mat'];
cohfname = [eegSubjectDir, 'Coherence_AllTrials.mat'];

emgDir = [dataDir, 'EMG', filesep];
emgSubjectDir = [emgDir, sprintf('Subject%d', subject), filesep];
emgfname = [emgSubjectDir, 'EMG_AllTrials.mat'];

resfname = [resDir, sprintf('%s_%s_subject%d_AllTrials.mat', eegMode, emgMode, subject)];

[premove_buf, postmove_buf] = getMovementBuffers(SUBJECTS, 'OAF');


switch eegMode
    case 'eegSpectrogram'
        if (exist(spectfname, 'file') == 2)
            fprintf('Loading pre-computed EEG spectrogram.\n');
            load(spectfname);
            eeg_flag = 1;
        else 
            Xdata_group = cell(1, ntrials);
        end
    case 'eegCoherence'
        if (exist(cohfname, 'file') == 2)
            fprintf('Loading pre-computed EEG coherence.\n');
            load(cohfname);
            eeg_flag = 1;
        else 
            Xdata_group = cell(1, ntrials);
        end
end


if (exist(emgfname, 'file') == 2)
    fprintf('Loading pre-computed EMG.\n');
    load(emgfname);
    emg_flag = 1;
else
    Ydata_group = cell(1, ntrials);
end

EEG = loadEEG(subject, exptype, 'all', eegsrate);

% find epoch sample length
Event = loadEvent(9, exptype, 1);
emgdata = extractEMGSignal(subject, exptype, 'TRAPS', 1);
emg = extractEMGSignalFast(emgdata, Event, 2, 'uniform length', premove_buf, postmove_buf);
epochsamples = size(emg, 1);

if (~emg_flag || ~eeg_flag)
    for i=1:ntrials
        trial = trials(i);
        fprintf('trial %d\n', i);
        Event = loadEvent(subject, exptype, trial);
        eegdata = EEG(trial).data';
        emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
        srate_new = EEG(trial).srate;
        eegData = [];
        emgData = [];
        for j=1:ntargets
            %tic;
            target = targets(j);
            eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, ...
                                   'uniform length', premove_buf, postmove_buf);
            emg = extractEMGSignalFast(emgdata, Event, target, ...
                                   'uniform length', premove_buf, postmove_buf);
            eegData = [eegData, eeg'];
            emgData = [emgData, emg'];

        end
        
        if (subject == 6)
            eegData = resample(double(eegData'), 1000,2048)';
            eegsrate = 1000;
        end

        if (~emg_flag)
            switch emgMode % extract temporal features of EMG data
                case 'imfAmp'  % extract amplitudes of EMD components
                    imfData = [];
                    for ch = 1:size(emgData, 1)
                        imf = emd2(emgData(ch, :), 'maxmodes', 5, 'display', 0);
                        imfData = [imfData; imf(1:3, :)];
                    end

                    imfAmp = zeros(size(imfData));
                    for ch = 1:size(imfData)
                        imfAmp(ch, :) = est_amp( imfData(ch, :), 2, 60 );
                    end
                    emgComp = imfAmp;

                case 'emgAmp' % extract amplitude of raw EMG data
                    fprintf('Calculating EMG amplitude...  ');
                    emgAmp = zeros(size(emgData));
                    for ch = 1:size(emgData)
                        emgAmp(ch, :) = est_amp( emgData(ch, :), 2, 60 );
                    end      
                    emgComp = emgAmp;
                    fprintf('done.\n');
                case 'emgPre'
                    emgComp = emgData;
            end   
        end

        if (~eeg_flag)
            switch eegMode % extract temporal features of EEG data
                case 'eegSpectrogram' % spectrogram of preprocessed EEG data
                    fprintf('Calculating EEG spectrogram...  \n');
                    eegSData = [];
                    WL = eegsrate*1.5;  % Could use differnt value based on ricampionas
                    for ch = 1:size(eegData)
                        tic;
                        fprintf('Channel %d\n', ch);
                        [S,F,Tspect] = spectrogram(double(eegData(ch, :)), WL, ...
                            WL-20, 3:0.5:50, eegsrate, 'yaxis');

                        eegSData(:, :, ch) = (abs((S).'));
                        fprintf('Channel spectrogram took %2.2f seconds.\n', toc);
                    end

                case 'eegCoherence' % coherence between EEG channels
                    fprintf('Calculating EEG coherence...  ');
                    WL = 500;
                    T = [];
                    curpos = 0; ll = 0;
                    while curpos+WL < size(eegData, 2)
                        T(ll+1)= curpos/eeg_sample_rate;
                        curpos = curpos + shift_len;
                        ll = ll +1; 
                    end
                    eegCoh = NaN(ll, 100, nchoosek(size(eegData, 1), 2));
                    counter = 1; 
                    for m = 1:size(eegData, 1)
                        for n = m+1:size(eegData, 1)   
                            for z = 1:ll
                                index = (z-1)*shift_len + 1 :  ((z-1)*shift_len + WL);
                                [Cxy,F] = mscohere(eegData(m, index),eegData(n, index), ...
                                    hamming(100), 100-shift_len, ...
                                    eeg_sample_rate*2,eeg_sample_rate); 
                                eegCoh(z, :, counter) = Cxy((F>0)& (F<=50));
                                F = F((F>0)& (F<=50));
                            end
                            counter = counter+1;
                        end
                    end
                    fprintf('done.\n');
            end
        end

        if (~eeg_flag)
            switch eegMode
                case 'eegSpectrogram'
                    Xpls=nprocess(eegSData, [1 0 0], [1 0 0], [], [], 1, 1);
                case 'eegCoherence'
                    Xpls=nprocess(eegCoh, [1 0 0], [1 0 0], [], [], 1, 1);
            end
            xpls_samples = size(Xpls, 1);
            fprintf('writing to data group %d\n', i);
            Xdata_group{i} = Xpls;
        else
            xpls_samples = size(Xdata_group{i},1);
        end

        if (~emg_flag)
            Ypls=nprocess(resample(emgComp.', xpls_samples, size(emgComp,2)), [1 0], [1, 0], [], [], 1, 1);
            Ydata_group{i} = Ypls;
        end
        fprintf('done.\n');
    end
    

    if (~eeg_flag)
        fprintf('saving Xdata_group\n');
        save(spectfname, 'Xdata_group', 'F');
    end
    if (~emg_flag)
        fprintf('saving Ydata_group\n');
        save(emgfname, 'Ydata_group');
    end
end



% perform MBPLS decomposition
numFactor = 6;
fprintf('running mbpls\n');
if (exist(resfname, 'file') == 2)
    load(resfname);
else
    [Tb,W1b,W2b,Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, W1b_temp, W2b_temp, Tu_orig] = ...
        MBtriPLS2_cmmnW1b(Xdata_group, Ydata_group, numFactor, size(Xdata_group{1}, 1)*1e-12);
    save(resfname, 'Tb', 'W1b', 'W2b', 'Wt', 'Tt', 'Ub', 'Qb', 'Wu', 'Tu',...
                   'ssx', 'ssy', 'W1b_temp', 'W2b_temp', 'Tu_orig');
end


% plot componenets
if (strcmp(musclestr, 'all'))
    for f=1:numFactor
        figure; 

        subplot(3, 2, 1);
        plot(F, abs(W1b{1}(:, f))); title('Spectral Signature');
        xlim([min(F) max(F)]); 
        xlabel('Frequency (Hz)'); 

        subplot(3, 2, 3);
        nsamples = size(Tt,1);
        plot(1:nsamples, unitVariance( Tt(:, f).'), 'r', ... 
             1:nsamples, unitVariance(Tu(:, f).'), 'b'); %, ...
             %1:nsamples, unitVariance(Tu_orig(:,f).'), 'g'); 
        legend('EEG', 'EMG');%, 'original EMG activation');
        title('Temporal Signature');
        xlim([1 nsamples]); 

        subplot(3, 2, 5);
        %topoplot(W2b{1}(:,f)+.9, EEG(1).chanlocs, 'style', 'blank');
        topoplot(W2b{1}(:,f), EEG(1).chanlocs);
        suptitle(sprintf('Component %i', f));

        subplot(3,2, 2), stem(Wu(:,f)), title(sprintf('Response super block weights', f));

        [Qb_dotprods, Qb_colors] = rankQb(Qb, f);
        [~, ind] = max(Qb_dotprods);
        subplot(3,2,4), bar(Qb{ind}(:,f), 'FaceColor', Qb_colors{ind}), title('Best fit synergy');
        ax = gca;
        set(ax, 'XTickLabel', correct_labels, 'FontSize', 6);
        subplot(3,2,6), bar(Qb_dotprods), title('Dot product with best fit synergy');
        xlabel('trial');
    end
else
    for f = 1:numFactor
        figure; 

        subplot(3, 2, [1 2]);
        plot(F, abs(W1b{1}(:, f))); title('Spectral Signature');
        xlim([min(F) max(F)]); 
        xlabel('Frequency (Hz)'); 

        subplot(3, 2, [3 4]);
        nsamples = size(Tt,1);
        plot(1:nsamples, unitVariance( Tt(:, f).'), 'r', ... 
             1:nsamples, unitVariance(Tu(:, f).'), 'b', ...
             1:nsamples, unitVariance(Tu_orig(:,f).'), 'g'); 
        legend('EEG', 'EMG', 'original EMG activation');
        title('Temporal Signature');
        xlim([1 nsamples]); 
        xlabel('sample index'); 

        subplot(3, 2, [5,6]);
        topoplot(W2b{1}(:,f), EEG(1).chanlocs);
        suptitle(sprintf('Component %i', f));

    end
end
